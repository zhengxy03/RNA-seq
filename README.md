# RNA-seq workflow
# 1 mkdir
```mkdir biosoft
mkdir -p project/rat
cd project/rat
mkdir annotation genome sequence output script
```
# 2 install tools
cd ~/biosoft
basic workfolw:
* download
* decompress
* change dir
* config
* make
* export PATH
## 2.1 sra
download and convert data
```cd ~/biosoft
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-win64.zip
unzip -d sratoolkit.2.11.3-win64.zip
mv sratoolkit.2.11.3-win64 sratoolkit.2.11.3
cd sratoolkit.2.11.3/bin

#set PATH
export PATH="$(pwd):$PATH"

#or:brew install sratoolkit
```
## 2.2 fastqc
```
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
cd FastQC
export PATH="$(pwd):$PATH"

#brew install fastqc
```
## 2.3 multiqc
```
pip install multiqc
```
## 2.4 cutadapt
```
pip install cutadapt
```
## 2.5 quality control and trim adapter
* Trim Galore
```
wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.3.tar.gz -O TrimGalore.gz
gzip -d TrimGalore.gz
```
* fastp
* trimmomatic
```
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

cd Trimmomatic-0.39
export PATH="$(pwd):$PATH"
```
## 2.6 hisat2
```
wget http://www.di.fc.ul.pt/~afalcao/hisat2.1/hisat2.1_Windows.zip
unzip hisat2.1_Windows.zip

cd hisat2.1
export PATH="$(pwd):$PATH"
source ~/.bashrc
hisat2 -h
```
## 2.7 sortmerna[optional]
```
wget https://github.com/sortmerna/sortmerna/releases/download/v4.3.7/sortmerna-4.3.7-Linux.tar.gz -O sortmerna-4.3.7.tar.gz
tar -xvzf sortmerna-4.3.7.tar.gz
cd sortmerna-4.3.7

echo 'export PATH=~/biosoft/sortmerna-4.3.7:$PATH' >> ~/.bashrc
source ~/.bashrc

sortmerna -h

wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
tar -xvzf database.tar.gz
mkdir -p ~/database/sortmerna_db/rRNA_databases
mv ./*.fasta ~/database/sortmerna_db/rRNA_databases

cd ~/database/sortmerna_db
#sortmerna_ref_data=$(pwd)/rRNA_databases/smr_v4.3_default_db.fasta


```
## 2.8 samtools
```
wget "https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2"
tar jxvf samtools-1.21.tar.bz2
cd samtools-1.21
./configure --prefix=$(pwd)
make -j 4
export PATH="$(pwd):PATH"
```
## 2.9 HTseq
```
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple HTseq
```
## 2.10 parallel
```
brew install parallel
```
## 2.11 StringTie[optional]
```
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.3.tar.gz
```
## 2.12 Ballgown[optional]
```
wget https://bioconductor.org/packages/release/bioc/bin/windows/contrib/4.4/ballgown_2.36.0.zip
```
# 3 data download
## 3.1 ref data
### 3.1.1 genomes data(.fasta)
* download from Ensembl
```
wget https://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
```
* decompress
* rename
```
mv Rattus_norvegicus.Rnor_7.2.dna.toplevel.fa rn7.raw.fa
```
* check(chr...)
* clear some additional text
```
cat rn7.raw.fa | perl -ne 'if(m/^>(.+?)(?:\s|$)/){ print ">$1\n";}else{print}' > rn6.fa
rm rn7.raw.fa
```
* count per chr lengths
```
cat rn7.fa | perl -n -e '
  s/\r?\n//;
  if(m/^>(.+?)\s*$/){
    $title = $1;
    push @t, $title;
  }elsif(defined $title){
    $title_len{$title} += length($_);
  }
  END{
    for my $title (@t){
      print "$title","\t","$title_len{$title}","\n";
    }
  }
'
```
chr1:
```
cat rn7.fa | perl -n -e '
  if(m/^>/){
    if(m/>1$/){
      $title = 1;
    }else{
      $title = 0;
    }
  }else{
    push @s, $_ if $title;
  }
  END{
    printf ">1\n%s", join("", @s);
  }
' > rn7.chr1.fa
```
### 3.1.2 Genome Index File[optional]
* hisat2
```
wget https://genome-idx.s3.amazonaws.com/hisat/rn6_genome.tar.gz
gzip -d rn6_genome.tar.gz
```
### 3.1.3 annotation(.gff)
```
cd ~/project/rat/annotation
wget https://ftp.ensembl.org/pub/release-112/gff3/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.112.gff3.gz
gzip -d Rattus_norvegicus.mRatBN7.2.112.gff3.gz
mv Rattus_norvegicus.mRatBN7.2.112.gff3 rn7.gff
head rn7.gff
```

## 3.2 Experimental Data(.sra)
cd ~/project/rat/sequence
* NCBI-GEO >> GEO accession
* SRA Run Selector
* download
```
nohup prefetch SRR2190795 SRR224018{2..7} SRR2240228 -O . &
```
* convert format(.sar > .fastq > .gz)
```
array=(SRR2190795 SRR224018{2..7} SRR2240228)
for i in "${array[@]}"; do
  dir="$HOME/project/rat/sequence/$i"
  cd "${dir}"
  mv ${dir}/* $HOME/project/rat/sequence/srr
done

parallel -j 4 "
  fastq-dump --split-3 --gzip {1}
" ::: $(ls *.sra)

rm *.sra

gzip -d -c SRR2190795.fastq.gz | head -n 20
```
# 4 quality control(.fastq)
cd ~/project/rat/output
## 4.1 quality assessment(fastqc)
input:sequence (.fastq.gz)
output:fastqc (fastqc.html;fastqc.zip)
* fastqc [选项][测序文件]
```
cd ~/project/rat/sequence
mkdir ../output/fastqc
# -o 指定输出文件夹
# *.gz 表示这个目录下以 .gz 的所有文件
fastqc -t 6 -o ../output/fastqc *gz

cd ~/project/rat/output/fastqc
ls
```
* multiqc
  * per seq GC content
```
multiqc .
```
解读：https://www.jianshu.com/p/db2100baabb5
  * seq quality histograms(phred score < 30)
  * count numbers of per seq quality scores
  * adapter content
## 4.2 cut adapter and low quality bases(trimmomatic)
input:sequence (.fastq.gz)
output:adapter (.fastq.gz)
```
mkdir -p ../output/adapter

cd ~/project/rat/sequence
for i in $(ls *.fastq.gz); do
  cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
  --minimum-length 30 --overlap 4 --trim-n \
  -o ../output/adapter/${i} ${i}
done
```
## 4.3 trim low quality regions(trimmomatic)
input adapter (.fastq.gz)
output:trim (.fastq.gz)
```
mkdir ../trim

parallel -j 4 "
  java -jar ~/biosoft/Trimmomatic-0.39/trimmomatic-0.39.jar \
  SE -phred33 {1} ../trim/{1} \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:5:15 MINLEN:30 \
" ::: $(ls *.gz)
```
## 4.4 quality assessment again(fastqc)
input:trim (.fastq.gz)
output:fastqc_trim (fastqc.html;fastqc.zip)
```
mkdir ~/project/rat/output/fastqc_trim
parallel -j 4 "
  fastqc -t 4 -o ../fastqc_trim {1}
" ::: $(ls *.gz)

cd ../fastqc_trim
multiqc .
```
# 5 rm rRNA[optional](.fastq)
cd ~/project/rat/output/rRNA
input:trim (.fastq.gz > .fastq)
```
mkdir -p ~/project/rat/output/rRNA/discard

cd trim
parallel -j 4 "
  gzip -d {1}*.fq.gz

  # euk_rNRA_ref_data就是之前安装sortmerna的时候定义的数据库文件
  # --reads  : 测序文件
  # --aligned: 与rRNA数据库能比对上的序列(后续不需要的)
  # --other  : 与rRNA数据库不能比对上的序列(后续需要的)
  # --fastx  : 输出fastq文件
  # --log    : 生成日志文件
  # -a       : 线程数
  # -v       : 吵闹模式
  
  # 注意--aligned和--other后接文件名前缀，不用在加什么 .fq 或者 .fastq之类的，否则将会生成 xxx.fq.fq
  sortmerna \
  sortmerna \
    --ref $euk_rNRA_ref_data \
    --reads {1}.fastq \
    --aligned ../rRNA/discard/{1} \
    --other ../rRNA/{1} \
    --fastx \
    --log \
    -a 4 \
    -v 

    gzip ../rRNA/discard/{1}.fastq
    gzip ../rRNA/{1}.fastq
" ::: $(ls *.fastq.gz | perl -n -e 'print $1. "\n" if m/(.+?)_/')
```
# 6 seq alignment
hisat2 samtools
## 6.1 build index(.ht2)
input:genome (.fa)<br>
output:index (.ht2)
```
hisat2-build [选项] [基因组序列(.fa)] [索引文件的前缀名]

cd ~/project/rat/genome
mkdir index
cd index

hisat2-build -p 6 ../rn7.chr1.fa rn7.chr1
```
## 6.2 alignment(.sam)
input:trim/rRNA (.fastq.gz) -- index()<br>
output:align (.sam;.log)
```
hisat2 [选项] -x [索引文件] [ -1 1测序文件 -2 2测序文件 -U 未成对测序文件 ] [ -S 输出的sam文件 ]

cd ~/project/rat/output
mkdir align
cd trim
parallel -k -j 4 "
  hisat2 -t -x ../../genome/index/rn7.chr1 \
  -U {1}.fastq.gz -S ../align/{1}.sam \
  2>../align/{1}.log
" ::: $(ls *.gz | perl -p -e 's/.fastq.gz$//')

cd ../align
cat SRR2190795.log
```
* summarize alignment rate and time
input:align (.log)
```
file_list=$(ls *.log)
for i in ${file_list[@]};
do
  prefix=$( echo ${i} | perl -p -e "s/\.log//" )
  echo -n -e "${prefix}\t"
  perfix = $(echo ${i} | perl -p -e 's/\.log//')
  echo -n -e '${perfix}\t'
  cat ${i} |
    grep -E "(overall alignment rate)|(Overall time)" |
    perl -n -e '
      if(m/alignment/){
        $hash{precent}=$1 if m/([\d.]+)%/;
      }elsif(m/time/){
        if(m/(\d\d):(\d\d):(\d\d)/){
          my $time = $1 * 60 + $2 + $3 / 60;
          $hash{time} = $time;
        }
      }
      END{
        $hash{precent} = "NA" if not exists $hash{precent};
        $hash{time} = "NA" if not exists $hash{time};
        printf "%.2f\t%.2f\n", $hash{precent}, $hash{time};
      }
    '
done
```
## 6.3 convert format&sort(.sam->.bam)
samtools
input:align (.sam)
```
parallel -k -j 4 "
  samtools sort -@ 2 {1}.sam > {1}.sort.bam && samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

rm *.sam
ls
```
# 7 expression level count
> HTseq-count:determine if the RNA reads belong to one gene(three models:union,intersection_strict,intersection_nonempty)

input:align (.bam);annotation (.gff) <br>
output:HTseq (.count)(id + count)
```
htseq-count [options] <alignment_files> <gff_file>

cd ~/project/rat/output
mkdir HTseq

cd align
parallel -j -4 "
  htseq-count -s no -f bam {1}.sort.bam ../../annotation/rn7.gtf \
    >../HTseq/{1}.count 2>../HTseq/{1}.log
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

cd ../HTseq
cat SRR2190795.count | head -n 10
```
# 8 merge&standardization
```
rm(list=ls())

setwd=("//wsl.localhost/Ubuntu/home/zxy0303/project/rat/output/HTseq")
#perfix-file
files <- list.files(".", "*.count")
f_lists <- list()
for (i in files){
  perfix = gsub("(_\\w+)?\\.count", "", i ,perl=TRUE)
  f_lists[[perfix]]= i
}

id_list <- names(f_lists)
count <- 0
data <- list()
for (i in id_list){
  count <- count +1
  a <- read.table(f_lists[[i]], sep="\t", col.names=c("gene_id", i))
  data[[count]] <- a
}

data_merge <- data[[1]]
for (i in seq(2,length(id_list))){
  data_merge <- merge(data_merge, data[[i]], by="gene_id)>
}

write.csv(data_merge, "merge.csv", quote= FALSE, row.names = FALSE)
```

