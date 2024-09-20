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
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip

cd Trimmomatic-0.38
export PATH="$(pwd):$PATH"
```
## 2.6 hisat2
```
wget http://www.di.fc.ul.pt/~afalcao/hisat2.1/hisat2.1_Windows.zip
unzip hisat2.1_Windows.zip

cd hisat2.1
export PATH="$(pwd):$PATH"
```
## 2.7 sortmerna[optional]
```
wget https://github.com/sortmerna/sortmerna/releases/download/v4.3.7/sortmerna-4.3.7-win64.7z -O sortmerna-4.3.7.7z
7z x sortmerna-4.3.7.7z
cd sortmerna-4.3.7

./configure --prefix=$PWD
make -j 4
./sortmerna --help
export PATH="$(pwd):$PATH"
mv ./rRNA_databases/ ~/database/sortmerna_db/rRNA_databases
```
## 2.8 samtools
```
wget -c https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
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
mv Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa rn6.raw.fa
```
* check(chr...)
* clear some additional text
```
cat rn6.raw.fa | perl -ne 'if(m/^>(.+?)(?:\s|$)/){ print ">$1\n";}else{print}' > rn6.fa
rm rn6.raw.fa
```
* count per chr lengths
```
cat rn6.fa | perl -n -e '
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
cat rn6.fa | perl -n -e '
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
' > rn6.chr1.fa
```
### 3.1.2 Genome Index File[optional]
* hisat2
### 3.1.3 annotation(.gff)

## 3.2 Experimental Data(.sra)
* NCBI-GEO >> GEO accession
* SRA Run Selector
* download
```
nohup prefetch SRR2190795 SRR224018{2..7} SRR2240228 -o . &
```
* convert format(.sar > .fastq > .gz)
```
parallel -j 4 "
    fastq-dump --split-3 --gzip {1}
" ::: $(ls *.sra)
rm *.sra
```
# 4 quality control
cd ~/project/rat/output
## 4.1 quality assessment
* fastqc [选项][测序文件]
```
mkdir -p ../output/fastqc
fastqc -t 6 -o ../output/fastqc *.gz

cd ../output/fastqc
ls
```
* multiqc
    ** per seq GC content
    ```
    cd ../output/fastqc
    multiqc .
    ```
    ** seq quality histograms(phred score < 30)
    ** count numbers of per seq quality scores
    ** adapter content
## 4.2 cut adapter and low quality bases
```
mkdir -p ../output/adapter
for i in $(ls *.fastq.gz);
do
    cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT --minimum-length 30 --overlap 4 --trim-n \
    -o ../output/adapter/${i} ${i}
done
```
## 4.3 trim low quality regions again
```
cd  ~/project/rat/output/adapter/
mkdir ../trim

parallel -j 4 "
    java -jar ~/biosoft/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 {1} ../trim/{1} \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:5:15 MINLEN:30 \
" ::: $(ls *.gz)
```
## 4.4 quality assessment again
```
mkdir ~/project/rat/output/fastqc_trim
parallel -j 4 "
    fastqc -t 4 -o ../fastqc_trim {1}
" ::: $(ls *.gz)
cd ../fastqc_trim
multiqc .
```
# 5 rm rRNA[optional]

# 6 seq alignment
hisat2
## 6.1 build index(.ht2)
```
hisat2-build [选项] [基因组序列(.fa)] [索引文件的前缀名]

cd ~/project/rat/genome
mkdir index
cd index

hisat2-build -p 6 ../rn6.chr1.fa rn6.chr1
```
## 6.2 alignment(.sam--.bam)
cd ~/project/rat/output/align
```
hisat2 [选项] -x [索引文件] [ -1 1测序文件 -2 2测序文件 -U 未成对测序文件 ] [ -S 输出的sam文件 ]

cd ~/project/rat/output
mkdir align
cd rRNA
parallel -k -j 4 "
    hisat2 -t -x ../../genome/index/rn6.chr1 \ 
    -U {1}.fq.gz -S ../align/{1}.sam \ 
    2>../align/{1}.log
" ::: $(ls *.gz | perl -p -e 's/.fq.gz$//')

cat SRR2190795.log
```
* summarize alignment rate and time
```
file_list=($(ls *.log))

echo -e "sample\tratio\ttime"
for i in ${file_list[@]};
do
    prefix=$(echo ${i} | perl -p -e 's/\.log//')
    echo -n -e "${prefix}\t"
    cat ${i} | 
        grep -E "(overall alignment rate)|(Overall time)" |
        perl -n -e '
            if(m/alignment/){
                $hash{precent} = $1 if m/([\d.]+)%/;
            }elseif(m/time/){
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
* convert format and sort(.bam)
samtools
```
parallel -k -j 4 "
    samtools sort -@ 4 {1}.sam > {1}.sort.bam
    samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
rm *.sam
ls
```
# 7 expression level
cd ~/project/rat/output 
```
htseq-count [options] <alignment_files> <gff_file>

mkdir HTseq
cd align
parallel -j 4 "
    htseq-count -s no -f bam {1}.sort.bam ../../annotation/rn6.gff \ 
      >../HTseq/{1}.count 2>../HTseq/{1}.log
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

cd HTseq
cat SRR2190795.count | head -n 10
```
# 8 
## 8.1 merge
```
rm(list=ls())
setwd("~/project/rat/output/HTseq")

files <- list.files(".", "*.count")
f_lists <- list()
for(i in files){
  prefix = gsub("(_\\w+)?\\.count","",i,perl=TRUE)
  f_lists[[prefix]] = i
}

id_list <- names(f_lists)
data <- list()
count <- 0
for(i in id_list){
  count <- count + 1
  a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))
  data[[count]] <- a
}

data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
  data_merge <- merge(data_merge, data[[i]], by = "gene_id")
}

write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)
```
## 8.2 read counts Normalization
calculate
* gene length
```
library(GenomicFeatures)
# Granges对象
txdb <- makeTxDbFromGFF("rn6.gff")
# exon
exons_gene <- exonsBy(txdb, by = "gene")
# full length
gene_len <- list()
for(i in names(exons_gene)){
  range_info = reduce(exons_gene[[i]])
  width_info = width(range_info)
  sum_len = sum(width_info)
  gene_len[[i]] = sum_len
}
# gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

data <- t(as.data.frame(gene_len))
#outfile
write.table(data, file = "rn6_gene_len.tsv", row.names =TRUE, sep = "\t", quote = FALSE, col.names = FALSE)