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
unzip sratoolkit.2.11.3-win64.zip
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
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
~/biosoft/TrimGalore-0.6.10/trim_galore
export PATH="$(pwd):PATH"
source ~/.bashrc
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
tar -xzvf stringtie-2.2.3.tar.gz
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
wget https://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
```
* decompress
* rename
```
mv Rattus_norvegicus.Rnor_7.2.dna.toplevel.fa rn7.raw.fa
```
* check(chr...)
* clear some additional text
```
cat rn7.raw.fa | perl -ne 'if(m/^>(.+?)(?:\s|$)/){ print ">$1\n";}else{print}' > rn7.fa
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
* NCBI-GEO >> GEO accession (GSE72960)
* SRA Run Selector (SRP063345)
* download
```
nohup prefetch SRR2190795 SRR224018{2..7} SRR2240228 -O . &
```
* convert format(.sra > .fastq > .gz)
```
mkdir srr
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
fastqc -t 6 -o ../output/fastqc *.gz

cd ~/project/rat/output/fastqc
ls
```
* multiqc
```
multiqc .
```
解读：https://www.jianshu.com/p/db2100baabb5
  * per seq GC content
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
> Hisat2是一种用于转录组数据比对的软件
* Hisat2使用两种类型的索引，全局索引和局部索引。全局索引基于Burrows-Wheeler Transform (BWT)和Ferragina-Manzini (FM) index的方法，可以快速定位reads在基因组上的大致位置；局部索引基于哈希表的方法，可以对reads进行精确的扩展和比对。这样可以提高比对的速度和敏感性，同时节省内存空间。
* BWT是一种线性时间复杂度的排序算法，它通过重新排列输入字符串的字符来创建一个新的字符串，这个新字符串保留了原始字符串的排序信息。
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
## 7.1 HTseq
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
## 7.2 StringTie
> StringTie是一个用于转录组组装和定量的软件，它结合了有参考基因组的转录本拼接和无参考基因组的从头组装方法。
* 首先将reads比对到参考基因组，然后根据定位的坐标信息和跨越内含子的spliced reads中反映的连接关系建立备选的拓扑图。接着，设计相应的算法在拓扑图中选择合理的转录本形成最终的转录组数据集。
* StringTie提供了多种表达量的估计方法，包括raw count、FPKM（每百万reads的fragments每千个碱基的比率）和TPM（每百万reads的transcripts每千个碱基的比率）。
* StringTie还提供了一个合并模式（--merge），它可以将多个样本的转录本组装结果合并成一个非冗余的转录本集合，这在处理多个RNA-Seq样本时非常有用。

input:align (.bam);annotation (.gff) <br>
output:assembly(.gtf)
```
cd ~/project/rat/output
mkdir assembly

cd align
parallel -j 4 "
  stringtie -p 4 -G ../../annotation/rn7.gff -l {1} {1}.sort.bam -o ../assembly/{1}.gtf
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

#merge all samples
cd ../assembly
ls *.gtf > mergelist.txt

stringtie --merge -p 8 -G ../../annotation/rn7.gff -o merge.gtf mergelist.txt

#abundance
cd ~/project/rat/output
mkdir abundance

cd align
parallel -j 4 "
  mkdir ../abundance/{1}
  stringtie -e -B -p 4 -G ../assembly/merge.gtf -l {1} {1}.sort.bam -o ../assembly/{1}/{1}.gtf
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

#gffcompare

#expression matrix
cd ~/project/rat/script
wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py
python2 prepDE.py --help

cd ~/project/rat/output/
mkdir matrix

python2 ~/project/rat/script/prepDE.py \
   -i ./abundance \
   -g ./matrix/gene_count_matrix.csv \
   -t ./matrix/transcript_count_matrix.csv \
   -l 100
```
# 8 merge & normalize
## 8.1 merge HTseq-count results
input:HTseq (*.count) <br>
output:merge.csv
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
## 8.2 normalize
* same gene from different samples
CPM
* different genes from one sample
RPKM, FPKM, TPM <br>

CPM = (10^6 *nr) / N <br>
RPKM = (10^6 *nr) / (L * N) <br>
  * RPKM: Reads Per Kilobase per Million
  * nr : 比对至目标基因的read数量
  * L : 目标基因的外显子长度之和除以1000(因此，要注意这里的L单位是kb，不是bp)
  * N : 是总有效比对至基因组的read数量
TPM = (nr / g_r) * 10^6 / ∑(ni / gi)
  * TPM : Transcripts Per Million
  * nr : 比对至目标基因的read数量
  * read_r: 是比对至基因r的平均read长度
  * g_r : 是基因r的外显子长度之和（这里无需将其除以1000）
### 8.2.1 cufflinks
### 8.2.2 manual calculation
* gene_length
input:annotation (.gff) <br>
output:HTseq (rn7_gene_len.tsv)
```
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("rn7.gff")

exons_gene <- exonsBy(txdb, by= "gene")
gene_len <- list()
for (i in names(exons_gene)){
  range_info = reduce(exons_gene[[i]])
  width_info = width(raneg_info)
  sum_len = sum(width_info)
  gene_len[[i]] = sum_len
}
#lapply
gene_len <- lapply(exons_gene, function(x{sum(width(reduce(x)))}))

data <- t(as.data.fram(gene_len))
write.table(data, file = "rn7_gene_len.tsv", row.names = TRUE, sep="\t", quote = FALSE, col.names = FALSE)
```
* RPKM
RPKM = (10^6 *nr) / (L * N) <br>
input:HTseq (rn7_gene_len.tsv);(*.count) <br>
output: count["RPKM"]
```
gene_len_file <- "rn7_gene_len.tsv"
gene_len <- read.table(gene_len_file, header =FALSE, row.name = 1)
colnames(gene_len) <- c("length")

count_file <- "SRR2190795.count"
count <- read.table(count_file, header = FALSE, row.name = 1)
colnames(count) <- c("count")
all_count <- sum(count["count"])

RPKM <- c()
for (i in row.names(count)){
  count = 0
  rpkm = 0
  exon_kb =1
  count_ = count[i, ]
  exon_kb = gene_len[i, ] / 1000
  rpkm = (10 ^ 6 * count_ ) / (exon_kb * all_count )
  RPKM = c(RPKM, rpkm)
}
count["RPKM"] <- RPKM
```
* TPM
TPM = (nr / g_r) * 10^6 / ∑(ni / gi) <br>
input:gene-len;*.count <br>
output:count["TPM"];"*.normalize.count"
```
sum_ <- 0
for(i in row.names(count)){
    count_ = 0
    exon = 1
    count_ = count[i, ]
    exon  = gene_len[i, ]
    value = count_ / exon
    na_values <- is.na(value)
    if(is.na(value)){
        print(paste(i, "is error! please check"))
    }else{
        sum_ = sum_ + value
    }
}

TMP <- c()
for (i in row.names(count)){
    count_ = 0
    tpm = 0
    exon = 1
    count_ = count[i, ]
    exon = gene_len[i, ]
    tpm = (10 ^ 6 * count_ / exon ) / sum_
    TPM = c(TPM, tpm)
}
count["TPM"] <- TPM

write.table(count, "SRR2190795.normalize.count", col.names = TRUE, row.names = TRUE, sep="\t", quote = FALSE)
```
# 9 differential expression
## 9.1 pre-treatment
```
dataframe <- read.csv("merge.csv, header = TRUE, row.name = 1)

#delete head 5 lines
countdata <- dataframe[-(1:5), ]
head(countdata)

# delete ID version-number
row_names <- row.names(countdata)
new_row_names <- gsub("\\.\\w+","", row.names(countdata))
row.names(countdata) <- new_row_names

# remove low-expression gene
countdata <- countdata[rowSum(countdata) > 0, ]
```
## 9.2 differential expression
DEseq2(raw data--HTseq--*.count)
> 构建一个dds(DESeqDataSet)的对象
> 利用DESeq函数进行标准化处理
> 用result函数来提取差异比较的结果
### 9.2.1 download packages
```
#download
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("biomaRt")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("clusterProfiler")

library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Rn.eg.db)
library(clusterProfiler)
```
### 9.2.2 create DESeqDataSet(dds)
```
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ batch + condition)
```
* countData
merge.csv--pre-treatment(countdata)
* colData
input:sequence (phenotype.csv)
* design
phenotype.csv--treatment
```
#bash
cat <<EOF >./phenotype/phenotype.csv
"ids","state","condition","treatment"
"SRR2240185","Liver cirrhosis","DEN","treatment"
"SRR2240186","Liver cirrhosis","DEN","treatment"
"SRR2240187","Healthy control","PBS","control"
"SRR2240228","Healthy control","PBS","control"
EOF
```
R
```
countdata

coldata <- read.table("../phenotype/phenotype.csv", row.names = 1, header = TRUE, sep = "," )
head(coldata)
countdata <- countdata[row.names(coldata)]

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
dds
```
### 9.2.3 sample correlation
* PCA
rlog vst<br>
input:dds(countdata,coldata,design)
output:PCA
```
vsdata <- rlog(dds, blind=FALSE)
plotPCA(vsdata, intgroup="treatment") + ylim(-10, -10)
```
* sample-to-sample distances heat map
input:vsdata(dds)
output:heat map
```
library("RColorBrewer")
gene_data_transform <- assay(vsdata)
sampleDists <- dist(t(gene_data_transform))
sampleDistsMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")
# colnames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blue")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
### 9.2.4 differential expression gene
difference between different samples
input:dds
```
dds$treatment <- factor(as.vector(dds$treatment), levels = c("control","treatment"))
dds <- DESeq(dds)
```
> 估计样本大小（消除测序深度的影响）
> 对每一个基因的表达量的离散度做计算
> 拟合广义的线性模型（generalized linear model）
```
result <- results(dds, pAdjustMethod = "fdr", alpha = 0.05)
head(result)
result_order <- result[order(result$pvalue),]
head(result_order)
```
output:
```
log2 fold change (MLE): treatment treatment vs control 
Wald test p-value: treatment treatment vs control 
DataFrame with 6 rows and 6 columns
                           baseMean    log2FoldChange             lfcSE              stat               pvalue                 padj
                          <numeric>         <numeric>         <numeric>         <numeric>            <numeric>            <numeric>
ENSRNOG00000011250 208.046881231885 -7.50369356010508  0.44485821990427 -16.8676068562245 7.78886122158816e-64 1.14472893373681e-59
ENSRNOG00000047945 3799.51961509786  4.50434780195392 0.358671277660941  12.5584290755837 3.57357467823096e-36 2.62604135229802e-32
ENSRNOG00000017897 1130.41206772961  4.41361091204456 0.353924586455456  12.4704840549416 1.08166978868575e-35 5.29910029477147e-32
ENSRNOG00000001466 542.805654192746  4.87418957369525 0.412058420866664  11.8288799035913 2.76810877295631e-32 1.01707236590347e-28
ENSRNOG00000014013 400.690803761036  2.83690340404308 0.246440071910237  11.5115345570764 1.15406329271928e-30 3.39225364261904e-27
ENSRNOG00000018371 705.943312960284  4.65111741552834  0.41021741987017  11.3381762700384  8.4895191640983e-30 2.07950771924588e-26
```
analysis:
* upregulate/downregulate
```
summary(result_order)
```
output:
```
out of 19962 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2248, 11%
LFC < 0 (down)     : 1148, 5.8%
outliers [1]       : 31, 0.16%  #异常值
low counts [2]     : 5234, 26%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```
* significant gene number
```
table(result_order$padj<0.05)
```
output:
```
FALSE  TRUE 
11301  3396
```
* save results
```
dir.create("../DESeq2") #output
write.csv(result, file="../DESeq2/results.csv", quote = F)
```
# 10 find diff_gene
## 10.1 diff_gene
input:result_order
```
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)
write.csv(diff_gene, file="../DESeq2/difference.csv", quote = F)
```
## 10.2 convert ensembl gene ID(ClusterProfiler)
```
ensembl_gene_id <- row.names(diff_gene)
ensembl_id_transform <- function(ENSEMBL_ID){
    # geneID是输入的基因ID，fromType是输入的ID类型，toType是输出的ID类型，OrgDb注释的db文件，drop表示是否剔除NA数据
    a = bitr(ENSEMBL_ID, fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Rn.eg.db")
    return(a) 
}
ensembl_if_transform(ensembl_gene_id)
```
## 10.3 annotation(biomaRt)
```
#choose database & get symbols
mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
ensembl_gene_id <- row.names(diff_gene)
rat_symbols <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)

#merge diff_gene & symbols
diff_gene$ensembl_gene_id <- ensembl_gene_id
diff_gene_dataframe <- as.data.frame(diff_gene)
diff_gene_symbols <- merge(diff_gene_dataframe, rat_symbols, by= c("ensembl_gene_id"))

#save data
write.table(result, "../stat/all_gene.tsv", sep="\t", quote = FALSE)
write.table(diff_gene_symbols, "../stat/diff_gene.tsv", row.names = F,sep="\t", quote = FALSE)

#count diff_gene
echo -e "sample\t\num" > all.samples.csv
for i in ${ls};
do
  if [-d ${i}];then
    prefix=$i
    diff_num=$(cat $i/diff_gene.tsv | tail -n+2 | wc -l)
    echo -e "${predix}\t${diff_num}" >> all_samples.tsv
    fi
done

#plot
library(ggplot2)
data <- read.table("all_sampels.tsv", head=T)

pdf("samples_diff_gene_num.pdf")
  ggplot(data=data, aes(x=sample, y=num, fill=samples)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x= "samples", y= "num", title= "different gene number")
dev.off()
```
# 11 visualize
* MA plot
```
plotMA(result_order, ylim=c(-10,10))
```
* heat map
# 12 enrichment analysis
* ClusterProfiler
```
ensembl_gene_id <- row.names(diff_gene)
rat_symbols <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"), filters='ensembl_gene_id', values=ensembl_gene_id, mart=mart)
diff_gene_ensembl_id <- rat_symbols$ensembl_gene_id
```
## 12.1 GO analysis
molecular function(MF), biological process(BP), cellular component(CC)
```
for(i in c("MF", "BP", "CC")){
    ego <- enrichGO(gene       = rat_symbols$entrezgene_id,
                    OrgDb      = org.Rn.eg.db,
                    keyType    = 'ENSEMBL',
                    ont        = i,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05)
    dotplot(ego, showCategory = 30, title = paste("The GO ", i, " enrichment analysis", sep = ""))
}
```
## 12.2 KEGG
```
kk <- enrichKEGG(gene = gene, 
                 organism ='rno',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data = FALSE)
```
## 12.3 GSEA
GSEA针对所有基因，KEGG针对差异基因。<br>
输入文件不是一个基因列表，而是除了基因之外还有它的表达量（目前样本中所有的非0的基因的倍数的改变的数值）。
## 12.4 DO（Disease Ontology）
## 12.5 ReactomePA
## 12.6 website
[metascape](https://metascape.org/gp/index.html#/main/step1)
[webgenstal](http://www.webgestalt.org/)
[DAVID](https://david-d.ncifcrf.gov/)
