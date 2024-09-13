# mkdir
```mkdir biosoft
mkdir -p project/rat
cd project/rat
mkdir annotation genome sequence output script
```
# install tools
## sra
download and convert data
```cd ~/biosoft
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
tar xvfz sratoolkit.2.11.3-ubuntu64.tar.gz
mv sratoolkit.2.11.3-ubuntu6 sratoolkit.2.11.3
cd sratoolkit.2.11.3/bin

#set PATH
export PATH="$(pwd):$PATH"

#or:brew install sratoolkit
```
## fastqc
```
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
uzip fastqc_v0.12.1.zip

#brew install fastqc
```
