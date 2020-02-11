## FASTQ Files Location

`/home/courses/BMI535/data/VariantCalling/fastq/`

## Exome Capture Target Location

`/home/courses/BMI535/data/VariantCalling/1000GenomesExomeCapture/`

## GATK4 Best Practices

https://github.com/gatk-workflows/gatk4-data-processing

## GATK Bundle

`/home/courses/BMI535/data/VariantCalling/GATKBundle`

## FastQC

- https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Download: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip

```
# add my FastQC directory to the PATH variable
export PATH="$HOME/FastQC/:$PATH"
# create a symbolic link called ~/var_data to easily access it
ln -s /home/courses/BMI535/data/VariantCalling/ ~/var_data
fastqc -o ~/var_data/fastq/SRR702072_1.filt.fastq 
```

## BWA

Main page: http://bio-bwa.sourceforge.net/

```
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
bzip2 -d bwa-0.7.17.tar.gz
tar -xvf bwa-0.7.17.tar
cd bwa-0.7.17
more README.md
# compile the code
make
```

```
# Add my bwa-0.7.17/ to PATH
export PATH="$HOME/bwa-0.7.17:$PATH"
# run bwa index
bwa index ~/var_data/GATKBundle/ucsc.hg19.fasta
```
