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
# Add my bwa-0.7.17/ folder to PATH
export PATH="$HOME/bwa-0.7.17:$PATH"
# run bwa index to create index on hg19 genome file
# You don't need to run this - it's already created
# Just included so you know that this is a step
bwa index ~/var_data/GATKBundle/ucsc.hg19.fasta
```

```
# run bwa mem on our fastq file
# run this in your student folder
bwa mem -Y -t 12 -R '@RG\tID:SRR702072\tPL:Illumina\tLB:SRR702072\tSM:SRR702072' ~/var_data/GATKBundle/ucsc.hg19.fasta ~/var_data/fastq/SRR702072_1.filt.fastq ~/var_data/fastq/SRR702072_2.filt.fastq -o SRR702072.sam
```

# Samtools

Install samtools and add to your `$PATH` - more instructions here: http://www.htslib.org/download/

```
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
```

Run `samtools flagstaff` on our `.sam` file:

```
samtools flagstat SRR702072.sam
```


# GATK

Download, unzip with `unzip`, and add to your `$PATH`

```
wget https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip
```

```
gatk --list
gatk MarkDuplicates --help

```

GATK Best Practices Tutorials: https://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk-1

HaplotypeCaller: https://docs.google.com/file/d/0B2dK2q40HDWeYzVTUGs0bjM3M1E/preview
