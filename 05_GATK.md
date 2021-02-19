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
## Download FASTQC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 711 [fastqc directory]

# add my FastQC directory to the PATH variable
# you can add this to your .bash_profile
export PATH="$HOME/FastQC/:$PATH"
# create a symbolic link called ~/var_data to easily access the data
ln -s /home/courses/BMI535/data/VariantCalling/ ~/var_data

# create symbolic link
ln -s /home/courses/BMI535/students/laderast ~/BMI535/

#make an output directory
mkdir ~/qc_analyses

fastqc -o ~/qc_analyses ~/var_data/fastq/SRR702072_1.filt.fastq ~/var_data/fastq/SRR702072_2.filt.fastq
```

## BWA

Main page: http://bio-bwa.sourceforge.net/

```
conda install -c bioconda bwa
which bwa
```

```

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

Install samtools via conda

```
conda install -c bioconda samtools
```

Run `samtools flagstat` on our `.sam` file:

```
samtools flagstat SRR702072.sam
```


# GATK

Download, unzip with `unzip`, and add to your `$PATH`

```
wget https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip
```

```
#list all utilities
gatk --list
#convert SAM to BAM file
gatk SamFormatConverter -I SRR702072.sam -O SRR702072.bam
#sort BAM file using SortSam
gatk SortSam -I SRR702072.bam -O sortedSRR702072.bam -SO coordinate
```

Try running `gatk MarkDuplicates` on `sortedSRR702072.bam` by [reading the man page!](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)

```
gatk MarkDuplicates --help
```

## Where is GATK on exacloud?

```
/opt/installed/gatk-4.0.0.0/gatk
```

## Useful GATK Links

- GATK Command-line Syntax: https://gatk.broadinstitute.org/hc/en-us/articles/360035531892
- GATK Data Preprocessing for Variant Discovery (start here, thanks Mike!): https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
- GATK Best Practices Tutorials: https://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk-1
- HaplotypeCaller: https://docs.google.com/file/d/0B2dK2q40HDWeYzVTUGs0bjM3M1E/preview
