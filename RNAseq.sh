#!/bin/bash
# Usage: ./workflow.sh <prefix length> <data1_1> <data1_1> <data2_2> <data2_2> ...
# Add --readFilesCommand zcat in STAR argument options and -8 to -11 if data1_1 is of type .fastq.gz
cd /your/working/dir

# GFF="$1"
# shift

while [[ $# > 1 ]]
do

        FILEONE="$1"
        FILETWO="$2"
        shift 2

        prefix=${FILEONE:${#FILEONE}-16:7}
        echo ${#FILEONE}
        SAM=${prefix}Aligned.out.sam
        outfile=${prefix}.txt
        #apply star alignment
        /usr/local/star/latest/bin/STAR --genomeDir ../index --genomeLoad LoadAndKeep --runThreadN 14 --readFilesIn ${FILEONE} ${FILETWO} -$
        #sorting by using samtools
        /usr/local/samtools/latest/bin/samtools view -b $SAM > ${prefix}Aligned.out.bam
        /usr/local/samtools/latest/bin/samtools sort -m 1000000000 ${prefix}Aligned.out.bam $prefix
        /usr/local/samtools/latest/bin/samtools view -h ${prefix}.bam > ${prefix}Sorted.sam
        #apply HT-seq to do the seq count
        /usr/local/htseq/latest/bin/htseq-count --mode=intersection-nonempty --stranded=no --order=name --type=exon --idattr=gene_name -m u$

done
