#!/bin/bash
cd /you/working/dir

/your/star/dir/STAR --runMode genomeGenerate --runThreadN 14 --genomeDir /place/you/want/to/put/genome/index --genomeFastaFiles your.fasta --sjdbGTFfile your.gtf
