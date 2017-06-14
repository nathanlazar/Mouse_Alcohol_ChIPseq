#!/bin/bash

bin=/home/exacloud/lustre1/users/lazar/bin/bedtools2/bin
bin2=/home/exacloud/lustre1/users/lazar/bin
bin3=/opt/installed

PATH=$PATH/:$bin/:$bin2/:$bin3; export PATH

# Usage: make_BigWig.sh <chrom_sizes> <input.bam> <fragment length>
# .density (.bed format) and .bw (BigWig) are created

# Note: bam file should be sorted and indexed

# Example: make_BigWig.sh MOUSE_mm10/chrom_lens.txt MAPPED/DNA151222MG_1_131_K4_S6.q30.sort.bam 150

chrom_sizes=$1
f=$2
extend=$3
sample=${f/.bam/}

librarySize=$(samtools idxstats ${sample}.bam | awk '{total+=$3}END{print total}')
bamToBed -i ${sample}.bam > ${sample}.bed

awk -vCHROM=$chrom_sizes -vEXTEND=$extend -vOFS='\t' \
 'BEGIN{while(getline<CHROM){chromSize[$1]=$2}} \
 {chrom=$1;start=$2;end=$3;strand=$6; \
 if(strand=="+"){end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}}; \
 if(strand=="-"){start=end-EXTEND;if(start<1){start=1}}; \
 print chrom,start,end}' ${sample}.bed | \
sort -k1,1 -k2,2n | \
genomeCoverageBed -i stdin -g $chrom_sizes -d | \
awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | \
gzip > ${sample}.density.gz

#awk -vCHROM=$chrom_sizes -vEXTEND=$extend -vOFS='\t' \
# 'BEGIN{while(getline<CHROM){chromSize[$1]=$2}} \
# {chrom=$1;start=$2;end=$3;strand=$6;if(strand=="+") \
#   {end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}}; \
# if(strand=="-"){start=end-EXTEND;if(start<1){start=1}};print chrom,start,end}' ${sample}.bed | \

# Create WIG file
gunzip -c ${sample}.density.gz | \
awk -vOFS='\t' '($4!=0) {if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | \
gzip > ${sample}.wig.gz

# Convert to BigWig file
wigToBigWig ${sample}.wig.gz MOUSE_mm10/chrom_lens.txt ${sample}.bw

# Clean up
mv ${sample}.density.gz COVERAGE/
mv ${sample}.bw COVERAGE/

rm ${sample}.wig.gz
