#!/bin/bash

bin=/opt/installed
bin2=/home/exacloud/lustre1/users/lazar/bin/bedtools2/bin
cores=24
REF=$1
F1=$2
NM=$3

# This script runs bwa-mem and filters for quality of 30 (which excludes multi-mapping reads)
# Results are placed in MAPPED
# Example outputs are:
#  reads.bam
#  reads.q30.sort.bam
#  reads.q30.sort.bam.bai
#  reads.mapped_reads   # Read counts at each step and summary measures


# Get basic stats from the fasta file
echo 'total unique per_unique most_common_seq most_com_seq_count per_most_common' > $NM.mapped_reads
zcat $F1 | awk '(NR%2==0){read=$1;total++;count[read]++} \
  END{for(read in count){if(!max||count[read]>max) \
    {max=count[read];maxRead=read};\
    if(count[read]==1){unique++}}; \
  print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}' >> $NM.mapped_reads

# Map reads, filter out secondary alignments and unmapped reads
$bin/bwa mem -M -t $cores $REF $F1 > $NM.sam 
$bin/samtools view -bShF0x0104 $NM.sam > $NM.bam
# -M marks shorter split hits as secondary
# -t is the number of threads
# -b output bam
# -S input is sam
# -h print sam header
# -F 0x0100 skip secondary hits 0x0004 skip unmapped 

# Record number of mapped reads
echo -n "Uniquely mapped reads: " >> $NM.mapped_reads
$bin/samtools view $NM.bam | wc -l >>  $NM.mapped_reads

# Filter by mapping quality scores over 30
$bin/samtools view -b -q 30 $NM.bam  > $NM.q30.bam
echo -n "Passing Q 30 filter: " >> $NM.mapped_reads
$bin/samtools view $NM.q30.bam | wc -l >> $NM.mapped_reads

# Sort
$bin/samtools sort -@ $cores $NM.q30.bam $NM.q30.sort

# Index the filtered and sorted .bam file
$bin/samtools index $NM.q30.sort.bam

echo "#########################" >> $NM.mapped_reads

# Count mapping coordinates, unique mapping coordinates & maximum number of reads
# mapped to the same position
echo "Mapping positions:" >> $NM.mapped_reads
echo "Total Unique PerUnique MaxCoor CountMaxCoor PerMaxCoor" >> $NM.mapped_reads
$bin2/bamToBed -i $NM.q30.sort.bam | \
awk '{coordinates=$1":"$2"-"$3; total++; count[coordinates]++}\
  END{for(coordinates in count){if(!max||count[coordinates]>max) \
    {max=count[coordinates];maxCoor=coordinates}; \
    if(count[coordinates]==1){unique++}}; \
  print total,unique,unique*100/total,maxCoor,count[maxCoor],count[maxCoor]*100/total}' >> $NM.mapped_reads

# The non-redundant read fraction (NRF in ENCODE) is the third field.
# This is the ratio between the number of positions in the genome that 
# uniquely mappable reads map to and the total number of uniquely mappable reads

# Clean up intermediate files
rm $NM.sam
rm $NM.q30.bam

###############################################################################
# Depreciated
# echo -n "Number of unique mapping start positions:" >> $NM.mapped_reads
# $bin/samtools view $NM.q30.sort.bam | awk -vOFS='\t' '{print $3,$4,$6}' | sort -u | wc -l >> $NM.mapped_reads
#
# echo "## AFTER samtools rmdup -s ##" >> $NM.mapped_reads
#
# Remove duplicates w/ samtools rmdup
# $bin/samtools rmdup -s $NM.q30.sort.bam  $NM.q30.sort.rmdup.bam
# echo -n "After samtools rmdup: " >> $NM.mapped_reads
# raw=$($bin/samtools view $NM.q30.sort.rmdup.bam | wc -l)
# echo $raw >> $NM.mapped_reads
#
# Count the final mapped reads, unique read coordinates & maximum number of reads
# mapped to the same position
# echo "Total Unique PerUnique MaxCoor CountMaxCoor PerMaxCoor" >> $NM.mapped_reads
# $bin2/bamToBed -i $NM.q30.sort.rmdup.bam | \
# awk '{coordinates=$1":"$2"-"$3; total++; count[coordinates]++}\
#   END{for(coordinates in count){if(!max||count[coordinates]>max) \
#     {max=count[coordinates];maxCoor=coordinates}; \
#     if(count[coordinates]==1){unique++}}; \
#   print total,unique,unique*100/total,maxCoor,count[maxCoor],count[maxCoor]*100/total}' >> $NM.mapped_reads
#
# Calculate NRF (non-redundant read fraction) as in ENCODE. This is the ratio between
# the number of positions in the genome that uniquely mappable reads map to and
# the total number of uniquely mappable reads
# echo -n "Number of unique mapping start positions:" >> $NM.mapped_reads
# $bin/samtools view $NM.q30.sort.rmdup.bam | awk -vOFS='\t' '{print $3,$4,$6}' | sort -u | wc -l >> $NM.mapped_reads

