# nathan dot lazar at gmail dot com

# Runs both an occupancy analysis and differential binding analysis
# on samples listed in the input file.
# A number of plots are generated which are stored in the supplied 
# output directory

# Usage: run_DiffBind.R <sample file> <fragment length> <output directory>
# The sample sheet should be a comma separated file with a row for 
# each sample. Fields include Sample ID, Tissue, Factor,
# Condition, Treatment, Replicate, Read .bam file ,ControlID,
# Control .bam file, Peak .bed.gz file

# Example: samples.csv 150 all_out

library(DiffBind)

args <- commandArgs(TRUE)
samples.csv <- args[1]
frag.size <- args[2]
out.dir <- args[3]
# samples.csv <- 'H3K4_analysis_MINUS12/H3K4me3_minus12.csv'
# frag.size <- 150
# out.dir <- 'H3K4_analysis/'

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

# Read in sample sheet (which was written by hand)
samples <- read.csv(samples.csv)

# Merge peak sets to get a master list of peaks
# and get a binding matrix with a value for each peak
# and sample. The values are the confidence of the peak
# call and range from 0 to 1. If a peak is not called in a 
# sample it has a value of -1.
samples.dba = dba(sampleSheet=samples.csv, peakCaller='raw', minOverlap=1)
samples2.dba = dba(sampleSheet=samples.csv, peakCaller='raw')

############################################################
# 1) Occupancy analysis: Look at differences in called peaks
#    between samples and conditions

# Get counts of peaks shared in >=n samples
olap.rate <- dba.overlap(samples.dba, mode=DBA_OLAP_RATE)
print(olap.rate)
pdf(paste0(out.dir, 'olap_rate.pdf'))
plot(olap.rate, type='b', ylab='# peaks', 
     xlab='Overlap at least this many peaksets')
dev.off()

olap.alc <- dba.overlap(samples.dba, samples.dba$masks$Alcohol, mode=DBA_OLAP_RATE)
pdf(paste0(out.dir, 'olap_rate_alc.pdf'))
plot(olap.alc, type='b', ylab='# peaks', 
     xlab='Overlap at least this many peaksets')
dev.off()
olap.con <- dba.overlap(samples.dba, samples.dba$masks$NonAlc, mode=DBA_OLAP_RATE)
pdf(paste0(out.dir, 'olap_rate_con.pdf'))
plot(olap.con, type='b', ylab='# peaks', 
     xlab='Overlap at least this many peaksets')
dev.off()

# Venn diagrams for the two conditions separately
pdf(paste0(out.dir, 'Alc_venn.pdf'))
dba.plotVenn(samples.dba, samples.dba$masks$Alcohol, main='Alcohol')
dev.off()
pdf(paste0(out.dir, 'Con_venn.pdf'))
dba.plotVenn(samples.dba, samples.dba$masks$NonAlc, main='NonAlc')
dev.off()

# Generate a peakset for each condition where peaks are in at least
# two samples 
consensus.dba <- dba.peakset(samples.dba, consensus = -DBA_REPLICATE,
                             minOverlap=1)
consensus2.dba <- dba.peakset(samples.dba, consensus = -DBA_REPLICATE,
                              minOverlap=2)

# Venn diagram between conditions (the DiffBind Venn is off by a bit
# I think due to merging).
pdf(paste0(out.dir, 'Alc_vs_Con_Venn_consensus.pdf'))
dba.plotVenn(consensus.dba, mask = consensus.dba$masks$Consensus)
dev.off()
pdf(paste0(out.dir, 'Alc_vs_Con_Venn_consensus2.pdf'))
dba.plotVenn(consensus2.dba, mask = consensus2.dba$masks$Consensus)
dev.off()

# Do some counting by hand.
# The number of peaks in two or more samples that are in two 
# alcohol samples and NOT in two control samples
all.in2 <- which((consensus2.dba$allvectors$`3_151` != -1) +
                 (consensus2.dba$allvectors$`5_171` != -1) +
                 (consensus2.dba$allvectors$`7_211` != -1) +
                 (consensus2.dba$allvectors$`4_161` != -1) +
                 (consensus2.dba$allvectors$`6_201` != -1) +
                 (consensus2.dba$allvectors$`8_221` != -1) >=2)
print(sprintf('There are %d peaks found in 2 or more samples', 
              length(all.in2)))

alc.in2 <- which((consensus2.dba$allvectors$`3_151` != -1) +
                 (consensus2.dba$allvectors$`5_171` != -1) +
                 (consensus2.dba$allvectors$`7_211` != -1) >=2)
print(sprintf('There are %d peaks found in 2 or more alcohol samples', 
              length(alc.in2)))

con.in2 <- which((consensus2.dba$allvectors$`4_161` != -1) +
                 (consensus2.dba$allvectors$`6_201` != -1) +
                 (consensus2.dba$allvectors$`8_221` != -1) >=2)
print(sprintf('There are %d peaks found in 2 or more control samples', 
              length(con.in2)))

print(sprintf('of the %d peaks found in 2 or more alcohol samples,
 %d of them are found in two or more control samples', length(alc.in2),
 sum(alc.in2 %in% con.in2)))

print(sprintf('of the %d peaks found in 2 or more control samples,
 %d of them are found in two or more alcohol samples', length(con.in2),
 sum(alc.in2 %in% con.in2)))

# The number of peaks in two or more samples that are in two 
# control samples and NOT in two or more alcohol samples.

# Get the number of peaks shared by 2 or more samples for each sample
apply(consensus.dba$allvectors, 2, function(x) sum(x!=-1))
apply(consensus2.dba$allvectors, 2, function(x) sum(x!=-1))

# Get a list of GRanges objects of peaks unique to each group
# samples.GRlist = dba.overlap(consensus.dba, consensus.dba$masks$Consensus)
# samples.GRlist$onlyA
# samples.GRlist$onlyB
# samples.GRlist$inAll

# Using only peak caller data make a correlation heatmap.
# This gives an initial clustering of the samples based on 
# "occupancy scores" using the cross-correlations of each 
# row of the binding matrix:
pdf(paste0(out.dir, 'peak_cc_heatmap_samples.pdf'))
plot(samples2.dba)
dev.off()

pdf(paste0(out.dir, 'peak_cc_heatmap_consensus.pdf'))
plot(consensus.dba)
dev.off()

############################################################
# 2) Differential binding analysis

# Count reads in each peak for each sample
counts.dba <- dba.count(samples.dba, minOverlap=1, fragmentSize=frag.size)
# Only count reads for peaks that occur in 2 or more samples (minOverlap)
counts2.dba <- dba.count(samples.dba, minOverlap=2, fragmentSize=frag.size)

# Make correlation heatmaps using count data "affinity scores"
pdf(paste0(out.dir, 'count_cc_heatmap.pdf'))
plot(counts.dba)
dev.off()

pdf(paste0(out.dir, 'count_cc_heatmap2.pdf'))
plot(counts2.dba)
dev.off()

# Set up contrasts between conditions
counts.dba <- dba.contrast(counts.dba, categories=DBA_CONDITION)
counts2.dba <- dba.contrast(counts2.dba, categories=DBA_CONDITION)

# Make condition a factor so 'NonAlc' will be the reference state
# counts2.dba$samples$Condition <- factor(counts2.dba$samples$Condition, 
#   levels=c('NonAlc', 'Alcohol')
# counts2.dba$samples$groups <- factor(counts2.dba$samples$Condition, 
#   levels=c('NonAlc', 'Alcohol')

# Perform differential binding analysis 
counts.dba <- dba.analyze(counts.dba, method=DBA_ALL_METHODS)
counts2.dba <- dba.analyze(counts2.dba, method=DBA_ALL_METHODS)

# Make correlation heatmaps using only differentially bound peaks
# with default threshold of FDR >= 0.1
pdf(paste0(out.dir, 'diff_cc_heatmap.pdf'))
plot(counts.dba, contrast=1)
dev.off()

pdf(paste0(out.dir, 'diff_cc_heatmap2.pdf'))
plot(counts2.dba, contrast=1)
dev.off()

# Retrieve differentially bound sites in GRanges objects
# with an FDR cutoff of 0.1
diff.gr <- dba.report(counts.dba)
diff2.gr <- dba.report(counts2.dba)

# Retrieve all sites bound in 2 sorted by FDR
# th=1 sets the FDR threshold to 1 giving all regions
# bCalled=T show which samples contain this peak
# bCounts=T tells to report counts for all samples
# bCalledDetail? peak caller status for each sample...
# bAll include peaksets combining peaks w/ both + and - fold changes

diff.full.gr <- dba.report(counts.dba, th=1,
  bCalled=T, bCounts=T, bCalledDetail=T)
diff2.full.gr <- dba.report(counts2.dba, th=1,
  bCalled=T, bCounts=T, bCalledDetail=T)

# Add in names and location of the nearest gene and distance to 
# the nearest gene (0 if in a gene).
genes <- read.table('mm10_refseq.bed')
genes.gr <- makeGRangesFromDataFrame(genes, keep.extra.columns=T,
  seqnames.field='V1', start.field='V2', end.field='V3')

# Add info on the nearest gene to peaks
add_gene_info <- function(gr, genes.gr) {
  gr.dist <- distanceToNearest(gr, genes.gr)
  dist.df <- as.data.frame(gr.dist)
  gr$gene.id <- NA
  gr$gene.chr <- NA
  gr$gene.start <- NA
  gr$gene.end <- NA

  gr$dist.to.gene[dist.df$queryHits] <- mcols(gr.dist)$distance
  gr$gene.id[dist.df$queryHits] <- as.vector(genes.gr$V4[dist.df$subjectHits])
  gr$gene.chr[dist.df$queryHits] <- 
    as.vector(seqnames(genes.gr)[dist.df$subjectHits])
  gr$gene.start[dist.df$queryHits] <- 
    as.vector(start(genes.gr)[dist.df$subjectHits])
  gr$gene.end[dist.df$queryHits] <- 
    as.vector(end(genes.gr)[dist.df$subjectHits])

  return(gr)
}

# Add in gene info for each type of overlap
diff.gr <- add_gene_info(diff.gr, genes.gr)
diff2.gr <- add_gene_info(diff2.gr, genes.gr)
diff.full.gr <- add_gene_info(diff.full.gr, genes.gr)
diff2.full.gr <- add_gene_info(diff2.full.gr, genes.gr)

# Write all peaks to csv
write.table(as.data.frame(diff.full.gr), file=paste0(out.dir, 'H3K4_peaks_minus12.csv'), sep=',', 
  quote=F, row.names=F, col.names=T)
write.table(as.data.frame(diff2.full.gr), file=paste0(out.dir, 'H3K4_peaks_in2_minus12.csv'), sep=',', 
  quote=F, row.names=F, col.names=T)

# Write differentially bound GRanges objects to csv
write.table(as.data.frame(diff.gr), file=paste0(out.dir, "H3K4_diff_bound_minus12.csv"), 
  quote=F, sep=",", row.names=F, col.names=T)
write.table(as.data.frame(diff2.gr), file=paste0(out.dir, "H3K4_diff_bound_in2_minus12.csv"), 
  quote=F, sep=",", row.names=F, col.names=T)

# Plots
####################################

# Make MA plots
pdf(paste0(out.dir, 'MA_plot.pdf'))
dba.plotMA(counts.dba)
dev.off()

pdf(paste0(out.dir, 'MA_plot2.pdf'))
dba.plotMA(counts2.dba)
dev.off()

# PCA plots
pdf(paste0(out.dir, 'pca.pdf'))
dba.plotPCA(counts.dba, DBA_CONDITION, label=DBA_REPLICATE)
dev.off()

pdf(paste0(out.dir, 'pca2.pdf'))
dba.plotPCA(counts2.dba, DBA_CONDITION, label=DBA_REPLICATE)
dev.off()

# PCA using just differentially bound sites
pdf(paste0(out.dir, 'pca_diff.pdf'))
dba.plotPCA(counts.dba, DBA_CONDITION, contrast=1, label=DBA_REPLICATE)
dev.off()

pdf(paste0(out.dir, 'pca_diff2.pdf'))
dba.plotPCA(counts2.dba, DBA_CONDITION, contrast=1, label=DBA_REPLICATE)
dev.off()

# Heatmaps showing differentially bound loci
pdf(paste0(out.dir, 'DB_heatmap.pdf'))
dba.plotHeatmap(counts.dba, contrast=1, correlations=FALSE)
dev.off()

pdf(paste0(out.dir, 'DB_heatmap2.pdf'))
dba.plotHeatmap(counts2.dba, contrast=1, correlations=FALSE)
dev.off()

# RPKM heatmaps
pdf(paste0(out.dir, 'DB_heatmapRPKM.pdf'))
dba.plotHeatmap(counts.dba, score=DBA_SCORE_RPKM_FOLD)
dev.off()

pdf(paste0(out.dir, 'DB_heatmapRPKM2.pdf'))
dba.plotHeatmap(counts2.dba, score=DBA_SCORE_RPKM_FOLD)
dev.off()

save.image(paste0(out.dir, 'all.Rdata'))
