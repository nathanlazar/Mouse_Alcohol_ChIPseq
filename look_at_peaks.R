# nathan dot lazar at gmail dot com

# Read in bed files of peaks, make histograms of peak sizes 
# bar plot of the chromosomal distribution.

library(dplyr)
library(ggplot2)
library(vioplot)

alc <- read.table('IGV/Alc_pool_peaks.bed', stringsAsFactors = F,
                  col.names=c('chr', 'start', 'end', 'strand', 'idr_score'),
                  colClasses = c(character(0), numeric(0), numeric(0), 
                                 character(0), numeric(0)))
con <- read.table('IGV/Con_pool_peaks.bed', stringsAsFactors = F,
                  col.names=c('chr', 'start', 'end', 'strand', 'idr_score'),
                  colClasses = c(character(0), numeric(0), numeric(0), 
                                 character(0), numeric(0)))

alc$size <- alc$end-alc$start
con$size <- con$end-con$start

pdf('Peak_size_hist.pdf')
par(mfrow=c(2,2))
hist(alc$size, breaks=50)
hist(con$size, breaks=50)
hist(alc$size, breaks=500, xlim=c(550,650))
hist(con$size, breaks=500, xlim=c(500,600))
summary(alc$size)
summary(con$size)
# Test for a difference in the sizes of peaks 
wilcox.test(alc$size, con$size)

alc2 <- alc[!grepl('_', alc$chr),]
alc2$chr <- as.factor(alc2$chr)
boxplot(size~chr, data=alc2)

vioplot(alc2$size[alc2$chr=='chr1'], alc2$size[alc2$chr=='chr2'],
        alc2$size[alc2$chr=='chr3'], alc2$size[alc2$chr=='chr4'],
        alc2$size[alc2$chr=='chr5'], alc2$size[alc2$chr=='chr6'],
        alc2$size[alc2$chr=='chr7'], alc2$size[alc2$chr=='chr8'],
        alc2$size[alc2$chr=='chr9'], alc2$size[alc2$chr=='chr10'],
        alc2$size[alc2$chr=='chr11'], alc2$size[alc2$chr=='chr12'],
        alc2$size[alc2$chr=='chr13'], alc2$size[alc2$chr=='chr14'],
        alc2$size[alc2$chr=='chr15'], alc2$size[alc2$chr=='chr16'],
        alc2$size[alc2$chr=='chr17'], alc2$size[alc2$chr=='chr18'],
        alc2$size[alc2$chr=='chr19'], alc2$size[alc2$chr=='chr19'],
        alc2$size[alc2$chr=='chrX'])

boxplot

