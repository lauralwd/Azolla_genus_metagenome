#setwd
setwd("/stor/azolla_metagenome/Azolla_genus_metagenome/")
library(data.table)
library(pheatmap)
library(ggplot2)

dm <- fread(input = "./data/assembly_singles_hostfiltered_binningsignals/azca1_SRR6480231/azca1_SRR6480231.hostfiltered.depthmatrix",
            sep = '\t',
            #fill = TRUE,
            header = T,
            stringsAsFactors = F,
            data.table = T,
            #col.names = header,
            select=c(1,2,3,4,6,8,10,12,14,16,18,20,22,24)
            #colClasses = classes
            )
summary(dm)

# cool that's the data. There is a lot of very short contigs, let's get those out before we continue.
dm5000 <- dm[contigLen > 5000]
summary(dm5000)
# drop non-depthmatrix columns
mat <- dm5000[,4:length(dm5000)]
# tranform to log10 for better visualisation
mat <- log10(mat+1)
# Differential abundance is better in log10
pheatmap::pheatmap(mat = mat)

# melt dataframe for convenient plotting in ggplot
dm_melt <- melt(data = dm5000,
                id.vars = colnames(dm5000)[1:3],
                measure.vars = colnames(dm5000)[4:length(dm5000)],
                value.name = 'depth',
                variable.name = 'binsignal',
                variable.factor = T
                )
summary(dm_melt)

# now plot
freqplot <- ggplot(dm_melt,aes(x=depth,col=binsignal))
freqplot <- freqplot + geom_freqpoly()
freqplot <- freqplot + scale_y_log10()
freqplot <- freqplot + scale_x_log10()
freqplot <- freqplot + theme_classic()
freqplot <- freqplot + labs(y='frequency',x='depth of coverage per binsignal')
freqplot

# now filter for everything lower than 10
dm_melt_subset <- dm_melt[depth < 10]
dm_low <- data.table::dcast(data = dm_melt_subset, contigName ~ binsignal,value.var= 'depth')
rm(dm_melt_subset)
summary(dm_low)
# prep for heatmap
mat <- dm_low[,2:length(dm_low)]
mat <- log10(mat+1)
pheatmap::pheatmap(mat = mat)

# do the same for high dpeth contigs
dm_melt_subset <- dm_melt[depth > 10]
dm_high <- data.table::dcast(data = dm_melt_subset, contigName ~ binsignal,value.var= 'depth')
rm(dm_melt_subset)
summary(dm_high)
# prep for heatmap
mat <- dm_high[,2:length(dm_high)]
mat <- log10(mat+1)
pheatmap::pheatmap(mat = mat)

