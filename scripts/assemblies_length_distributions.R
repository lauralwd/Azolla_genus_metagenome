#setwd
setwd("/stor/azolla_metagenome/Azolla_genus_metagenome/analyses/")
library(data.table)
library(ggplot2)
header <- c('assembly','sample','node','length','coverage','root',  'classified', 'life','kingdom','FCGgroup','somegroup','phylum','class','order','family','genus','species')
classes <- c('factor','factor','numeric','numeric','numeric','factor',  'factor', 'factor','factor','factor','factor','factor','factor','factor','factor','factor','factor')
contigs <- fread(input = "/stor/azolla_metagenome/Azolla_genus_metagenome/analyses/assembly_basicstats_and_CAT.tab",
                 sep = '\t',
#                 nrows = 100,
                 fill = T,
                 header = F,
                 stringsAsFactors = T,
                 data.table = T,
                 col.names = header,
                 colClasses = classes)
rm(header,classes)
head(contigs)

#subset

contigs_5000 <- contigs[length > 5000]
contigs_2500 <- contigs[length > 2500]

# clean
contigs_5000[1,7:17 := tstrsplit(life,' ')[1]]
head(contigs_5000)

# have a length_dist plot
length_dist <- ggplot(contigs_2500,aes(x=node,y=length,size=coverage,alpha=.01,col=assembly))
length_dist <- length_dist + facet_grid(~ sample)
length_dist <- length_dist + geom_point()
length_dist <- length_dist + scale_y_log10()
length_dist <- length_dist + labs(x='contig length rank',y='contig length')
length_dist <- length_dist + theme_classic()
length_dist
ggsave(filename = "./length_distributions_per_assembly_overlay.png",
       plot = length_dist,
       device = 'png',
       dpi = 500,
       units = 'cm',
       width = 42,
       height = 29.7 )
rm(length_dist)
#double filtering helps most in improving lenghts of lower length contigs.

#clean
contigs_5000$kingdom <- as.character(contigs_5000$kingdom)
contigs_5000$kingdom <- tstrsplit(x = contigs_5000$kingdom,split = ' ')[1]
contigs_5000$kingdom <- as.factor(contigs_5000$kingdom)

contigs_2500$kingdom <- as.character(contigs_2500$kingdom)
contigs_2500$kingdom <- tstrsplit(x = contigs_2500$kingdom,split = ' ')[1]
contigs_2500$kingdom <- as.factor(contigs_2500$kingdom)


contigs_5000$phylum <- as.character(contigs_5000$phylum)
contigs_5000$phylum <- tstrsplit(x = contigs_5000$phylum,split = ' ')[1]
contigs_5000$phylum <- as.factor(contigs_5000$phylum)

contigs_2500$phylum <- as.character(contigs_2500$phylum)
contigs_2500$phylum <- tstrsplit(x = contigs_2500$phylum,split = ' ')[1]
contigs_2500$phylum <- as.factor(contigs_2500$phylum)

contigs_5000$somegroup <- as.character(contigs_5000$somegroup)
contigs_5000$somegroup <- tstrsplit(x = contigs_5000$somegroup,split = ' ')[1]
contigs_5000$somegroup <- as.factor(contigs_5000$somegroup)

contigs_2500$somegroup <- as.character(contigs_2500$somegroup)
contigs_2500$somegroup <- tstrsplit(x = contigs_2500$somegroup,split = ' ')[1]
contigs_2500$somegroup <- as.factor(contigs_2500$somegroup)



# have a length_dist plot
length_dist <- ggplot(contigs_2500,aes(x=node,y=length,size=coverage,alpha=.01,col=kingdom))
length_dist <- length_dist + facet_grid(assembly ~ sample)
length_dist <- length_dist + geom_point()
length_dist <- length_dist + scale_y_log10()
length_dist <- length_dist + labs(x='contig length rank',y='contig length')
length_dist <- length_dist + theme_classic()
length_dist
ggsave(filename = "./length_distributions_per_assembly_kingdom.png",
       plot = length_dist,
       device = 'png',
       dpi = 500,
       units = 'cm',
       width = 42,
       height = 29.7 )

length_dist <- ggplot(contigs_2500,aes(x=length,y=coverage,size=length,alpha=.001,col=kingdom))
length_dist <- length_dist + facet_grid(assembly ~ sample,scales = 'free_x')
length_dist <- length_dist + geom_point()
length_dist <- length_dist + scale_x_log10()
length_dist <- length_dist + scale_y_log10()
length_dist <- length_dist + labs(x='contig length',y='contig coverage')
length_dist <- length_dist + theme_classic()
length_dist
ggsave(filename = "./coverage_overlength_distributions_per_assembly_kingdom.png",
       plot = length_dist,
       device = 'png',
       dpi = 500,
       units = 'cm',
       width = 42,
       height = 29.7 )

summary(contigs_2500$somegroup)
length_dist <- ggplot(contigs_2500,aes(x=length,y=coverage,size=length,alpha=.001,col=phylum))
length_dist <- length_dist + facet_grid(assembly ~ sample)
length_dist <- length_dist + geom_point()
length_dist <- length_dist + scale_x_log10()
length_dist <- length_dist + scale_y_log10()
length_dist <- length_dist + labs(x='contig length',y='contig coverage')
length_dist <- length_dist + theme_classic()
length_dist
ggsave(filename = "./coverage_overlength_distributions_per_assembly_order.png",
       plot = length_dist,
       device = 'png',
       dpi = 500,
       units = 'cm',
       width = 42,
       height = 29.7 )

