#setwd
setwd("/stor/azolla_metagenome/Azolla_genus_metagenome/analyses/")
library(data.table)
library(ggplot2)
header <-  c('assembly', 'sample', 'scaffolded', 'node',    'length',  'coverage', 'classified', 'root',   'life',  'kingdom', 'FCBgroup', 'somegroup', 'phylum', 'class',  'order',  'family', 'genus',  'species')
classes <- c('factor',   'factor', 'factor',     'numeric', 'numeric', 'numeric',  'factor',     'factor', 'factor','factor',  'factor',   'factor',    'factor', 'factor', 'factor', 'factor', 'factor', 'factor')
metrics <- fread(input = "/stor/azolla_metagenome/Azolla_genus_metagenome/analyses/test.tab",
                 sep = '\t',
                 fill = TRUE,
                 header = F,
                 stringsAsFactors = T,
                 data.table = T,
                 col.names = header,
                 select=seq(1,18,1),
                 colClasses = classes)
rm(header,classes)
summary(metrics)

#subset
metrics_5000 <- metrics[length > 5000]
metrics_2500 <- metrics[length > 2500]

# clean
metrics_5000[1,7:17 := tstrsplit(life,' ')[1]]
head(metrics_5000)

# have a length_dist plot
length_dist <- ggplot(metrics_2500,aes(x=node,y=length,size=coverage,alpha=.01,col=assembly))
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
#double filtering helps most in improving lenghts of lower length contigs/scaffolds

#clean
metrics_5000$kingdom <- as.character(metrics_5000$kingdom)
metrics_5000$kingdom <- tstrsplit(x = metrics_5000$kingdom,split = ' ')[1]
metrics_5000$kingdom <- as.factor(metrics_5000$kingdom)

metrics_2500$kingdom <- as.character(metrics_2500$kingdom)
metrics_2500$kingdom <- tstrsplit(x = metrics_2500$kingdom,split = ' ')[1]
metrics_2500$kingdom <- as.factor(metrics_2500$kingdom)


metrics_5000$phylum <- as.character(metrics_5000$phylum)
metrics_5000$phylum <- tstrsplit(x = metrics_5000$phylum,split = ' ')[1]
metrics_5000$phylum <- as.factor(metrics_5000$phylum)

metrics_2500$phylum <- as.character(metrics_2500$phylum)
metrics_2500$phylum <- tstrsplit(x = metrics_2500$phylum,split = ' ')[1]
metrics_2500$phylum <- as.factor(metrics_2500$phylum)

metrics_5000$somegroup <- as.character(metrics_5000$somegroup)
metrics_5000$somegroup <- tstrsplit(x = metrics_5000$somegroup,split = ' ')[1]
metrics_5000$somegroup <- as.factor(metrics_5000$somegroup)

metrics_2500$somegroup <- as.character(metrics_2500$somegroup)
metrics_2500$somegroup <- tstrsplit(x = metrics_2500$somegroup,split = ' ')[1]
metrics_2500$somegroup <- as.factor(metrics_2500$somegroup)



# have a length_dist plot
length_dist <- ggplot(metrics_2500,aes(x=node,y=length,size=coverage,alpha=.01,col=kingdom))
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

length_dist <- ggplot(metrics_2500,aes(x=length,y=coverage,size=length,alpha=.001,col=kingdom))
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

summary(metrics_2500$somegroup)
length_dist <- ggplot(metrics_2500,aes(x=length,y=coverage,size=length,alpha=.001,col=phylum))
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

