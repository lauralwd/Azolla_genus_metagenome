# Azolla genus metagenome
***This project aims do describe and explore the metagenomes of fern species of the genus Azolla.***

## Current workflow
Current workflow specifically filters plant DNA out of a specific set of Azolla genomic sequencing runs to assembly microbial genomes scatered within this dataset. 
The ambition is to generalise the code so it can take any other set of sequencing files and look within these for microbial genomes as well. 
Feel free to fork or contribute, I would be happy to work together.

A consise visual summary of the snakemake workflow would look like this:

![Snakemake rule graph](https://github.com/lauralwd/Azolla_genus_metagenome/blob/binning/rulegraph.svg)

## Students
Students willing to travel to-, or currently studying at Utrecht University are most welcome to discuss starting an internship for ECTS as part of this project. 
Naturally, at least some skills in coding and knowledge of metagenomics are required. 
Please find my contact details on my [personel page](https://www.uu.nl/medewerkers/LWDijkhuizen).

## Biological context
Floating ferns of the genus _Azolla_ are known for their productive symbiosis with cyanobacterium _Nostoc azollae_. 
Residing inside specialised _Azolla sp._ leaf pockets, _N. azollae_ fixes N<sub>2</sub> using energy of its own photosynthesis thereby providing the host fern with sufficient nitrogen to maintain optimum growth rates in absence of any nitrogen fertilisation. 
The symbiosis is unique to the _Azolla_ genus dating back to 90M years [(Metzgar,2007)](http://doi.org/doi.org/10.1086/519007); the only other genus within this family of ferns lacks the symbiosis. 
The symbiotic cyanobacteria are systematically transfered through generations of ferns via their megaspores and have eroded genomes; Ran et al. (2010) therefore argued that the boundary between symbiotic partner and plant organelle is thin or even absent in the case of Azolla. 
Co-speciation patterns of the fern (chloroplast) and _N. azollae_ support the close association between the symbiotic partners (Li et al 2018). While early research indicated that _Azolla sp._ has only one symbiont, electron microscopy [(Carrapico 1990)](http://doi.org/10.1007/BF02187448) revealed bacterial cells accompany the cyanobacteria both in the hosts leaf pockets and megaspores. 
My own study [(Dijkhuizen 2018)](https://doi.org/10.1111/nph.14843) describes one of several bacterial genomes which were found as an efficient by-product of the _Azolla_ genome project and confirmed in plants taken from the environment. 
The paper shows that similar bactera are likely present in other _Azolla_ species as well albeit in lower relative abundances than the cyanobacteria. 
Given 
(1) the transfer mechanism of (cyano)bacteria over _Azolla_ generations, 
(2) the high level of specialisation of the Cyanobacteria and the host leaf pocket, and 
(3) the reproducible presence of other bacteria of similar genomic content in the _Azolla_ leaf pockets, 
I theorise that these third parties in the symbiosis have some fittness benefit to the plant-microbe consortium.
Previous work and past student projects suggest that these microbes are endophytic (i.e. living inside the plant) and may share an evolutionary origin (Dijkhuizen 2018). 
A metagenomic approach like the one that is under construction here, may answer questions of similarity of microbes in genome sequence or pathways encoded. 
Some microbes may be passengers, some may be persistent endophytes. 
This approach may shed light on the evolutionary history of these microbes and if these are shared or not.
Finally it could indicate conserved types of functions encoded in these genomes.
Their function however remains elusive.

## Technical context
Recently, a wealth of data was made available to the community including WGS sequencing of multiple _Azolla_ species, the genome of _Azolla filiculoides_, and genomes of several strains of _N. azollae_; strains which are symbionts to different species of host plant (Li et al. 2018). 
Here I aim to repurpose this data to assemble the individual metagenomes of the host plants, specifically those hypothetical third partners in the symbiosis. 
Figures and results may be included at a later stage when confidence of their correctness allows so. 
These metagenomes may elucidate what, if any, microbes from the different host plants share similarity, origin, and/or functions.

### Reproducibility
This study aims to achieve full reproducibility by using the [anaconda](https://anaconda.org/) and [snakemake](https://snakemake.readthedocs.io/en/stable/) frameworks for reproducible science. 
All versions of the workflow will be version controlled, here on github. 
Data used in the study is already freely available in the Sequencing Read Archive. 
Accesion numbers are detailed in the respository.
