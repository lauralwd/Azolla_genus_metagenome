######## Host genome download
if	[ ! -d ./host_genome ]
then	mkdir  ./host_genome
fi
if	[ ! -f ./host_genome/Azolla_filiculoides.genome_v1.2.fasta ]
then	wget ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.genome_v1.2.fasta -O ./host_genome/Azolla_filiculoides.genome_v1.2.fasta
fi

######## NCBI SRA // EBI ENA study accession: PRJNA430527
# caroliniana
#SRA: SRS2841097
#ENA sample: SAMN08375393
#ENA run: SRR6480231
mkdir ./sequencing_genomic 2> /dev/null
if	[ ! -f ./sequencing_genomic/azca1_SRR6480231_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480231/SRR6480231_1.fastq.gz -O ./sequencing_genomic/azca1_SRR6480231_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/azca1_SRR6480231_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480231/SRR6480231_2.fastq.gz -O ./sequencing_genomic/azca1_SRR6480231_2.fastq.gz
fi
# caroliniana 2
#SRA: SRS2841067
#ENA sample: SAMN08375392
#ENA run: SRR6480201
if	[ ! -f ./sequencing_genomic/azca2_SRR6480201_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480201/SRR6480201_1.fastq.gz -O ./sequencing_genomic/azca2_SRR6480201_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/azca2_SRR6480201_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480201/SRR6480201_2.fastq.gz -O ./sequencing_genomic/azca2_SRR6480201_2.fastq.gz
fi
# nilotica
#SRA: SRS2841062
#ENA sample: SAMN08375391
#ENA run: SRR6480196
if	[ ! -f ./sequencing_genomic/aznil_SRR6480196_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/006/SRR6480196/SRR6480196_1.fastq.gz -O ./sequencing_genomic/aznil_SRR6480196_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/aznil_SRR6480196_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/006/SRR6480196/SRR6480196_2.fastq.gz -O ./sequencing_genomic/aznil_SRR6480196_2.fastq.gz
fi
#ENA run: SRR6482158
if	[ ! -f ./sequencing_genomic/aznil_SRR6482158_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6482158/SRR6482158_1.fastq.gz -O ./sequencing_genomic/aznil_SRR6482158_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/aznil_SRR6482158_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6482158/SRR6482158_2.fastq.gz -O ./sequencing_genomic/aznil_SRR6482158_2.fastq.gz
fi
# mexicana
#SRA: SRS2841025
#ENA sample: SAMN08375390
#ENA run: SRR6480159
if	[ ! -f ./sequencing_genomic/azmex_SRR6480159_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/009/SRR6480159/SRR6480159_1.fastq.gz -O  ./sequencing_genomic/azmex_SRR6480159_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/azmex_SRR6480159_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/009/SRR6480159/SRR6480159_2.fastq.gz -O ./sequencing_genomic/azmex_SRR6480159_2.fastq.gz
fi
# microphylla
#SRA: SRS2841027
#ENA sample: SAMN08375389
#ENA run: SRR6480161
if	[ ! -f ./sequencing_genomic/azmic_SRR6480161_1.fastq.g ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480161/SRR6480161_1.fastq.gz -O ./sequencing_genomic/azmic_SRR6480161_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/azmic_SRR6480161_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480161/SRR6480161_2.fastq.gz -O ./sequencing_genomic/azmic_SRR6480161_2.fastq.gz
fi
# rubra
#SRA: SRS2841026
#ENA sample: SAMN08375388
#ENA run: SRR6480160
if	[ ! -f ./sequencing_genomic/azrub_SRR6480160_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/000/SRR6480160/SRR6480160_1.fastq.gz -O ./sequencing_genomic/azrub_SRR6480160_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/azrub_SRR6480160_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/000/SRR6480160/SRR6480160_2.fastq.gz -O ./sequencing_genomic/azrub_SRR6480160_2.fastq.gz
fi
# filiculoides
#SRA: SRS2841024
#ENA sample: SAMN08375372
## short reads
#ENA experiment: SRX3570018
#ENA run: SRR6480158
if	[ ! -f ./sequencing_genomic/azfil_SRR6480158_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6480158/SRR6480158_1.fastq.gz -O ./sequencing_genomic/azfil_SRR6480158_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/azfil_SRR6480158_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6480158/SRR6480158_2.fastq.gz -O ./sequencing_genomic/azfil_SRR6480158_2.fastq.gz
fi
#ENA experiment: SRX3878102
#ENA run: SRR6932851
if	[ ! -f ./sequencing_genomic/azfil_SRR6932851_1.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR693/001/SRR6932851/SRR6932851_1.fastq.gz -O ./sequencing_genomic/azfil_SRR6932851_1.fastq.gz
fi
if	[ ! -f ./sequencing_genomic/azfil_SRR6932851_2.fastq.gz ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR693/001/SRR6932851/SRR6932851_2.fastq.gz -O ./sequencing_genomic/azfil_SRR6932851_2.fastq.gz
fi
## long reads
#ENA experiment: SRX4088716
#HTTP: https://www.ebi.ac.uk/ena/data/view/SRX4088716
