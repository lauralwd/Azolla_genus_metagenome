#study accession: PRJNA430527
# caroliniana
#SRA: SRS2841097
#ENA sample: SAMN08375393
#ENA run: SRR6480231
mkdir ./sequencing_genomic
mkdir ./sequencing_genomic/caroliniana1
mkdir ./sequencing_genomic/caroliniana2
mkdir ./sequencing_genomic/nilotica
mkdir ./sequencing_genomic/mexicana
mkdir ./sequencing_genomic/microphylla
mkdir ./sequencing_genomic/rubra
mkdir ./sequencing_genomic/filiculoides
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480231/SRR6480231_1.fastq.gz -O ./sequencing_genomic/caroliniana1/azca1.SRR6480231_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480231/SRR6480231_2.fastq.gz -O ./sequencing_genomic/caroliniana1/azca1.SRR6480231_2.fastq.gz
fi
# caroliniana 2
#SRA: SRS2841067
#ENA sample: SAMN08375392
#ENA run: SRR6480201
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480201/SRR6480201_1.fastq.gz -O ./sequencing_genomic/caroliniana2/azca2.SRR6480201_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480201/SRR6480201_2.fastq.gz -O ./sequencing_genomic/caroliniana2/azca2.SRR6480201_2.fastq.gz
fi
# nilotica
#SRA: SRS2841062
#ENA sample: SAMN08375391
#ENA run: SRR6480196
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/006/SRR6480196/SRR6480196.fastq.gz -O ./sequencing_genomic/nilotica/aznil.SRR6480196.fastq.gzif
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/006/SRR6480196/SRR6480196_1.fastq.gz -O ./sequencing_genomic/nilotica/aznil.SRR6480196_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/006/SRR6480196/SRR6480196_2.fastq.gz -O ./sequencing_genomic/nilotica/aznil.SRR6480196_2.fastq.gz
fi
#ENA run: SRR6482158
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6482158/SRR6482158_1.fastq.gz -O ./sequencing_genomic/nilotica/aznil.SRR6482158_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6482158/SRR6482158_2.fastq.gz -O ./sequencing_genomic/nilotica/aznil.SRR6482158_2.fastq.gz
fi
# mexicana
#SRA: SRS2841025
#ENA sample: SAMN08375390
#ENA run: SRR6480159
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/009/SRR6480159/SRR6480159_1.fastq.gz -O  ./sequencing_genomic/mexicana/azmex.SRR6480159_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/009/SRR6480159/SRR6480159_2.fastq.gz -O ./sequencing_genomic/mexicana/azmex.SRR6480159_2.fastq.gz
fi
# microphylla
#SRA: SRS2841027
#ENA sample: SAMN08375389
#ENA run: SRR6480161
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480161/SRR6480161.fastq.gz -O  ./sequencing_genomic/microphylla/azmic.SRR6480161.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480161/SRR6480161_1.fastq.gz -O ./sequencing_genomic/microphylla/azmic.SRR6480161_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480161/SRR6480161_2.fastq.gz -O ./sequencing_genomic/microphylla/azmic.SRR6480161_2.fastq.gz
fi
# rubra
#SRA: SRS2841026
#ENA sample: SAMN08375388
#ENA run: SRR6480160
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/000/SRR6480160/SRR6480160_1.fastq.gz -O ./sequencing_genomic/rubra/azrub.SRR6480160_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/000/SRR6480160/SRR6480160_2.fastq.gz -O ./sequencing_genomic/rubra/azrub.SRR6480160_2.fastq.gz
fi
# filiculoides
#SRA: SRS2841024
#ENA sample: SAMN08375372
## short reads
#ENA experiment: SRX3570018
#ENA run: SRR6480158
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6480158/SRR6480158_1.fastq.gz -O ./sequencing_genomic/filiculoides/azfil.SRR6480158_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6480158/SRR6480158_2.fastq.gz -O ./sequencing_genomic/filiculoides/azfil.SRR6480158_2.fastq.gz 
fi
#ENA experiment: SRX3878102
#ENA run: SRR6932851
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR693/001/SRR6932851/SRR6932851_1.fastq.gz -O ./sequencing_genomic/filiculoides/azfil.SRR6932851_1.fastq.gz
fi
if	[ ! -f ]
then	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR693/001/SRR6932851/SRR6932851_2.fastq.gz -O ./sequencing_genomic/filiculoides/azfil.SRR6932851_2.fastq.gz
fi
## long reads
#ENA experiment: SRX4088716
#HTTP: https://www.ebi.ac.uk/ena/data/view/SRX4088716
