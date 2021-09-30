#!/bin/bash

# this script collects read counts, assembly sizes and max RAM usage for assemblies

# Fastq sizes by FASTQC reports
samples=( $(find analyses/analyses_reads -maxdepth 1 -type d | sed 's,analyses/analyses_reads/,,') )


echo -e "sample\tdirection\thost\tstage\tmetric\tvalue"
# FastQC reports for genomic input, or "raw" files
for s in "${samples[@]}"
do  host=$(echo "$s" | cut -d '_' -f 1,2)
    unzip -p ./analyses/analyses_reads/"$s"/"$s"_fastqc.zip 			\
	"$s"_fastqc/fastqc_data.txt 						\
	| grep -E "Total Sequences|Sequence length"				\
	| tr ' ' _								\
	| sed -E 's/[0-9]+-//g'							\
	| sed "s/^/$s\t$host\traw\t/g"						\
	| sed -E "s/_1\t/\t1\t/g"						\
	| sed -E "s/_2\t/\t2\t/g"
done

# FastQC reports for genomic trimmed files
for s in "${samples[@]}"
do  host=$(echo "$s" | cut -d '_' -f 1,2)
    unzip -p ./analyses/analyses_reads_trimmed/"$s"/"$s"_fastqc.zip 		\
	"$s"_fastqc/fastqc_data.txt 						\
	| grep -E "Total Sequences|Sequence length"				\
	| tr ' ' _								\
	| sed -E 's/[0-9]+-//g'							\
	| sed "s/^/$s\t$host\ttrimmed\t/g"					\
	| sed -E "s/_1\t/\t1\t/g"						\
	| sed -E "s/_2\t/\t2\t/g"
done

# FastQC reports for filtered files
for s in "${samples[@]}"
do  host=$(echo "$s" | cut -d '_' -f 1,2)
    unzip -p ./analyses/analyses_reads_trimmed_filtered/"$s"/"$s"_fastqc.zip 	\
	"$s"_fastqc/fastqc_data.txt 						\
	| grep -E "Total Sequences|Sequence length"				\
	| tr ' ' _								\
	| sed -E 's/[0-9]+-//g'							\
	| sed "s/^/$s\t$host\tfiltered\t/g"					\
	| sed -E "s/_1\t/\t1\t/g"						\
	| sed -E "s/_2\t/\t2\t/g"
done

# FastQC reports for double filtered files
for s in "${samples[@]}"
do  host=$(echo "$s" | cut -d '_' -f 1,2)
    unzip -p ./analyses/analyses_reads_doublefiltered/"$s"/"$s"_fastqc.zip 	\
	"$s"_fastqc/fastqc_data.txt 						\
	| grep -E "Total Sequences|Sequence length"				\
	| tr ' ' _								\
	| sed -E 's/[0-9]+-//g'							\
	| sed "s/^/$s\t$host\ttrimmed\t/g"					\
	| sed -E "s/_1\t/\t1\t/g"						\
	| sed -E "s/_2\t/\t2\t/g"
done


