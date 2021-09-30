#!/bin/bash

# this script collects read counts, assembly sizes and max RAM usage for assemblies

# Fastq sizes by FASTQC reports
samples=( $(find analyses/analyses_reads -maxdepth 1 -mindepth 1 -type d | sed 's,analyses/analyses_reads/,,') )
#echo ${samples[@]}

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
	| sed -E "s/_1\t$host/\t1\t$host/g"					\
	| sed -E "s/_2\t$host/\t2\t$host/g"
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
	| sed -E "s/_1\t$host/\t1\t$host/g"					\
	| sed -E "s/_2\t$host/\t2\t$host/g"
done

# FastQC reports for filtered files
for s in "${samples[@]}"
do  host=$(echo "$s" | cut -d '_' -f 1,2)
    S=$(echo "$s" | sed 's/_1$/\.1/g' | sed 's/_2$/\.2/g' | sed 's/$/.fastq.00.0_0.cor/'  )
    unzip -p ./analyses/analyses_reads_trimmed_filtered/"$s"/"$S"_fastqc.zip 	\
	"$S"_fastqc/fastqc_data.txt 						\
	| grep -E "Total Sequences|Sequence length"				\
	| tr ' ' _								\
	| sed "s/^/$s\t$host\tfiltered\t/g"					\
	| sed -E 's/[0-9]+-//g'							\
	| sed -E "s/_1\t$host/\t1\t$host/g"					\
	| sed -E "s/_2\t$host/\t2\t$host/g"
done

## FastQC reports for double filtered files
for s in "${samples[@]}"
do  host=$(echo "$s" | cut -d '_' -f 1,2)
    S=$(echo "$s" | sed 's/_1$/\.1/g' | sed 's/_2$/\.2/g'  )
    unzip -p ./analyses/analyses_reads_doublefiltered/"$s"/"$S"_fastqc.zip 	\
	"$S"_fastqc/fastqc_data.txt 						\
	| grep -E "Total Sequences|Sequence length"				\
	| tr ' ' _								\
	| sed -E 's/[0-9]+-//g'							\
	| sed "s/^/$s\t$host\tdoublefiltered\t/g"				\
	| sed -E "s/_1\t$host/\t1\t$host/g"					\
	| sed -E "s/_2\t$host/\t2\t$host/g"
done
unset samples


# Get assembly sizes and RAM usage
assemdirs=( data/assembly_*filtered )
#echo ${assemdirs[@]}
for assemdir in "${assemdirs[@]}"
do  #assemdir='data/assembly_singles_hostfiltered'
    samples=( $(find "$assemdir" -maxdepth 1 -mindepth 1 -type d | sed "s,$assemdir/,,g" ) )
    for s in "${samples[@]}"
    do  host=$(echo "$s" | cut -d '_' -f 1,2)
        assemname=$(echo "$assemdir" | rev | cut -f 1 -d '/' | rev)
        cat "$assemdir/$s"/spades.log		\
	| grep -oE '[0-9]+G +/ +[0-9]+G'	\
	| cut -d '/' -f 1			\
	| tr -d 'G'				\
	| sort -n				\
	| tail -n 1				\
	| sed -E "s/^/$s\t$host\t$assemname-RAM\t/g"
        cat "$assemdir/$s"/contigs.fasta	\
	| grep -v '>'				\
	| wc -c					\
	| sed -E "s/^/$s\t$host\t$assemname-size\t/g"


    done
done

