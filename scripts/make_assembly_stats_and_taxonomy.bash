#!/bin/bash
# This script takes tabular CAT taxonomy files derived from both contigs and scaffolds ($files) of
# several assemblies ($assemblytype) based on several sequencing libraries($hostcode).
# with these tabular files, it creates one summary containing stats coded in the SPAdes header
# after assembly, and the taxonomy added to these long sequences by CAT (Contig Annotation Tool)
# This script seamlessly fits in with the snakemake workflow it came with.

# Let's get the variables from the commandline
assemblytype=$1
hostcodes=$2
files=$3
CPU=$4
MEM=$5
output=$6

echo "running assemblytypes $assemblytype for hostcodes $hostcodes"
echo "Considering assembly files: $files"
echo "Using $CPU CPUs and $MEM RAM and writing to $output"

# run a loop over the input variables given to this script by snakemake wildcards
for   a in ${assemblytype[@]}
do    for   h in ${hostcodes[@]}
      do    for   f in ${files[@]}
            do    grep -v '#' ./data/assembly_"$a"/"$h"/CAT_"$h"_"$f"_taxonomy.tab 	\
			| tr '_' "\t" 							\
			| cut -f 2-  							\
			| sort -k1n --parallel $CPU -s $mem 				\
			| sed  "s/^/$a\t$h\t$f\t/g"
	    done
      done
done | cut -f 1,2,3,4,6,8,10,11,14-20 \
	> "$output"
