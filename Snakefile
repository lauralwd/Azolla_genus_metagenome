HOSTCODES=["azca1_SRR6480231", "azca2_SRR6480201", "azfil_SRR6480158", "azfil_SRR6932851", "azmex_SRR6480159", "azmic_SRR6480161", "aznil_SRR6480196", "aznil_SRR6482158", "azrub_SRR6480160"]
DIRECTIONS=["1","2"]

rule fastqc_raw_data:
  input:
    "data/sequencing_genomic/{hostcode}_{PE}.fastq.gz"
  output:
    "analyses/analyses_reads/{hostcode}_{PE}"
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule allfastqc:
  input:
    expand("analyses_reads/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)

rule CAT_download:
  output:
    db="references/CAT_prepare_20190108/2019-01-08_CAT_database",
    tf="references/CAT_prepare_20190108/2019-01-08_taxonomy"
  shell:
    "cd ./references && wget -qO - http://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20190108.tar.gz | tar -xz "

rule CAT_classify_host:
  input:
    c="references/host_genome/Azolla_filiculoides.genome_v1.2.fasta",
    db="references/CAT_prepare_20190108/2019-01-08_CAT_database",
    tf="references/CAT_prepare_20190108/2019-01-08_taxonomy"
  output: 
    "references/host_genome/CAT/"
  threads: 12
  log: 
    stdout="logs/CAT_classify_host.stdout",
    stderr="logs/CAT_classify_host.stderr"
  shell:
    "CAT contigs -c {input.c} -d {input.db} -t {input.tf} --out_prefix 'host' -n {threads} 2> {log.stderr} > {log.stdout} && mv host.* references/host_genome/CAT/"
# shall I remove host.alignment.diamond, it's huge and does not serve a purpose further on.
rule CAT_add_names:
  input:
    i="references/host_genome/CAT/host.contig2classification.txt",
    t="references/CAT_prepare_20190108/2019-01-08_taxonomy"
  output:
     "references/host_genome/contig_taxonomy.tab"
  log: 
    stdout="logs/CAT_addnames_host.stdout",
    stderr="logs/CAT_addnames_host.stderr"
  threads: 1
  shell:
    "CAT add_names -i {input.i} -t {input.t} -o {output} > {log.stdout} 2> {log.stderr}"

rule CAT_filter_contignames:
  input:
     "references/host_genome/contig_taxonomy.tab"
  output:
     "references/host_genome/contig_filterlist.txt"
  threads: 1
  log:
    "logs/CAT_contigfilterlist.stderr"
  shell:
    "cat {input} | grep -v Bacteria | grep -v Fungi | grep -v Opisthokonta | grep -v Alveolata | grep Eukaryota | cut -f 1  > {output} 2> {log}"

rule create_host_filter_fasta:
  input:
    n="references/host_genome/contig_filterlist.txt",
    f="references/host_genome/Azolla_filiculoides.genome_v1.2.fasta"
  output:
    "references/host_genome/host_filter.fasta"
  threads: 1
  log: 
    stdout="logs/CAT_createhostfilter.stdout",
    stderr="logs/CAT_createhostfilter.stderr"
  shell:
    "samtools faidx {input.f} -o {output} -r {input.n} > {log.stdout} 2> {log.stderr}"

rule create_host_filter_bt2_index:
  input:
    "references/host_genome/host_filter.fasta"
  output:
    output1="references/host_genome/host_filter_bt2index/host_filter.1.bt2",
    output2="references/host_genome/host_filter_bt2index/host_filter.2.bt2",
    output3="references/host_genome/host_filter_bt2index/host_filter.3.bt2",
    output4="references/host_genome/host_filter_bt2index/host_filter.4.bt2",
    outputrev1="references/host_genome/host_filter_bt2index/host_filter.rev1.bt2",
    outputrev2="references/host_genome/host_filter_bt2index/host_filter.rev2.bt2",
    base="references/host_genome/host_filter_bt2index/host_filter"
  threads: 12
  log:
    stdout="logs/CAT_createhostfilterbt2index.stdout",
    stderr="logs/CAT_createhostfilterbt2index.stderr"
  shell:
    "bowtie2-build --threads {threads} {input} {output.base} > {log.stdout} 2> {log.stderr}"

rule trimmomatic_genomic_sequencing:
  input:
    expand("data/sequencing_genomic/{{hostcode}}_{PE}.fastq.gz", PE=DIRECTIONS)
  params:
    outputbase="data/sequencing_genomic_trimmed/{hostcode}"
  output:
    "p1={params}_1.fastq.gz",
    "p2={params}_2.fastq.gz",
    "s1={params}_s1.fastq.gz",
    "s2={params}_s2.fastq.gz"
  threads: 2
  log:
    log="logs/trimmomatic_genomicsequencing.log",
    summary="logs/trimmomatic_genomicsequencing_summary.log"
  shell:
    "trimmomatic PE -threads {threads} -trimlog {log.log} -summary {log.summary} {input} {output.p1} {output.s1} {output.p2} {output.s2}"

rule all_rimmed:
  input:
    expand("data/sequencing_genomic_trimmed/{hostcode}_1.fastq.gz", hostcode=HOSTCODES)

rule filter_for_host:
  input:
    s1="data/sequencing_genomic_trimmed/{hostcode}_1.fastq.gz",
    s2="data/sequencing_genomic_trimmed/{hostcode}_2.fastq.gz",
    i="references/host_genome/host_filter_bt2index"
  output:
    "data/sequencing_genomic_trimmed_mapped/{hostcode}_{PE}.bam"
  threads: 12
  log:
    stdout="logs/bowtie2filterforhost.stdout",
    stderr="logs/bowtie2filterforhost.stderr"
  shell:
    "bowtie2 -x {input.i} -1 {input.s1} -2 {input.s2}   > {log.stdout} 2> {log.stderr}"
