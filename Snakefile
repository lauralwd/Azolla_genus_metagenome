HOSTCODES=["azca1_SRR6480231", "azca2_SRR6480201", "azfil_SRR6480158", "azfil_SRR6932851", "azmex_SRR6480159", "azmic_SRR6480161", "aznil_SRR6480196", "aznil_SRR6482158", "azrub_SRR6480160"]
DIRECTIONS=["1","2"]

rule fastqc_raw_data:
  input:
    "data/sequencing_genomic/{hostcode}_{PE}.fastq.gz"
  output:
    "analyses/analyses_reads/{hostcode}_{PE}"
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule all:
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
    stderr="logs/Cat_classify_host.stderr"
  shell:
    "CAT contigs -c {input.c} -d {input.db} -t {input.tf} --out_prefix 'host' -n {threads} 2> {log.stderr} > {log.stdout}"
 

