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

rule CAT_prepare:
  output:
    db="references/cat_database",
#    tf="references/cat_taxonomy"
#  threads: 12
  log: "logs/CAT_prepare.log"
  shell:
#    "CAT prepare --fresh --database_folder {output.db} --taxonomy_folder {output.tf} --nproc {threads} > {log}"
    "wget tbb.bio.uu.nl/bastiaan/CAT_prepare -O {output.db}"


