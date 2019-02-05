HOSTCODES=["azca1_SRR6480231", "azca2_SRR6480201", "azfil_SRR6480158", "azfil_SRR6932851", "azmex_SRR6480159", "azmic_SRR6480161", "aznil_SRR6480196", "aznil_SRR6482158", "azrub_SRR6480160"]
DIRECTIONS=["1","2"]

rule fastqc_raw_data:
  input:
    "sequencing_genomic/{hostcode}_{PE}.fastq.gz"
  output:
    "analyses_reads/{hostcode}_{PE}"
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule all:
  input:
    expand("analyses_reads/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)
