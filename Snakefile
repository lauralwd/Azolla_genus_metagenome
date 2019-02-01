rule fastqc:
  input:
    "{sequencing_folder}/{host}/{hostcode}.{accession}_{PE}.fastq.gz"
  output:
    "analyses_reads/fastqc/{sequencing_folder}/{host}/{hostcode}.{accession}_{PE}"
  threads: 
    1
  shell:
    "fastqc -o {output} {input} "

rule all:
  input:
    "analyses_reads/fastqc/{sequencing_folder}/{host}/{hostcode}.{accession}_{PE}/"
  shell:
    "ls {input}"

