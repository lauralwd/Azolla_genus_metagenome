HOSTCODES=["azca1_SRR6480231", "azca2_SRR6480201", "azfil_SRR6480158", "azfil_SRR6932851", "azmex_SRR6480159", "azmic_SRR6480161", "aznil_SRR6480196", "aznil_SRR6482158", "azrub_SRR6480160"]
DIRECTIONS=["1","2"]

## 'All'-rules
rule allfastqc:
  input:
    expand("analyses/analyses_reads/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)
rule allfastqc_trimmed:
  input:
    expand("analyses/analyses_reads_trimmed/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)
rule allfiltered:
  input:
    expand("data/sequencing_genomic_trimmed_filtered/{hostcode}.{PE}",hostcode=HOSTCODES,PE=DIRECTIONS)
rule allfirstassemblies:
  input:
    expand("data/hostfiltered_assembly_{hostcode}/contigs.fasta",hostcode=HOSTCODES)
rule allsecondassemblies:
  input:
    expand("data/doublefiltered_assembly_{hostcode}/contigs.fasta",hostcode=HOSTCODES)
    
## analyses rules
rule fastqc_raw_data:
  input:
    "data/sequencing_genomic/{hostcode}_{PE}.fastq.gz"
  output:
    "analyses/analyses_reads/{hostcode}_{PE}"
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule fastqc_trimmed_data:
  input:
    "data/sequencing_genomic_trimmed/{hostcode}_{PE}.fastq.gz"
  output:
    "analyses/analyses_genomic_trimmed/{hostcode}_{PE}"
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

## reference rules
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
  params:
    "references/host_genome/host_filter_bt2index/host_filter"
  output:
    expand("references/host_genome/host_filter_bt2index/host_filter.{i}.bt2",i=range(1,4)),
    expand("references/host_genome/host_filter_bt2index/host_filter.rev.{i}.bt2",i=range(1,2))
  threads: 12
  log:
    stdout="logs/CAT_createhostfilterbt2index.stdout",
    stderr="logs/CAT_createhostfilterbt2index.stderr"
  shell:
    "bowtie2-build --threads {threads} {input} {params} > {log.stdout} 2> {log.stderr}"

## data processing rules
rule trimmomatic_genomic_sequencing:
  input:
    expand("data/sequencing_genomic/{{hostcode}}_{PE}.fastq.gz", PE=DIRECTIONS)
  params:
    "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36"
  output:
    p1="data/sequencing_genomic_trimmed/{hostcode}_1.fastq.gz",
    p2="data/sequencing_genomic_trimmed/{hostcode}_2.fastq.gz",
  threads: 2
  log:
    log="logs/trimmomatic_genomicsequencing.{hostcode}.log",
    summary="logs/trimmomatic_genomicsequencing_summary{hostcode}.log"
  shell:
    "trimmomatic PE -threads {threads} -trimlog {log.log} -summary {log.summary} {input} {output.p1} /dev/null {output.p2} /dev/null {params}"

rule filter_for_host:
  input:
    expand("references/host_genome/host_filter_bt2index/host_filter.{i}.bt2",i=range(1,4)),
    expand("references/host_genome/host_filter_bt2index/host_filter.rev.{i}.bt2",i=range(1,2)),
    s1=expand("data/sequencing_genomic_trimmed/{{hostcode}}_{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed/{{hostcode}}_{PE}.fastq.gz",PE=2)
  params:
    opts="--very-sensitive -N 1",
    i="references/host_genome/host_filter_bt2index/host_filter",
    outbase="data/sequencing_genomic_trimmed_filtered/{hostcode}"
  output:
       expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=DIRECTIONS)
  threads: 12
  log:
    stderr="logs/bowtie2filterforhost{hostcode}.stderr"
  shell:
    "bowtie2 {params.opts} --threads {threads} --un-conc-gz {params.outbase} -x {params.i} -1 {input.s1} -2 {input.s2}   > /dev/null 2> {log.stderr}"

rule spades_first_assembly:
  input:
    reads=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=DIRECTIONS),
    s1=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=1),
    s2=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=2)
  params:
    "--meta"
  output:
    basedir="data/hostfiltered_assembly_{hostcode}",
    contigs="data/hostfiltered_assembly_{hostcode}/contigs.fasta",
    r1="data/hostfiltered_assembly_{hostcode}/corrected/{hostcode}.R1.00.0_0.fastq.gz",
    r2="data/hostfiltered_assembly_{hostcode}/corrected/{hostcode}.R2.00.0_0.fastq.gz"
  threads: 12
  resources:
    mem_mb=450000
  log:
    stdout="logs/SPADES_first_assembly_{hostcode}.stdout",
    stderr="logs/SPADES_first_assembly_{hostcode}.stderr" 
  shell:
    "spades.py {params} -t {threads} -m {resources.mem_mb} -1 {input.s1} -2 {input.s2} -o {output.basedir} > {log.stdout} 2> {log.stderr}"

## process first assembly for second filter
rule CAT_first_spades_assembly:
  input:
    contigs="data/hostfiltered_assembly_{hostcode}/contigs.fasta",
    db="references/CAT_prepare_20190108/2019-01-08_CAT_database",
    tf="references/CAT_prepare_20190108/2019-01-08_taxonomy"
  params:
  output:
    "data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}",
    i="data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}.contig2classification.txt"
  threads: 12
  log:
    stdout="logs/CAT_first_assembly_taxonomy{hostcode}.stdout",
    stderr="logs/CAT_first_assembly_taxonomy{hostcode}.stderr" 
  shell:
    "CAT contigs -c {input.contigs} -d {input.db} -t {input.tf} --out_prefix {output} -n {threads} 2> {log.stderr} > {log.stdout}"

rule CAT_add_names_first_spades_assembly:
  input:
    i="data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}.contig2classification.txt",
    t="references/CAT_prepare_20190108/2019-01-08_taxonomy"
  output:
     "data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_contig_taxonomy.tab"
  log: 
    stdout="logs/CAT_first_assembly_taxonomy_{hostcode}.stdout",
    stderr="logs/CAT_first_assembly_taxonomy_{hostcode}.stderr"
  threads: 1
  shell:
    "CAT add_names -i {input.i} -t {input.t} -o {output} > {log.stdout} 2> {log.stderr}"

rule CAT_filter_contignames_first_spades_assembly:
  input:
     "data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_contig_taxonomy.tab"
  output:
    "data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_contig_filterlist.txt"
  threads: 1
  log:
    "logs/CAT_first_assembly_taxonomy_{hostcode}_contigfilterlist.stderr"
  shell:
    "cat {input} | grep -v Bacteria | grep -v Fungi | grep -v Opisthokonta | grep -v Alveolata | grep Eukaryota | cut -f 1  > {output} 2> {log}"

rule create_filter_fasta_first_spades_assembly:
  input:
    n="data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_contig_filterlist.txt",
    f="data/hostfiltered_assembly_{hostcode}/contigs.fasta"
  output:
    "data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_filter.fasta"
  threads: 1
  log: 
    stdout="logs/CAT_create_{hostcode}_first_filter.stdout",
    stderr="logs/CAT_create_{hostcode}_first_filter.stderr"
  shell:
    "samtools faidx {input.f} -o {output} -r {input.n} > {log.stdout} 2> {log.stderr}"

rule create_first_spades_assembly_filter_bt2_index:
  input:
    "data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_filter.fasta"
  params:
    "data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_filter_bt2index/{hostcode}_filter"
  output:
    expand("data/hostfiltered_assembly_taxonomy_{{hostcode}}/CAT_{{hostcode}}_filter_bt2index/{{hostcode}}_filter.{i}.bt2",i=range(1,4)),
    expand("data/hostfiltered_assembly_taxonomy_{{hostcode}}/CAT_{{hostcode}}_filter_bt2index/{{hostcode}}_filter.rev.{i}.bt2",i=range(1,2)),
  threads: 12
  log:
    stdout="logs/CAT_create_{hostcode}_filterbt2index.stdout",
    stderr="logs/CAT_create_{hostcode}_filterbt2index.stderr"
  shell:
    "bowtie2-build --threads {threads} {input} {params} > {log.stdout} 2> {log.stderr}"

rule all:
  input:
    expand("data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_filter.fasta",hostcode=HOSTCODES)

rule filter_for_assembly:
  input:
    expand("data/hostfiltered_assembly_taxonomy_{{hostcode}}/CAT_{{hostcode}}_filter_bt2index/{{hostcode}}_filter.{i}.bt2",i=range(1,4)),
    expand("data/hostfiltered_assembly_taxonomy_{{hostcode}}/CAT_{{hostcode}}_filter_bt2index/{{hostcode}}_filter.rev.{i}.bt2",i=range(1,2)),
    s1="data/hostfiltered_assembly_{hostcode}/corrected/{hostcode}.R1.00.0_0.fastq.gz",
    s2="data/hostfiltered_assembly_{hostcode}/corrected/{hostcode}.R2.00.0_0.fastq.gz"
  params:
    opts="--very-fast",
    i="data/hostfiltered_assembly_taxonomy_{hostcode}/CAT_{hostcode}_filter_bt2index/{hostcode}_filter",
    outbase="data/assemblyfiltered_reads/{hostcode}"
  output:
       expand("data/assemblyfiltered_reads/{{hostcode}}.{PE}",PE=DIRECTIONS)
  threads: 12
  log:
    stderr="logs/bowtie2filterforfirstassembly{hostcode}.stderr"
  shell:
    "bowtie2 {params.opts} --threads {threads} --un-conc-gz {params.outbase} -x {params.i} -1 {input.s1} -2 {input.s2}   > /dev/null 2> {log.stderr}"

rule spades_second_assembly:
  input:
    reads=expand("data/assemblyfiltered_reads/{{hostcode}}.{PE}",PE=DIRECTIONS),
    s1=expand("data/assemblyfiltered_reads/{{hostcode}}.{PE}",PE=1),
    s2=expand("data/assemblyfiltered_reads/{{hostcode}}.{PE}",PE=2)
  params:
    "--meta --no-correction"
  output:
    basedir="data/doublefiltered_assembly_{hostcode}",
    contigs="data/doublefiltered_assembly_{hostcode}/contigs.fasta"
  threads: 12
  resources:
    mem_mb=450000
  log:
    stdout="logs/SPADES_second_assembly_{hostcode}.stdout",
    stderr="logs/SPADES_second_assembly_{hostcode}.stderr" 
  shell:
    "spades.py {params} -t {threads} -m {resources.mem_mb} -1 {input.s1} -2 {input.s2} -o {output.basedir} > {log.stdout} 2> {log.stderr}"
