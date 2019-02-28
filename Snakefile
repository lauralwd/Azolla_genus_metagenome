#HOSTCODES=["azca1_SRR6480231", "azca2_SRR6480201", "azfil_SRR6480158", "azfil_SRR6932851", "azmex_SRR6480159", "azmic_SRR6480161", "aznil_SRR6480196", "aznil_SRR6482158", "azrub_SRR6480160"]
HOSTCODES= ['Azfil.lab.250', 'Azfil.lab.500', 'Azfil.lab.800', 'Azfil.minuscyano.170', 'Azfil.minuscyano.350', 'Azfil.minuscyano_HIC', 'Azfil.wild.galgw_E_1', 'Azfil.wild.galgw_E_2', 'Azfil.wild.galgw_E_3', 'Azfil.wild.galgw_P_2', 'Azfil.wild.galgw_P_3', 'Azfil.wild.galgw_P_4', 'Azmex.IRRI.486', 'Azmic.IRRI.456', 'Aznil.IRRI.479', 'Azrub.IRRI.479', 'Azspnov.IRRI_1.472', 'Azspnov.IRRI_2.489'] 
DIRECTIONS=["1","2"]

ASSEMBLYTYPES=['singles_doublefiltered','singles_hostfiltered'] # ,'hybrid_doublefiltered']
BINNINGSIGNALS=['dijkhuizen2018.E.1', 'dijkhuizen2018.E.2', 'dijkhuizen2018.E.3', 'dijkhuizen2018.P.2', 'dijkhuizen2018.P.3', 'dijkhuizen2018.P.4','ran2010.nostoc.SRR066216','ran2010.nostoc.SRR066217','ran2010.nostoc.SRR3923641','ran2010.nostoc.SRR3923645','ran2010.nostoc.SRR3923646']

## 'All'-rules
rule all:
  input:
    expand("analyses/analyses_reads/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS),
    expand("data/assembly_singles_doublefiltered/{hostcode}/contigs.fasta",hostcode=HOSTCODES),
    expand("analyses/analyses_reads_trimmed/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)

rule allfastqc:
  input:
    expand("analyses/analyses_reads/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)
rule allfastqc_trimmed:
  input:
    expand("analyses/analyses_reads_trimmed/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)
rule allfiltered:
  input:
    expand("data/sequencing_genomic_trimmed_filtered/{hostcode}.{PE}.fastq.gz",hostcode=HOSTCODES,PE=DIRECTIONS)
rule allfirstassemblies:
  input:
    expand("data/assembly_singles_hostfiltered/{hostcode}/contigs.fasta",hostcode=HOSTCODES)
rule allreadscorrected:
  input:
    expand("data/sequencing_genomic_trimmed_filtered_corrected/{hostcode}/corrected/{hostcode}.{PE}.fastq.00.0_0.cor.fastq.gz",hostcode=HOSTCODES, PE=DIRECTIONS)
rule allsecondassemblies:
  input:
    expand("data/assembly_singles_doublefiltered/{hostcode}/contigs.fasta",hostcode=HOSTCODES)

## analyses rules
rule fastqc_raw_data:
  input:
    "data/sequencing_genomic/{hostcode}.{PE}.fastq.gz"
  output:
    "analyses/analyses_reads/{hostcode}_{PE}"
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule fastqc_trimmed_data:
  input:
    "data/sequencing_genomic_trimmed/{hostcode}_{PE}.fastq.gz"
  output:
    "analyses/analyses_reads_trimmed/{hostcode}_{PE}"
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
  threads: 100
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
  threads: 100
  log:
    stdout="logs/CAT_createhostfilterbt2index.stdout",
    stderr="logs/CAT_createhostfilterbt2index.stderr"
  shell:
    "bowtie2-build --threads {threads} {input} {params} > {log.stdout} 2> {log.stderr}"

## data processing rules
rule trimmomatic_genomic_sequencing:
  input:
    expand("data/sequencing_genomic/{{hostcode}}.{PE}.fastq.gz", PE=DIRECTIONS)
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
    opts="--very-fast",
    i="references/host_genome/host_filter_bt2index/host_filter",
    outbase="data/sequencing_genomic_trimmed_filtered/{hostcode}"
  output:
       expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=DIRECTIONS)
  threads: 100
  log:
    stderr="logs/bowtie2filterforhost{hostcode}.stderr"
  shell:
    "bowtie2 {params.opts} --threads {threads} --un-conc-gz {params.outbase} -x {params.i} -1 {input.s1} -2 {input.s2}   > /dev/null 2> {log.stderr}"

rule filter_for_host_rename:
  input:
    reads=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=DIRECTIONS),
    s1=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=1),
    s2=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=2)
  output:
    reads=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}.fastq.gz",PE=DIRECTIONS),
    s1=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}.fastq.gz",PE=2)
  shell:
    "mv {input.s1} {output.s1} && mv {input.s2} {output.s2}"

rule spades_hammer:
  input:
    reads=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}.fastq.gz",PE=DIRECTIONS),
    s1=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}.fastq.gz",PE=2)
  params:
    "--only-error-correction"
  output:
    basedir="data/sequencing_genomic_trimmed_filtered_corrected/{hostcode}/",
    reads=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=DIRECTIONS)
  threads: 100
  log:
    stdout="logs/SPAdes_correct_sequencing{hostcode}.stdout",
    stderr="logs/SPAdes_correct_sequencing{hostcode}.stderr"
  shell:
    "spades.py {params} -t {threads} -1 {input.s1} -2 {input.s2} -o {output.basedir} > {log.stdout} 2> {log.stderr}"

rule spades_first_assembly:
  input:
    reads=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=DIRECTIONS),
    s1=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=2)
  params:
    "--meta --only-assembler"
  output:
    basedir=expand("data/assembly_{assemblytype}/{{hostcode}}/",assemblytype='singles_hostfiltered'),
    contigs=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs.fasta",assemblytype='singles_hostfiltered'),
    scaffolds=expand("data/assembly_{assemblytype}/{{hostcode}}/scaffolds.fasta",assemblytype='singles_hostfiltered')
  threads: 100
  resources:
    mem_mb=500
  log:
    stdout=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stdout",assemblytype='singles_hostfiltered'),
    stderr=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stderr",assemblytype='singles_hostfiltered')
  shell:
    "spades.py {params} -t {threads} -m {resources.mem_mb} -1 {input.s1} -2 {input.s2} -o {output.basedir} > {log.stdout} 2> {log.stderr}"

## process first assembly for second filter
rule CAT_prepare_ORFS:
  input:
    contigs=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs.fasta",assemblytype='singles_hostfiltered')
  output: 
    p=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs_predicted_proteins.fasta",assemblytype='singles_hostfiltered'),
    g=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs_predicted_proteins.gff",assemblytype='singles_hostfiltered')
  params:
    "-p meta -g 11 -q -f gff"
  threads: 1
  log:
    stdout=expand("logs/CAT_assembly_{assemblytype}_prodigal_{{hostcode}}.stdout",assemblytype='singles_hostfiltered'),
    stderr=expand("logs/CAT_assembly_{assemblytype}_prodigal_{{hostcode}}.stderr",assemblytype='singles_hostfiltered')
  shell:
    "prodigal -i {input.contigs} -a {output.p} -o {output.g} {params}"

rule CAT_first_spades_assembly:
  input:
    contigs=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs.fasta",assemblytype='singles_hostfiltered'),
    db="references/CAT_prepare_20190108/2019-01-08_CAT_database",
    tf="references/CAT_prepare_20190108/2019-01-08_taxonomy",
    p=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs_predicted_proteins.fasta",assemblytype='singles_hostfiltered')
  output:
    i=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.contig2classification.txt",assemblytype='singles_hostfiltered'),
#    g=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.predicted_proteins.gff",assemblytype='singles_hostfiltered'),
#    f=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.predicted_proteins.faa",assemblytype='singles_hostfiltered'),
    o=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.ORF2LCA.txt",assemblytype='singles_hostfiltered'),
    l=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.log",assemblytype='singles_hostfiltered')
  shadow: "shallow"
  params: b = lambda w: expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}",assemblytype='singles_hostfiltered',hostcode=w.hostcode)
  threads: 100
  resources:
    mem_mb=30000
  log:
    stdout=expand("logs/CAT_assembly_{assemblytype}_classification_{{hostcode}}.stdout",assemblytype='singles_hostfiltered'),
    stderr=expand("logs/CAT_assembly_{assemblytype}_classification_{{hostcode}}.stderr",assemblytype='singles_hostfiltered')
  shell:
    "CAT contigs -c {input.contigs} -d {input.db} -p {input.p} -t {input.tf} --out_prefix {params.b} -n {threads} 2> {log.stderr} > {log.stdout}"

rule CAT_add_names_first_spades_assembly:
  input:
    i=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.contig2classification.txt",assemblytype='singles_hostfiltered'),
    t="references/CAT_prepare_20190108/2019-01-08_taxonomy"
  output:
     expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}_contig_taxonomy.tab",assemblytype='singles_hostfiltered')
  log:
    stdout=expand("logs/CAT_assembly_{assemblytype}_classification_taxonomy_{{hostcode}}.stdout",assemblytype='singles_hostfiltered'),
    stderr=expand("logs/CAT_assembly_{assemblytype}_classification_taxonomy_{{hostcode}}.stderr",assemblytype='singles_hostfiltered')
  threads: 1
  shell:
    "CAT add_names -i {input.i} -t {input.t} -o {output} > {log.stdout} 2> {log.stderr}"

rule CAT_filter_contignames_first_spades_assembly:
  input:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}_contig_taxonomy.tab",assemblytype='singles_hostfiltered')
  output:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}_contig_filterlist.txt",assemblytype='singles_hostfiltered')
  threads: 1
  log:
    expand("logs/CAT_assembly_{assemblytype}_contigfilterlist_{{hostcode}}.stderr",assemblytype='singles_hostfiltered')
  shell:
    "cat {input} | grep -v Bacteria | grep -v Fungi | grep -v Opisthokonta | grep -v Alveolata | grep Eukaryota | cut -f 1  > {output} 2> {log}"

rule create_filter_fasta_first_spades_assembly:
  input:
    n=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}_contig_filterlist.txt",assemblytype='singles_hostfiltered'),
    f=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs.fasta",assemblytype='singles_hostfiltered')
  output:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_filter_{{hostcode}}.fasta",assemblytype='singles_hostfiltered')
  threads: 1
  log:
    stdout=expand("logs/CAT_create_assembly_{assemblytype}_filter_{{hostcode}}.stdout",assemblytype='singles_hostfiltered'),
    stderr=expand("logs/CAT_create_assembly_{assemblytype}_filter_{{hostcode}}.stderr",assemblytype='singles_hostfiltered')
  shell:
    "samtools faidx {input.f} -o {output} -r {input.n} > {log.stdout} 2> {log.stderr}"

rule create_first_spades_assembly_filter_bt2_index:
  input:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_filter_{{hostcode}}.fasta",assemblytype='singles_hostfiltered')
  params:
    filter = lambda w: expand("data/assembly_{assemblytype}/{hostcode}/CAT_filter_bt2index_{hostcode}/{hostcode}_filter",assemblytype='singles_hostfiltered', hostcode=w.hostcode)
  output:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.{i}.bt2",i=range(1,4),assemblytype='singles_hostfiltered'),
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.rev.{i}.bt2",i=range(1,2),assemblytype='singles_hostfiltered'),
  threads: 12
  log:
    stdout=expand("logs/CAT_create_assembly_filter_{assemblytype}_filterbt2index_{{hostcode}}.stdout",assemblytype='singles_hostfiltered'),
    stderr=expand("logs/CAT_create_assembly_filter_{assemblytype}_filterbt2index_{{hostcode}}.stderr",assemblytype='singles_hostfiltered')
  shell:
    "bowtie2-build --threads {threads} {input} {params.filter} > {log.stdout} 2> {log.stderr}"

rule filter_for_assembly:
  input:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.{i}.bt2",i=range(1,4),assemblytype='singles_hostfiltered'),
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.rev.{i}.bt2",i=range(1,2),assemblytype='singles_hostfiltered'),
    s1=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=1,),
    s2=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=2)
  params:
    opts="--very-fast",
    i = lambda w : expand("data/assembly_{assemblytype}/{hostcode}/CAT_filter_bt2index_{hostcode}/{hostcode}_filter",assemblytype='singles_hostfiltered', hostcode = w.hostcode),
    outbase= lambda w : expand("data/sequencing_doublefiltered/{hostcode}/{hostcode}", hostcode=w.hostcode)
  output:
       expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}",PE=DIRECTIONS)
  threads: 36
  log:
    stderr="logs/bowtie2_filter_for_first_assembly_{hostcode}.stderr"
  shell:
    "bowtie2 {params.opts} --threads {threads} --un-conc-gz {params.outbase} -x {params.i} -1 {input.s1} -2 {input.s2}   > /dev/null 2> {log.stderr}"

rule rename_filtered_sequencing_files:
  input:
    a1=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}",PE=1),
    a2=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}",PE=2)
  output:
    b1=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=1),
    b2=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=2)
  shell:
    "mv {input.a1} {output.b1} && mv {input.a2} {output.b2}"

rule spades_second_assembly:
  input:
    reads=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=DIRECTIONS),
    s1=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=2)
  params:
    "--meta --only-assembler"
  output:
    basedir=expand("data/assembly_{assemblytype}/{{hostcode}}/",assemblytype='singles_doublefiltered'),
    contigs=expand("data/assembly_{assemblytype}/{{hostcode}}/contigs.fasta",assemblytype='singles_doublefiltered'),
    scaffolds=expand("data/assembly_{assemblytype}/{{hostcode}}/scaffolds.fasta",assemblytype='singles_doublefiltered'),
    graph=expand("data/assembly_{assemblytype}/{{hostcode}}/assembly_graph.fastg",assemblytype='singles_doublefiltered'),
    graph_scaffolds=expand("data/assembly_{assemblytype}/{{hostcode}}/assembly_graph_with_scaffolds.gfa",assemblytype='singles_doublefiltered'),
    datasetyaml=expand("data/assembly_{assemblytype}/{{hostcode}}/input_dataset.yaml",assemblytype='singles_doublefiltered'),
    paramfile=expand("data/assembly_{assemblytype}/{{hostcode}}/params.txt",assemblytype='singles_doublefiltered')
  threads: 100
  shadow: "shallow"
  resources:
    mem_mb=500
  log:
    stdout=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stdout",assemblytype='singles_doublefiltered'),
    stderr=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stderr",assemblytype='singles_doublefiltered')
  shell:
    "spades.py {params} -t {threads} -m {resources.mem_mb} -1 {input.s1} -2 {input.s2} -o {output.basedir} > {log.stdout} 2> {log.stderr}"

rule shorten scaffold_names:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds.fasta"
  output:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds_short_names.fasta"
  shell:
    ""#shell("awk -F "_" '/>NODE/{$0=">NODE_"$2}1' {input} > {output}")

rule bwa_index_assembly_scaffolds:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds_short_names.fasta"
  params:
    "data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds"
  output:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwa','pac','ann','sa','amb'])
  threads: 1
  log:
    stdout="logs/bwa_index_{assemblytype}_{hostcode}.stdout",
    stderr="logs/bwa_index_{assemblytype}_{hostcode}.stderr"
  shell:
    "bwa index -p {params} {input} > {log.stdout} 2> {log.stderr}"

READTYPES= ['.trimmed_paired.R','.trimmed']

import os.path
def get_binning_reads(wildcards):
    pathpe=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed_paired.R1.fastq.gz")
    pathse=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed.fastq.gz")
    if os.path.isfile(pathpe) == True :
      return {'reads' : expand("data/sequencing_binning_signals/{binningsignal}.trimmed_paired.R{PE}.fastq.gz", PE=[1,2],binningsignal=wildcards.binningsignal) }
    elif os.path.isfile(pathse) == True :
      return {'reads' : expand("data/sequencing_binning_signals/{binningsignal}.trimmed.fastq.gz", PE=[1,2],binningsignal=wildcards.binningsignal) }

rule backmap_bwa_mem:
  input:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwa','pac','ann','sa','amb']),
    unpack(get_binning_reads)
  params:
    "data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.bam"
  threads: 100
  log:
    stdout="logs/bwa_backmap_{assemblytype}_{hostcode}.stdout",
    stderr="logs/bwa_backmap_{assemblytype}_{hostcode}.stderr"
  shell:
    "bwa mem -t {threads} -p {params} {input.reads} | samtools view -b -o {output}"

rule allbackmapped:
  input:
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.bam",binningsignal=BINNINGSIGNALS,assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)

rule backmap_samtools_sort:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.bam"
