HOSTCODES=["azca1_SRR6480231", "azca2_SRR6480201", "azfil_SRR6480158", "azfil_SRR6932851", "azmex_SRR6480159", "azmic_SRR6480161", "aznil_SRR6480196", "aznil_SRR6482158", "azrub_SRR6480160"]
#HOSTCODES= ['Azfil_lab_250', 'Azfil_lab_500', 'Azfil_lab_800', 'Azfil_minuscyano_170', 'Azfil_minuscyano_350', 'Azfil_minuscyano_HIC', 'Azfil_wild_galgw_E_1', 'Azfil_wild_galgw_E_2', 'Azfil_wild_galgw_E_3', 'Azfil_wild_galgw_P_2', 'Azfil_wild_galgw_P_3', 'Azfil_wild_galgw_P_4', 'Azmex_IRRI_486', 'Azmic_IRRI_456', 'Aznil_IRRI_479', 'Azrub_IRRI_479', 'Azspnov_IRRI_1_472', 'Azspnov_IRRI_2_489']
DIRECTIONS=["1","2"]

ASSEMBLYTYPES=['singles_doublefiltered','singles_hostfiltered'] # ,'hybrid_doublefiltered']
BINNINGSIGNALS=['dijkhuizen2018.E.1', 'dijkhuizen2018.E.2', 'dijkhuizen2018.E.3', 'dijkhuizen2018.P.2', 'dijkhuizen2018.P.3', 'dijkhuizen2018.P.4','ran2010.nostoc.SRR066216','ran2010.nostoc.SRR066217','ran2010.nostoc.SRR3923641','ran2010.nostoc.SRR3923645','ran2010.nostoc.SRR3923646']
ASSEMBLYFILES=['contigs','scaffolds']
## 'All'-rules

rule alltaxtab:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype='singles_doublefiltered',hostcode=HOSTCODES,assemblyfile=ASSEMBLYFILES)

rule allsecondcat:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype='singles_doublefiltered',hostcode=HOSTCODES,assemblyfile=ASSEMBLYFILES)
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
rule allbackmapped:
  input:
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.bam",       binningsignal=BINNINGSIGNALS,assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)
rule allsorted:
  input:
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.sorted.bam",binningsignal=BINNINGSIGNALS,assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)

rule half_fastq_file:
  input:
    "data/sequencing_genomic_trimmed_filtered_corrected/{hostcode}/corrected/{hostcode}.{PE}.fastq.00.0_0.cor.fastq.gz"
  output:
    "data/sequencing_genomic_trimmed_filtered_corrected_subset/subset_{hostcode}.{PE}.fastq.gz"
  threads: 5
  resources: io=1
  shell:
    """
    set +o pipefail
    zcat {input} | head -n $( echo "$(zcat {input} | wc -l) /2" / 2 | bc ) |  pigz -p 18 -c > {output}
    """

rule spades_subset_assembly:
  input:
    s1=expand("data/sequencing_genomic_trimmed_filtered_corrected_subset/subset_{{hostcode}}.{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed_filtered_corrected_subset/subset_{{hostcode}}.{PE}.fastq.gz",PE=2)
  params:
    basedir= lambda w: expand("data/assembly_{assemblytype}/{hostcode}/", hostcode=w.hostcode, assemblytype='singles_hostfiltered_subset'),
    options="--meta --only-assembler"
  output:
    contigs=expand("data/assembly_{assemblytype}/{{hostcode}}/{assemblyfile}.fasta",assemblytype='singles_hostfiltered_subset',assemblyfile='contigs'),
    scaffolds=expand("data/assembly_{assemblytype}/{{hostcode}}/{assemblyfile}.fasta",assemblytype='singles_hostfiltered_subset',assemblyfile='scaffolds')
  threads: 100
  shadow: "shallow"
  resources:
    mem_mb=500
  log:
    stdout=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stdout",assemblytype='singles_hostfiltered_subset'),
    stderr=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stderr",assemblyfile='contigs',assemblytype='singles_hostfiltered_subset')
  shell:
    "spades.py {params.options} -t {threads} --memory {resources.mem_mb} -1 {input.s1} -2 {input.s2} -o {params.basedir} > {log.stdout} 2> {log.stderr}"


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
    "analyses/analyses_reads_trimmed/{hostcode}_{PE}"
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

## reference rules
rule download_azolla_genome:
  output:
    "references/host_genome/host_genome.fasta"
  shell:
    "wget ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.genome_v1.2.fasta -O {output}"

rule download_azolla_proteins:
  output:
    "references/host_genome/host_proteins.fasta"
  shell:
    "wget ftp://ftp.fernbase.org/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.protein.highconfidence_v1.1.fasta -O {output}"

rule CAT_download:
  output:
    db="references/CAT_prepare_20190108/2019-01-08_CAT_database",
    tf="references/CAT_prepare_20190108/2019-01-08_taxonomy",
    nr="references/CAT_prepare_20190108/2019-01-08_CAT_database/2019-01-08.nr.gz"
  shell:
    "cd ./references && wget -qO - http://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20190108.tar.gz | tar -xz "

rule CAT_customise:
  input:
    db="references/CAT_prepare_20190108/2019-01-08_CAT_database/2019-01-08.nr.gz",
    tf="references/CAT_prepare_20190108/2019-01-08_taxonomy/",
    custom_proteins="references/host_genome/host_proteins.fasta"
  output:
    nr="references/CAT_customised_20190108/CAT_database_customised/2019-1-08.nr.gz",
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf=   "references/CAT_customised_20190108/taxonomy_customised",
    tf_id="references/CAT_customised_20190108/taxonomy_customised/2019-01-08.prot.accession2taxid.gz"
  shell:
    """
    cp	{input.db} {output.nr}
    cp -r {input.tf}/* {output.tf}
    pigz -c  {input.custom_proteins} >> {output.nr}
    grep '>' {input.custom_proteins} | tr -d '>' | awk -v OFS='\t' '{{print $0,  $0, 84609, 0}}' | pigz -c  >> {output.tf_id}
    """
rule CAT_build:
  input:
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
    "references/CAT_customised_20190108/CAT_database_customised/2019-03-13.nr.dmnd"
  log:
    stdout="logs/CAT_build_nr+host.stdout",
    stderr="logs/CAT_build_nr+host.stderr"
  threads: 100
  shell:
    "CAT prepare --existing -d {input.db} -t {input.tf} -n {threads} > {log.stdout} 2> {log.stderr}"

rule CAT_classify_host:
  input:
    c="references/host_genome/host_genome.fasta",
    dmnd="references/CAT_customised_20190108/CAT_database_customised/2019-03-13.nr.dmnd",
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
    i="references/host_genome/CAT/host.contig2classification.txt"
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
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
     "references/host_genome/contig_taxonomy.tab"
  log:
    stdout="logs/CAT_addnames_host.stdout",
    stderr="logs/CAT_addnames_host.stderr"
  threads: 1
  shell:
    "CAT add_names -i {input.i} -t {input.tf} -o {output} > {log.stdout} 2> {log.stderr}"

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
    f="references/host_genome/host_genome.fasta"
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
    expand("data/sequencing_genomic/{{hostcode}}_{PE}.fastq.gz", PE=DIRECTIONS)
  params:
    "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36"
  output:
    p1="data/sequencing_genomic_trimmed/{hostcode}_1.fastq.gz",
    p2="data/sequencing_genomic_trimmed/{hostcode}_2.fastq.gz",
  threads: 4
  resources: io=1
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
    basedir= lambda w: expand("data/assembly_{assemblytype}/{hostcode}/", hostcode=w.hostcode, assemblytype='singles_hostfiltered'),
    options="--meta --only-assembler"
  output:
    contigs=expand("data/assembly_{assemblytype}/{{hostcode}}/{assemblyfile}.fasta",assemblytype='singles_hostfiltered',assemblyfile='contigs'),
    scaffolds=expand("data/assembly_{assemblytype}/{{hostcode}}/{assemblyfile}.fasta",assemblytype='singles_hostfiltered',assemblyfile='scaffolds')
  threads: 100
  shadow: "shallow"
  resources:
    mem_mb=500
  log:
    stdout=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stdout",assemblytype='singles_hostfiltered'),
    stderr=expand("logs/SPADES_assembly_{assemblytype}_{{hostcode}}.stderr",assemblyfile='contigs',assemblytype='singles_hostfiltered')
  shell:
    "spades.py {params.options} -t {threads} -m {resources.mem_mb} -1 {input.s1} -2 {input.s2} -o {params.basedir} > {log.stdout} 2> {log.stderr}"

## process assembly for taxonomy including a taxonomy based second filter
rule CAT_prepare_ORFS:
  input:
    assembly="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}.fasta"
  output:
    p="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}_predicted_proteins.fasta",
    g="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}_predicted_proteins.gff"
  params:
    "-p meta -g 11 -q -f gff"
  threads: 1
  log:
    stdout="logs/CAT_assembly_{assemblytype}_{assemblyfile}_prodigal_{hostcode}.stdout",
    stderr="logs/CAT_assembly_{assemblytype}_{assemblyfile}_prodigal_{hostcode}.stderr"
  shell:
    "prodigal -i {input.assembly} -a {output.p} -o {output.g} {params} 2> {log.stderr} > {log.stdout}"

rule CAT_classify_contigs_assembly:
  input:
    assembly="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}.fasta",
    dmnd="references/CAT_customised_20190108/CAT_database_customised/2019-03-13.nr.dmnd",
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf="references/CAT_customised_20190108/taxonomy_customised",
    p="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}_predicted_proteins.fasta"
  output:
    i="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}2classification.txt",
#    g=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.predicted_proteins.gff",assemblytype='singles_hostfiltered'),
#    f=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.predicted_proteins.faa",assemblytype='singles_hostfiltered'),
    o="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}.ORF2LCA.txt",
    l="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}.log"
  shadow: "shallow"
  params: b = lambda w: expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}", assemblytype=w.assemblytype, hostcode=w.hostcode, assemblyfile=w.assemblyfile)
  threads: 100
  resources:
    mem_mb=30000
  log:
    stdout="logs/CAT_assembly_{assemblytype}_{assemblyfile}_classification_{hostcode}.stdout",
    stderr="logs/CAT_assembly_{assemblytype}_{assemblyfile}_classification_{hostcode}.stderr"
  shell:
    "CAT contigs -c {input.assembly} -d {input.db} -p {input.p} -t {input.tf} --out_prefix {params.b} -n {threads} 2> {log.stderr} > {log.stdout}"

rule CAT_add_names_assembly:
  input:
    i="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}.{assemblyfile}2classification.txt",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
     "data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab"
  params:
    "--only_official"
  log:
    stdout="logs/CAT_assembly_{assemblytype}_{assemblyfile}_classification_taxonomy_{hostcode}.stdout",
    stderr="logs/CAT_assembly_{assemblytype}_{assemblyfile}_classification_taxonomy_{hostcode}.stderr"
  threads: 1
  shell:
    "CAT add_names {params} -i {input.i} -t {input.t} -o {output} > {log.stdout} 2> {log.stderr}"

rule CAT_filter_contignames_first_spades_assembly:
  input:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_taxonomy.tab",assemblytype='singles_hostfiltered',assemblyfile='contigs')
  output:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_filterlist.txt",assemblytype='singles_hostfiltered',assemblyfile='contigs')
  threads: 1
  log:
    expand("logs/CAT_assembly_{assemblytype}_{assemblyfile}filterlist_{{hostcode}}.stderr",assemblytype='singles_hostfiltered',assemblyfile='contigs')
  shell:
    "cat {input} | grep -v Bacteria | grep -v Fungi | grep -v Opisthokonta | grep -v Alveolata | grep Eukaryota | cut -f 1  > {output} 2> {log}"

rule create_filter_fasta_first_spades_assembly:
  input:
    n=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_filterlist.txt",assemblytype='singles_hostfiltered',assemblyfile='contigs'),
    f=expand("data/assembly_{assemblytype}/{{hostcode}}/{assemblyfile}.fasta",assemblytype='singles_hostfiltered',assemblyfile='contigs')
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

## double filtered assemblies
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

## assembly analyses and diagnostigs
rule collect_assembly_stats:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES,assemblyfile=ASSEMBLYFILES)
  params:
    assemblytype= lambda w : expand("{assemblytype}",assemblytype=ASSEMBLYTYPES),
    hostcode= lambda w : expand("{hostcode}",hostcode=HOSTCODES),
    assemblyfile= lambda w : expand("{assemblyfile}",assemblyfile=ASSEMBLYFILES)
  output:
    "analyses/assembly_stats_and_taxonomy.tab.gz"
  threads: 12
  resources:
    mem_mb=1000
  shell:
    """
    scripts/make_assembly_stats_and_taxonomy.bash "{params.assemblytype}" "{params.hostcode}" "{params.assemblyfile}" {threads} {resources.mem_mb} {output}
    """

## assembly processing for binning an Anvi'o
rule shorten_scaffold_names:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds.fasta"
  output:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds_short_names.fasta"
  shell:
   """awk -F '_' '/>NODE/{{$0=">NODE_"$2}}1' {input} > {output}"""

rule bwa_index_assembly_scaffolds:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds_short_names.fasta"
  params:
    "data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds"
  output:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'])
  threads: 1
  log:
    stdout="logs/bwa_index_{assemblytype}_{hostcode}.stdout",
    stderr="logs/bwa_index_{assemblytype}_{hostcode}.stderr"
  shell:
    "bwa index -p {params} {input} > {log.stdout} 2> {log.stderr}"

import os.path
def get_binning_reads(wildcards):
    pathpe=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed_paired.R1.fastq.gz")
    pathse=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed.fastq.gz")
    if os.path.isfile(pathpe) == True :
      return {'reads' : expand("data/sequencing_binning_signals/{binningsignal}.trimmed_paired.R{PE}.fastq.gz", PE=[1,2],binningsignal=wildcards.binningsignal) }
    elif os.path.isfile(pathse) == True :
      return {'reads' : expand("data/sequencing_binning_signals/{binningsignal}.trimmed.fastq.gz", binningsignal=wildcards.binningsignal) }

rule backmap_bwa_mem:
  input:
    unpack(get_binning_reads),
    index=expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'])
  params:
    lambda w: expand("data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds",assemblytype=w.assemblytype,hostcode=w.hostcode)
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.bam"
  threads: 100
  log:
    stdout="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}.stdout",
    samstderr="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}.stdout",
    stderr="logs/bwa_backmap_{assemblytype}_{hostcode}.stderr"
  shell:
    "bwa mem -t {threads} {params} {input.reads} 2> {log.stderr} | samtools view -@ 12 -b -o {output}  2> {log.samstderr} > {log.stdout}"

rule backmap_samtools_sort:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.bam"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{binningsignal}.sorted.bam"
  threads: 6
  resources:
    mem_mb=5000
  shell:
    "samtools sort -@ {threads} -m {mem_mb}M -o {output} {input}"
