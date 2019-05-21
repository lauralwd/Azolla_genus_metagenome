# Ideally, format your hostcoes as 'species'_'accession'_'uniquenr/library/strain/treatment' for later automatic processing.
HOSTCODES= ['Azrub_IRRI_479','Azfil_lab_250', 'Azfil_lab_500', 'Azfil_lab_800', 'Azfil_minuscyano_170', 'Azfil_minuscyano_350', 'Azfil_wild_galgw_E_1', 'Azfil_wild_galgw_E_2', 'Azfil_wild_galgw_E_3', 'Azfil_wild_galgw_P_2', 'Azfil_wild_galgw_P_3', 'Azfil_wild_galgw_P_4', 'Azmex_IRRI_486', 'Azmic_IRRI_456', 'Aznil_IRRI_479', 'Azspnov_IRRI1_472', 'Azspnov_IRRI2_489']
DIRECTIONS=["1","2"]

# Assuming that host codes are formated as 'species'_'accession'_'uniquenr/library/strain/treatment' the following line of code
# autmatically creates an array of host accession combinations for cross-assembly of different strains,sequencing libraries, treatments.
# Likewise, the code below creates a list of species which have multiple sequencing libraries.

def get_species_cross_assembly():
  return_species=set([])
  for i in set([i.split('_',3)[0] for i in HOSTCODES]):
    species_sub=list(filter(lambda x:i in x, HOSTCODES))
    if len(species_sub) >= 2 :
      return_species.add(i)
  return(list(return_species))

def get_hosts_cross_assembly():
  return_hosts=set([])
  for i in set([i.split('_',3)[0] + '_' + i.split('_',3)[1] for i in HOSTCODES]):
    host_sub=list(filter(lambda x:i in x, HOSTCODES))
    if len(host_sub) >= 2 :
      return_hosts.add(i)
  return(list(return_hosts))

SPECIES=get_species_cross_assembly()
HOSTS=get_hosts_cross_assembly()

ASSEMBLYTYPES=['singles_doublefiltered','singles_hostfiltered'] # ,'hybrid_doublefiltered'] #,'species_doublefiltered']
# Ideally these should only contain letter, numbers and underscores. Exceptionally, these can contain points but they will be replaced by underscores for anvi'o
BINNINGSIGNALS=['dijkhuizen2018.E.1', 'dijkhuizen2018.E.2', 'dijkhuizen2018.E.3', 'dijkhuizen2018.P.2', 'dijkhuizen2018.P.3', 'dijkhuizen2018.P.4','ran2010.nostoc.SRR066216','ran2010.nostoc.SRR066217','ran2010.nostoc.SRR3923641','ran2010.nostoc.SRR3923645','ran2010.nostoc.SRR3923646']
ASSEMBLYFILES=['contigs','scaffolds']

## 'All'-rules
rule all:
  input:
    expand("data/bins_{assemblytype}_checkm/{hostcode}/{hostcode}.checkm_out",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES),
    expand("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE.db",assemblytype='singles_doublefiltered',hostcode=HOSTCODES),
    expand("data/assembly_{assemblytype}/{host}/contigs.fasta",assemblytype='hybrid_doublefiltered',host=HOSTS),
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype='hybrid_doublefiltered',hostcode=HOSTS,assemblyfile=ASSEMBLYFILES),
    expand("data/bins_{assemblytype}/{hostcode}.BAT.names.txt",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)
rule allhybridassemblies:
  input:
    contigs=expand("data/assembly_{assemblytype}/{host}/contigs.fasta",assemblytype='hybrid_doublefiltered',host=HOSTS)
rule alltaxtab:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype='singles_doublefiltered',hostcode=HOSTCODES,assemblyfile=ASSEMBLYFILES)
rule allcheckm:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/{hostcode}_depthmatrix.tab",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)
rule allbat:
  input:
    expand("data/bins_{assemblytype}/{hostcode}.BAT.bin2classification.txt",assemblytype='singles_hostfiltered',hostcode=HOSTCODES)

rule allsecondcat:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype='singles_doublefiltered',hostcode=HOSTCODES,assemblyfile=ASSEMBLYFILES)
rule allfirstcat:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_filterlist.txt",assemblyfile='scaffolds',hostcode=HOSTCODES,assemblytype='singles_hostfiltered')
rule allfastqc:
  input:
    expand("analyses/analyses_reads/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)
rule allfastqc_trimmed:
  input:
    expand("analyses/analyses_reads_trimmed/{hostcode}_{PE}", hostcode=HOSTCODES, PE=DIRECTIONS)
rule allfiltered:
  input:
    expand("data/sequencing_genomic_trimmed_filtered/{hostcode}.{PE}.fastq.gz",hostcode=HOSTCODES,PE=DIRECTIONS)
rule alldoublefiltered:
  input:
    expand("data/sequencing_doublefiltered/{hostcode}/{hostcode}.{PE}.fastq.gz",hostcode=HOSTCODES,PE=DIRECTIONS)
rule allfirstassemblies:
  input:
    expand("data/assembly_singles_hostfiltered/{hostcode}/contigs.fasta",hostcode=HOSTCODES)
rule allreadscorrected:
  input:
    expand("data/sequencing_genomic_trimmed_filtered_corrected/{hostcode}/corrected/{hostcode}.{PE}.fastq.00.0_0.cor.fastq.gz",hostcode=HOSTCODES, PE=DIRECTIONS)
rule allsecondassemblies:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/contigs.fasta",hostcode=HOSTCODES,assemblytype='singles_doublefiltered')
rule allbackmapped:
  input:
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.bam",       binningsignal=BINNINGSIGNALS,assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)

rule allsorted:
  input:
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.sorted.bam",binningsignal=BINNINGSIGNALS,assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES),
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)

rule allsourcemapped:
  input:
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.bam",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)
rule allsourcesorted:
  input:
    expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES)


## Read QC stuff
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
    db=temp("references/CAT_prepare_20190108/2019-01-08_CAT_database"),
    tf=temp("references/CAT_prepare_20190108/2019-01-08_taxonomy"),
    nr=temp("references/CAT_prepare_20190108/2019-01-08_CAT_database/2019-01-08.nr.gz")
  shell:
    "cd ./references && wget -qO - http://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20190108.tar.gz | tar -xz"

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
    "references/CAT_customised_20190108/CAT_database_customised/2019-03-27.nr.dmnd"
  log:
    stdout="logs/CAT_build_nr+host.stdout",
    stderr="logs/CAT_build_nr+host.stderr"
  threads: 100
  shell:
    "CAT prepare --existing -d {input.db} -t {input.tf} -n {threads} > {log.stdout} 2> {log.stderr}"

rule CAT_prepare_ORFS_host:
  input:
    "references/host_genome/host_genome.fasta"
  output:
    p="references/host_genome/host.predicted_proteins.faa",
    g="references/host_genome/host.predicted_proteins.gff",
  params:
    "-p meta -g 11 -q -f gff"
  threads: 1
  log:
    stdout="logs/CAT_host_prodigal.stdout",
    stderr="logs/CAT_host_prodigal.stderr"
  shell:
    "prodigal -i {input} -a {output.p} -o {output.g} {params} 2> {log.stderr} > {log.stdout}"

rule CAT_classify_host:
  input:
    prot="references/host_genome/host.predicted_proteins.faa",
    contigs="references/host_genome/host_genome.fasta",
    dmnd="references/CAT_customised_20190108/CAT_database_customised/2019-03-27.nr.dmnd",
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  params:
    prefix="references/host_genome/CAT/host",
    options="-r 5"
  output:
    i="references/host_genome/CAT/host.contig2classification.txt"
  threads: 100
  shadow: 'shallow'
  log:
    stdout="logs/CAT_classify_host.stdout",
    stderr="logs/CAT_classify_host.stderr"
  shell:
    "CAT contigs {params.options} -p {input.prot} -c {input.contigs} -d {input.db} -t {input.tf} --out_prefix {params.prefix} -n {threads} 2> {log.stderr} > {log.stdout}"
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
  params: "--only_official"
  threads: 1
  shell:
    "CAT add_names {params} -i {input.i} -t {input.tf} -o {output} > {log.stdout} 2> {log.stderr}"

rule CAT_filter_contignames:
  input:
     "references/host_genome/contig_taxonomy.tab"
  output:
     "references/host_genome/contig_filterlist.txt"
  threads: 1
  log:
    "logs/CAT_contigfilterlist.stderr"
  shell:
    "cat {input} | grep -v '#' | grep -v Bacteria | cut -f 1  > {output} 2> {log}"

rule CAT_get_bacterial_contignames:
  input:
     "references/host_genome/contig_taxonomy.tab"
  output:
     "references/host_genome/bacterial_contigs.txt"
  threads: 1
  log:
    "logs/CAT_contigfilterlist.stderr"
  shell:
    "cat {input} | grep -v '#' | grep Bacteria | cut -f 1  > {output} 2> {log}"

rule create_host_bacterialcontigs_fasta:
  input:
    n="references/host_genome/bacterial_contigs.txt",
    f="references/host_genome/host_genome.fasta"
  output:
    expand("references/host_genome/{host}_host-genome_bacterial_contigs.fasta",host='Azfil_lab')
  threads: 1
  log:
    stdout="logs/CAT_create_host_bacterial_contigs.stdout",
    stderr="logs/CAT_create_host_bacterial_contigs.stderr"
  shell:
    "samtools faidx {input.f} -o {output} -r {input.n} > {log.stdout} 2> {log.stderr}"

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
    p1=temp("data/sequencing_genomic_trimmed/{hostcode}_1.fastq.gz"),
    p2=temp("data/sequencing_genomic_trimmed/{hostcode}_2.fastq.gz")
  threads: 4
  resources: io=1
  log:
    stdout="logs/trimmomatic_genomicsequencing.{hostcode}.stdout",
    stderr="logs/trimmomatic_genomicsequencing.{hostcode}.stderr",
    log="logs/trimmomatic_genomicsequencing.{hostcode}.log",
    summary="logs/trimmomatic_genomicsequencing_summary{hostcode}.log"
  shell:
    "trimmomatic PE -threads {threads} -trimlog {log.log} -summary {log.summary} {input} {output.p1} /dev/null {output.p2} /dev/null {params} > {log.stdout} 2> {log.stderr}"

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
    temp(expand("data/sequencing_genomic_trimmed_filtered/{{hostcode}}.{PE}",PE=DIRECTIONS))
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
    options="--only-error-correction",
    basedir=lambda w: expand("data/sequencing_genomic_trimmed_filtered_corrected/{hostcode}/",hostcode=w.hostcode)
  output:
    reads=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=DIRECTIONS)
  threads: 100
  log:
    stdout="logs/SPAdes_correct_sequencing{hostcode}.stdout",
    stderr="logs/SPAdes_correct_sequencing{hostcode}.stderr"
  shell:
    "spades.py {params.options} -t {threads} -1 {input.s1} -2 {input.s2} -o {params.basedir} > {log.stdout} 2> {log.stderr}"

rule spades_first_assembly:
  input:
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
    dmnd="references/CAT_customised_20190108/CAT_database_customised/2019-03-27.nr.dmnd",
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf="references/CAT_customised_20190108/taxonomy_customised",
    p="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}_predicted_proteins.fasta"
  output:
    i="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}.contig2classification.txt",
#    g=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.predicted_proteins.gff",assemblytype='singles_hostfiltered'),
#    f=expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_{{hostcode}}.predicted_proteins.faa",assemblytype='singles_hostfiltered'),
    o="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}.ORF2LCA.txt",
    l="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}.log"
  shadow: 'shallow'
  params:
    b = lambda w: expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}", assemblytype=w.assemblytype, hostcode=w.hostcode, assemblyfile=w.assemblyfile),
    options = '-r 5'
  threads: 100
  resources:
    mem_mb=30000
  log:
    stdout="logs/CAT_assembly_{assemblytype}_{assemblyfile}_classification_{hostcode}.stdout",
    stderr="logs/CAT_assembly_{assemblytype}_{assemblyfile}_classification_{hostcode}.stderr"
  shell:
    "CAT contigs {params.options} -c {input.assembly} -d {input.db} -p {input.p} -t {input.tf} --out_prefix {params.b} -n {threads} 2> {log.stderr} > {log.stdout}"

rule CAT_add_names_assembly:
  input:
    i="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}.contig2classification.txt",
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
    "CAT add_names {params} -i {input.i} -t {input.tf} -o {output} > {log.stdout} 2> {log.stderr}"

rule CAT_filter_contignames_first_spades_assembly:
  input:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_taxonomy.tab",assemblyfile='scaffolds')
  output:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_filterlist.txt",assemblyfile='scaffolds')
  threads: 1
  log:
    expand("logs/CAT_assembly_{{assemblytype}}_{assemblyfile}filterlist_{{hostcode}}.stderr",assemblyfile='scaffolds')
  shell:
    "cat {input} | grep Eukaryota | cut -f 1 | cut -f 1,2 -d '_'  > {output} 2> {log}"

rule BAT_filter_contignames_bins:
  input:
    "data/bins_{assemblytype}/{hostcode}.BAT.names.txt"
  output:
    "data/assembly_{assemblytype}/{hostcode}/BAT_{hostcode}_filterlist.txt"
  params:
    binfolder = "data/bins_{assemblytype}/{hostcode}"
  log:
    "logs/BAT_filter_contignames_bins_{assemblytype}_{hostcode}"
  shell:
    """
    for f in  $(cat {input} | grep Eukaryota | cut -f 1)
    do  grep '>' {params.binfolder}/$f.fa
    done | tr -d '>' | cut -f 1 > {output} 2> {log}
    """

rule combine_filter_contignames:
  input:
    "data/assembly_{assemblytype}/{hostcode}/BAT_{hostcode}_filterlist.txt",
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_filterlist.txt",assemblyfile='scaffolds')
  output:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/combined_filterlist_{{hostcode}}_{assemblyfile}.tab",assemblyfile='scaffolds')
  shell:
    "cat {input} | cut -f 1 | sort | uniq > {output}"

rule create_filter_fasta_first_spades_assembly:
  input:
    n=expand("data/assembly_{{assemblytype}}/{{hostcode}}/combined_filterlist_{{hostcode}}_{assemblyfile}.tab",assemblyfile='scaffolds'),
    f=expand("data/assembly_{{assemblytype}}/{{hostcode}}/{assemblyfile}_short_names.fasta",assemblyfile='scaffolds')
  output:
    "data/assembly_{assemblytype}/{hostcode}/CAT_BAT_filter_{hostcode}.fasta"
  threads: 1
  log:
    stdout="logs/CAT_create_assembly_{assemblytype}_filter_{hostcode}.stdout",
    stderr="logs/CAT_create_assembly_{assemblytype}_filter_{hostcode}.stderr"
  shell:
    "samtools faidx {input.f} -o {output} -r {input.n} > {log.stdout} 2> {log.stderr}"

rule create_first_spades_assembly_filter_bt2_index:
  input:
    "data/assembly_{assemblytype}/{hostcode}/CAT_BAT_filter_{hostcode}.fasta"
  params:
    filter = lambda w: expand("data/assembly_{assemblytype}/{hostcode}/CAT_BAT_filter_bt2index_{hostcode}/{hostcode}_filter",assemblytype=w.assemblytype, hostcode=w.hostcode)
  output:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_BAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.{i}.bt2",i=range(1,4)),
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_BAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.rev.{i}.bt2",i=range(1,2))
  threads: 12
  log:
    stdout="logs/CAT_create_assembly_filter_{assemblytype}_filterbt2index_{hostcode}.stdout",
    stderr="logs/CAT_create_assembly_filter_{assemblytype}_filterbt2index_{hostcode}.stderr"
  shell:
    "bowtie2-build --threads {threads} {input} {params.filter} > {log.stdout} 2> {log.stderr}"

rule filter_for_assembly:
  input:
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_BAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.{i}.bt2",i=range(1,4),assemblytype='singles_hostfiltered'),
    expand("data/assembly_{assemblytype}/{{hostcode}}/CAT_BAT_filter_bt2index_{{hostcode}}/{{hostcode}}_filter.rev.{i}.bt2",i=range(1,2),assemblytype='singles_hostfiltered'),
    s1=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=1,),
    s2=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=2)
  params:
    opts="--very-sensitive",
    i = lambda w : expand("data/assembly_{assemblytype}/{hostcode}/CAT_BAT_filter_bt2index_{hostcode}/{hostcode}_filter",assemblytype='singles_hostfiltered', hostcode = w.hostcode),
    outbase= lambda w : expand("data/sequencing_doublefiltered/{hostcode}/{hostcode}", hostcode=w.hostcode)
  output:
    b1= expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=1),
    b2= expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=2),
    sin="data/sequencing_doublefiltered/{hostcode}/{hostcode}.singletons.fastq.gz"
  threads: 36
  log:
    stderr="logs/bowtie2_filter_for_assembly_doublefilter_{hostcode}.stderr",
    samstderr="logs/bowtie2_filter_for_assembly_doublefilter_{hostcode}_samtoolsfastq.stderr",
    samstdout="logs/bowtie2_filter_for_assembly_doublefilter_{hostcode}_samtoolsfastq.stdout"
  shell:
    """
    bowtie2 {params.opts} --threads {threads}	\
	-x {params.i}	\
	-1 {input.s1}	\
	-2 {input.s2}	\
		2> {log.stderr}	\
		| samtools fastq -f 13	\
			-@ {threads}	\
			-n -c 9		\
			-1 {output.b1}	\
			-2 {output.b2}	\
			-0 {output.sin}	\
			-		\
			2> {log.samstderr}	\
			> {log.samstdout}
    """

## double filtered assemblies
rule spades_second_assembly:
  input:
    reads=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=DIRECTIONS),
    s1=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_doublefiltered/{{hostcode}}/{{hostcode}}.{PE}.fastq.gz",PE=2)
  params:
    options="--meta --only-assembler",
    basedir=lambda w: expand("data/assembly_{assemblytype}/{hostcode}/",assemblytype='singles_doublefiltered',hostcode=w.hostcode)
  output:
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
    "spades.py {params.options} -t {threads} -m {resources.mem_mb} -1 {input.s1} -2 {input.s2} -o {params.basedir} > {log.stdout} 2> {log.stderr}"

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

ruleorder: shorten_scaffold_names_awk > shorten_scaffold_names_anvi

rule shorten_scaffold_names_awk:
  input:
    scaffolds=expand("data/assembly_{assemblytype}/{{hostcode}}/{{assemblyfile}}.fasta",assemblytype='singles_hostfiltered')
  output:
    scaffolds=expand("data/assembly_{assemblytype}/{{hostcode}}/{{assemblyfile}}_short_names.fasta",assemblytype='singles_hostfiltered')
  shell:
   """awk -F '_' '/>NODE/{{$0=">NODE_"$2}}1' {input} > {output}"""

rule shorten_scaffold_names_anvi:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}.fasta"
  output:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}_short_names.fasta"
  log:
    report="logs/anvi-script-reformat-fasta_{assemblytype}_{hostcode}_{assemblyfile}.report",
    stdout="logs/anvi-script-reformat-fasta_{assemblytype}_{hostcode}_{assemblyfile}.stdout",
    stderr="logs/anvi-script-reformat-fasta_{assemblytype}_{hostcode}_{assemblyfile}.stderr"
  conda:
    "envs/anvio.yaml"
  shell:
   "anvi-script-reformat-fasta -l 2500 --simplify-names -r {log.report} {input} -o {output} > {log.stdout} 2> {log.stderr} "

rule bwa_index_assembly_scaffolds:
  input:
    scaffolds=expand("data/assembly_{{assemblytype}}/{{hostcode}}/{assemblyfile}_short_names.fasta",assemblyfile='scaffolds')
  params:
    outbase=lambda w: expand("data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/{assemblyfile}",assemblyfile='scaffolds',assemblytype=w.assemblytype,hostcode=w.hostcode)
  output:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'])
  threads: 1
  log:
    stdout="logs/bwa_index_{assemblytype}_{hostcode}.stdout",
    stderr="logs/bwa_index_{assemblytype}_{hostcode}.stderr"
  shell:
    "bwa index -p {params.outbase} {input} > {log.stdout} 2> {log.stderr}"

import os.path
def get_binning_reads(wildcards):
    pathpe=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed_paired.R1.fastq.gz")
    pathse=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed.fastq.gz")
    if os.path.isfile(pathpe) ==  True :
      dict = {'reads' :  expand("data/sequencing_binning_signals/{binningsignal}.trimmed_paired.R{PE}.fastq.gz", PE=[1,2],binningsignal=wildcards.binningsignal) }
    elif os.path.isfile(pathse) == True :
      dict = {'reads' : expand("data/sequencing_binning_signals/{binningsignal}.trimmed.fastq.gz", binningsignal=wildcards.binningsignal) }
    return dict
    return dict

rule backmap_bwa_mem:
  input:
    unpack(get_binning_reads),
    index=expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'])
  params:
    index=lambda w: expand("data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds",assemblytype=w.assemblytype,hostcode=w.hostcode)
  output:
    temp("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.bam")
  threads: 100
  log:
    stdout="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    samstderr="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/bwa_backmap_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  shell:
    "bwa mem -t {threads} {params.index} {input.reads} 2> {log.stderr} | samtools view -F 4 -@ {threads} -b -o {output}  2> {log.samstderr} > {log.stdout}"

rule backmap_samtools_sort:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.bam"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.sorted.bam"
  threads: 100
  resources:
    mem_mb=2000
  log:
    stdout="logs/bwa_backmap_samtools_sort_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/bwa_backmap_samtools_sort_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  shell:
    "samtools sort -@ {threads} -m {resources.mem_mb}M -o {output} {input} > {log.stdout} 2> {log.stderr}"

rule backmap_bwa_mem_assemblysource:
  input:
    s1=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=2),
    index=expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'])
  params:
    index=lambda w: expand("data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds",assemblytype=w.assemblytype,hostcode=w.hostcode)
  output:
    temp("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.bam")
  threads: 100
  log:
    stdout="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}.stdout",
    samstderr="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}.stdout",
    stderr="logs/bwa_backmap_{assemblytype}_{hostcode}.stderr"
  shell:
    "bwa mem -t {threads} {params.index} {input.s1} {input.s2} 2> {log.stderr} | samtools view -F 4 -@ {threads} -b -o {output}  2> {log.samstderr} > {log.stdout}"

ruleorder: backmap_samtools_sort > backmap_samtools_sort_assemblysource

rule backmap_samtools_sort_assemblysource:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.bam"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam"
  threads: 100
  resources:
    mem_mb=2000
  log:
    stdout="logs/bwa_backmap_samtools_sort_{assemblytype}_{hostcode}_{hostcode}.stdout",
    stderr="logs/bwa_backmap_samtools_sort_{assemblytype}_{hostcode}_{hostcode}.stderr"
  shell:
    "samtools sort -@ {threads} -m {resources.mem_mb}M -o {output} {input} > {log.stdout} 2> {log.stderr}"

ruleorder: backmap_samtools_index > backmap_samtools_index_binningsignal

rule backmap_samtools_index:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam.bai"
  log:
    stdout="logs/bwa_backmap_samtools_index_{assemblytype}_{hostcode}_{hostcode}.stdout",
    stderr="logs/bwa_backmap_samtools_index_{assemblytype}_{hostcode}_{hostcode}.stderr"
  threads: 100
  shell:
    "samtools index -@ {threads} {input} > {log.stdout} 2> {log.stderr}"

rule backmap_samtools_index_binningsignal:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.sorted.bam"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.sorted.bam.bai"
  log:
    stdout="logs/bwa_backmap_samtools_index_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/bwa_backmap_samtools_index_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  threads: 100
  shell:
    "samtools index -@ {threads} {input} > {log.stdout} 2> {log.stderr}"

rule jgi_summarize_script:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam",
    expand("data/assembly_{{assemblytype}}_binningsignals/{{hostcode}}/{{hostcode}}_{binningsignal}.sorted.bam",binningsignal=BINNINGSIGNALS)
  output:
    "data/assembly_{assemblytype}/{hostcode}/{hostcode}_depthmatrix.tab"
  log:
    stdout="logs/jgi_summarize_script_{assemblytype}_{hostcode}.stdout",
    stderr="logs/jgi_summarize_script_{assemblytype}_{hostcode}.stdout"
  threads: 32
  shell:
    "jgi_summarize_bam_contig_depths --minContigLength 2500 --percentIdentity 80 --outputDepth {output} {input} > {log.stdout} 2> {log.stderr}"

rule metabat2:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds_short_names.fasta",
    depthmatrix="data/assembly_{assemblytype}/{hostcode}/{hostcode}_depthmatrix.tab"
  output:
    dir=directory("data/bins_{assemblytype}/{hostcode}")
  params:
    prefix=lambda w: expand("data/bins_{assemblytype}/{hostcode}/{hostcode}_bin",assemblytype=w.assemblytype,hostcode=w.hostcode)
  threads: 72
  log:
    stdout="logs/metabat2_{assemblytype}_{hostcode}.stdout",
    stderr="logs/metabat2_{assemblytype}_{hostcode}.stdout"
  shell:
    "metabat2 -t {threads} -i {input.scaffolds} -a {input.depthmatrix} -o {params.prefix} > {log.stdout} 2> {log.stderr}"

checkpoint dummy_metabat2:
  input:
    "data/bins_{assemblytype}/{hostcode}"
  output:
    "data/bins_{assemblytype}/{hostcode}/{hostcode}_bin.{bin_nr}.fa"

rule checkm:
  input:
    "data/bins_{assemblytype}/{hostcode}"
  output:
    table="data/bins_{assemblytype}_checkm/{hostcode}/{hostcode}.checkm_out"
  params:
    options="-x fa --pplacer_threads=12 --tab_table",
    dir=lambda w:expand("data/bins_{assemblytype}_checkm/{hostcode}",assemblytype=w.assemblytype,hostcode=w.hostcode)
  threads: 72
  log:
    stdout="logs/checkm_{assemblytype}_{hostcode}.stdout",
    stderr="logs/checkm_{assemblytype}_{hostcode}.stderr"
  conda:
    "envs/checkm.yaml"
  shell:
    "checkm lineage_wf -t {threads} {params.options} {input} {params.dir} -f {output.table} > {log.stdout} 2> {log.stderr}"

rule CAT_bins:
  input:
    bindir="data/bins_{assemblytype}/{hostcode}",
    dmnd="references/CAT_customised_20190108/CAT_database_customised/2019-03-27.nr.dmnd",
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
    "data/bins_{assemblytype}/{hostcode}.BAT.bin2classification.txt"
  shadow: 'shallow'
  params:
    options= " -s '.fa' ",
    prefix=lambda w : expand( "data/bins_{assemblytype}/{hostcode}.BAT" , assemblytype=w.assemblytype , hostcode=w.hostcode )
  threads: 72
  log:
    stdout="logs/BAT_{assemblytype}_{hostcode}.stdout",
    stderr="logs/BAT_{assemblytype}_{hostcode}.stderr"
  shell:
    "CAT bins -n {threads} -b {input.bindir} -d {input.db} -t {input.tf} {params.options} -o {params.prefix} > {log.stdout} 2> {log.stderr}"

rule BAT_add_names:
  input:
    i="data/bins_{assemblytype}/{hostcode}.BAT.bin2classification.txt",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
    "data/bins_{assemblytype}/{hostcode}.BAT.names.txt",
  params:
    "--only_official"
  log:
    stdout="logs/BAT_assembly_{assemblytype}_classification_taxonomy_{hostcode}.stdout",
    stderr="logs/BAT_assembly_{assemblytype}_classification_taxonomy_{hostcode}.stderr"
  threads: 1
  shell:
    "CAT add_names {params} -i {input.i} -t {input.tf} -o {output} > {log.stdout} 2> {log.stderr}"

rule anvi_gen_contigs_database:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds_short_names.fasta"
  output:
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db"
  log:
    stdout="logs/anvi-gen-contigs-database_{assemblytype}_{hostcode}.stdout",
    stderr="logs/anvi-gen-contigs-database_{assemblytype}_{hostcode}.stderr"
  params:
    "-n 'sample {hostcode} assembly {assemblytype}'"
  conda:
    "envs/anvio.yaml"
  shell:
    """
    anvi-gen-contigs-database -f {input} -o {output} {params} > {log.stdout} 2> {log.stderr}
    """

rule anvi_run_hmms:
  input:
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db"
  output:
    touch("data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs_db_run_hmms.done")
  threads: 5
  log:
    stdout="logs/anvi-run-hmms_{assemblytype}_{hostcode}.stdout",
    stderr="logs/anvi-run-hmms_{assemblytype}_{hostcode}.stderr"
  conda:
    "envs/anvio.yaml"
  shell:
    "anvi-run-hmms -c {input} -T {threads} > {log.stdout} 2> {log.stderr}"

ruleorder: anvi_profile_binningsignal > anvi_profile

rule anvi_profile:
  input:
    "data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs_db_run_hmms.done",
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db",
    bam="data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam",
    bai="data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{hostcode}.sorted.bam.bai"
  output:
    profile= "data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}_{hostcode}/PROFILE.db"
  params:
    length="--min-contig-length 2500",
    name=lambda w: expand(" -S 'assembly_{assemblytype}_sample_{hostcode}_binningsignal_{hostcode}' ", assemblytype=w.assemblytype , hostcode=w.hostcode.replace('.','_')),
    path= lambda w: expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}_{hostcode}", assemblytype=w.assemblytype , hostcode=w.hostcode)
  log:
    stdout="logs/anvi-profile_{assemblytype}_{hostcode}_{hostcode}.stdout",
    stderr="logs/anvi-profile_{assemblytype}_{hostcode}_{hostcode}.stderr"
  threads: 100
  conda:
    "envs/anvio.yaml"
  shell:
    """
    rmdir {params.path}
    anvi-profile -c {input.db} -i {input.bam} -o {params.path} -T {threads} {params.length} {params.name} > {log.stdout} 2> {log.stderr}
    """

rule anvi_profile_binningsignal:
  input:
    "data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs_db_run_hmms.done",
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db",
    bam="data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.sorted.bam",
    bai="data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}_{binningsignal}.sorted.bam.bai"
  output:
    profile="data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}_{binningsignal}/PROFILE.db"
  params:
    length="--min-contig-length 2500",
    name=lambda w: expand("-S assembly_{assemblytype}_sample_{hostcode}_binningsignal_{binningsignal}", assemblytype=w.assemblytype , hostcode=w.hostcode, binningsignal=w.binningsignal.replace('.','_')),
    path=lambda w: expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}_{binningsignal}", assemblytype=w.assemblytype , hostcode=w.hostcode, binningsignal=w.binningsignal)
  log:
    stdout="logs/anvi-profile_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/anvi-profile_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  threads: 100
  conda:
    "envs/anvio.yaml"
  shell:
    """
    rmdir {params.path}
    anvi-profile -c {input.db} -i {input.bam} -o {params.path} -T {threads} {params.length} {params.name} > {log.stdout} 2> {log.stderr}
    """

rule anvi_merge:
  input:
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db",
    source="data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}_{hostcode}/PROFILE.db",
    signal=expand("data/assembly_{{assemblytype}}_binningsignals_anvio/{{hostcode}}_{binningsignal}/PROFILE.db",binningsignal=BINNINGSIGNALS)
  output:
    profile="data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE.db"
  params:
    options="--enforce-hierarchical-clustering ",
    name= "-S 'assembly_{hostcode}_with_binningsignals'",
    path=lambda w:expand("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}",assemblytype=w.assemblytype , hostcode=w.hostcode)
  log:
    stdout="logs/anvi-merge-profile_{assemblytype}_{hostcode}.stdout",
    stderr="logs/anvi-merge-profile_{assemblytype}_{hostcode}.stderr"
  conda:
    "envs/anvio.yaml"
  shell:
    """
    rmdir {params.path}
    anvi-merge -c {input.db} -o {params.path} {params.options} {params.name} {input.source} {input.signal} > {log.stdout} 2> {log.stderr}
    """

def get_all_bins(wildcards):
    bins=checkpoints.metabat2.get(assemblytype='singles_doublefiltered',hostcode=wildcards.hostcode).output
    return bins

#rule prepare_anvi-import-metabat2:
#  input:
#    get_all_bins
#  output:
#    "data/bins_{assemblytype}/{hostcode}/{hostcode}_binlist.tab",
#  log:
#    "logs/prepare_anvi-import-metabat2.stderr
#  shell:
#    """
#    echo -e "bin\theader" > {output.binlist}
#    for   f in ( {input} )
#    do    cat $f | grep '>' | sed "s/^/bin_{bin_nr}\t/g" >> {output.binlist} 2> {log}
#    """

rule anvi_import_metabat2:
  input:
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db",
    profile="data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE.db",
    binlist="data/bins_{assemblytype}/{hostcode}/{hostcode}_binlist.tab"
  output:
    profile=touch("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE_db_imported-metabat2.done")
  params:
    "-C 'metabat2' ",
    "--contigs-mode"
  log:
    stdout="logs/anvi-merge-profile_{assemblytype}_{hostcode}.stdout",
    stderr="logs/anvi-merge-profile_{assemblytype}_{hostcode}.stderr"
  conda:
    "envs/anvio.yaml"
  shell:
    "anvi-import-collection {input.binlist} -c {input.db} {params} -p {input.profile} > {log.stdout} 2> {log.stderr}"


def get_input_hybrid_assemblies(wildcards):
    HOST=wildcards.host
    PE=wildcards.PE
    HOST_LIBRARIES=list(filter(lambda x:HOST in x, HOSTCODES))
    READS  = ['data/sequencing_doublefiltered/' + h + '/' + h +'.' + PE + '.fastq.gz' for h in HOST_LIBRARIES]
    return(READS)

#def get_input_hybrid_assemblies_commandline(wildcards):
#    HOST=wildcards.host
#    # these are the libraries that we'll assemble for this host, the HOSTS list should be difined such that this number is always bigger than 1
#    HOST_LIBRARIES=list(filter(lambda x:HOST in x, HOSTCODES))
#    # Count the number of libraries, SPAdes wants these explicitly numbered
#    COUNT = list(range(1,len(HOST_LIBRARIES)+1))
#    # construct paths for both pairs
#    LEFT  = ['data/sequencing_doublefiltered/' + h + '/' + h +'.1.fastq.gz' for h in HOST_LIBRARIES]
#    RIGHT = ['data/sequencing_doublefiltered/' + h + '/' + h +'.2.fastq.gz' for h in HOST_LIBRARIES]
#    # combine both lists into one list of lists
#    LEFTNUMBERED= list(zip(COUNT,LEFT ))
#    RIGHTNUMBERED=list(zip(COUNT,RIGHT))
#    # now use those lists of lists to construct commandlines
#    LEFTNAMED =['--pe' + str(h[0]) + '-1 ' + h[1] for h in LEFTNUMBERED ]
#    RIGHTNAMED=['--pe' + str(h[0]) + '-2 ' + h[1] for h in RIGHTNUMBERED]
#    READS=LEFTNAMED+RIGHTNAMED
#    return(READS)

rule concatenate_fastq_for_hybrid_assembly:
  input:
    get_input_hybrid_assemblies
  output:
    temp("data/sequencing_doublefiltered_concatenated/{host}.{PE}.fastq.gz")
  shell:
    "cat {input} > {output}"

rule SPADES_hybrid_assembly:
  input:
    s1=expand("data/sequencing_doublefiltered_concatenated/{{host}}.{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_doublefiltered_concatenated/{{host}}.{PE}.fastq.gz",PE=2),
    contigs="references/host_genome/{host}_host-genome_bacterial_contigs.fasta",
    pacbio="data/sequencing_genomic-longreads_trimmed_filtered/{host}_longreads-selfcorrected_trimmed_filtered.fasta.gz"
  output:
    contigs=expand("data/assembly_{assemblytype}/{{host}}/contigs.fasta",assemblytype='hybrid_doublefiltered'),
    scaffolds=expand("data/assembly_{assemblytype}/{{host}}/scaffolds.fasta",assemblytype='hybrid_doublefiltered'),
    graph=expand("data/assembly_{assemblytype}/{{host}}/assembly_graph.fastg",assemblytype='hybrid_doublefiltered'),
    graph_scaffolds=expand("data/assembly_{assemblytype}/{{host}}/assembly_graph_with_scaffolds.gfa",assemblytype='hybrid_doublefiltered'),
    datasetyaml=expand("data/assembly_{assemblytype}/{{host}}/input_dataset.yaml",assemblytype='hybrid_doublefiltered'),
    paramfile=expand("data/assembly_{assemblytype}/{{host}}/params.txt",assemblytype='hybrid_doublefiltered')
  params:
    options="--meta --only-assembler",
    basedir=lambda w: expand("data/assembly_{assemblytype}/{host}/",assemblytype='hybrid_doublefiltered',host=w.host)
  threads: 100
  shadow: "shallow"
  resources:
    mem_gb=500
  log:
    stdout=expand("logs/SPADES_assembly_{assemblytype}_{{host}}.stdout",assemblytype='hybrid_doublefiltered'),
    stderr=expand("logs/SPADES_assembly_{assemblytype}_{{host}}.stderr",assemblytype='hybrid_doublefiltered')
  shell:
    "spades.py {params.options} -t {threads} -m {resources.mem_gb} -1 {input.s1} -2 {input.s2} --pacbio {input.pacbio} --trusted-contigs {input.contigs} -o {params.basedir} > {log.stdout} 2> {log.stderr}"

rule unzip_long_reads_for_blasr:
  input:
    "data/sequencing_genomic-longreads_trimmed/{host}_longreads-selfcorrected_trimmed.fasta.gz"
  output:
    temp("data/sequencing_genomic-longreads_trimmed/{host}_longreads-selfcorrected_trimmed.fasta")
  shell:
    "unpigz {input} --keep "

rule BLASR_filter_long_reads:
  input:
    reads="data/sequencing_genomic-longreads_trimmed/{host}_longreads-selfcorrected_trimmed.fasta",
    genome="references/host_genome/host_filter.fasta"
  output:
    unaligned=temp("data/sequencing_genomic-longreads_trimmed_filtered/{host}_longreads-selfcorrected_trimmed_filtered.fasta"),
    bam="data/sequencing_genomic-longreads_trimmed_filtered/{host}_longreads-selfcorrected_trimmed_algined-to-host.bam"
  threads: 100
  params:
    "--hitPolicy allbest --bam"
  shadow: "shallow"
  conda:
    "envs/blasr.yaml"
  log:
    stdout="logs/BLASR_filter_long_reads_{host}.stdout",
    stderr="logs/BLASR_filter_long_reads_{host}.stderr"
  shell:
    "blasr {input.reads} {input.genome} {params} --nproc {threads} --out {output.bam} --unaligned {output.unaligned} > {log.stdout} 2> {log.stderr}"

rule compress_blasr_filtered_reads:
  input:
    "data/sequencing_genomic-longreads_trimmed_filtered/{host}_longreads-selfcorrected_trimmed_filtered.fasta"
  output:
    "data/sequencing_genomic-longreads_trimmed_filtered/{host}_longreads-selfcorrected_trimmed_filtered.fasta.gz"
  threads: 100
  shell:
    "pigz --best -j {threads} --keep data/sequencing_genomic-longreads_trimmed_filtered/{host}_longreads-selfcorrected_trimmed_filtered.fasta"
