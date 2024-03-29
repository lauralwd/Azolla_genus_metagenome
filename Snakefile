# Ideally, format your hostcoes as 'species'_'accession'_'uniquenr/library/strain/treatment' for later automatic processing.
HOSTCODES= ['Azrub_IRRI_479','Azfil_lab_250', 'Azfil_lab_500', 'Azfil_lab_800', 'Azfil_minuscyano_170', 'Azfil_minuscyano_350', 'Azfil_wild_galgw_E_1', 'Azfil_wild_galgw_E_2', 'Azfil_wild_galgw_E_3', 'Azfil_wild_galgw_P_2', 'Azfil_wild_galgw_P_3', 'Azfil_wild_galgw_P_4', 'Azmex_IRRI_486', 'Azmic_IRRI_456', 'Aznil_IRRI_479', 'Azspnov_IRRI1_472', 'Azspnov_IRRI2_489']
DIRECTIONS=["1","2"]
# refined is a manual subset of the final libraries/assemblies chosen for final analysis
REFINED  = ['Azrub_IRRI_479', 'Azfil_lab', 'Azfil_minuscyano', 'Azfil_wild', 'Azmex_IRRI_486', 'Azmic_IRRI_456', 'Aznil_IRRI_479', 'Azspnov_IRRI1_472', 'Azspnov_IRRI2_489']



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

def get_single_hosts():
  return_hosts=set([])
  for i in HOSTCODES:
    species  = i.split('_',3)[0] + '_' + i.split('_',3)[1]
    host_sub = list(filter(lambda x:species in x, HOSTCODES))
    if len(host_sub) == 1 :
      return_hosts.add(i)
  return(list(return_hosts))

SPECIES=get_species_cross_assembly()
HOSTS=get_hosts_cross_assembly()
SINGLEHOSTS=get_single_hosts()

ASSEMBLYTYPES=['singles_doublefiltered','singles_hostfiltered'] # ,'hybrid_doublefiltered'] #,'species_doublefiltered']
# Ideally these should only contain letter, numbers and underscores. Exceptionally, these can contain points but they will be replaced by underscores for anvi'o
BINNINGSIGNALS=['dijkhuizen2018.E.1', 'dijkhuizen2018.E.2', 'dijkhuizen2018.E.3', 'dijkhuizen2018.P.2', 'dijkhuizen2018.P.3', 'dijkhuizen2018.P.4','ran2010.nostoc.SRR066216','ran2010.nostoc.SRR066217','ran2010.nostoc.SRR3923641','ran2010.nostoc.SRR3923645','ran2010.nostoc.SRR3923646']
ASSEMBLYFILES=['contigs','scaffolds']

## 'All'-rules
rule all_assemblies_and_annotations:
  input:
    "analyses/assembly_stats_and_taxonomy.tab",
    expand("data/bins_{assemblytype}/{hostcode}.BAT.names.txt",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES),
    expand("data/bins_{assemblytype}_checkm/{hostcode}/{hostcode}.checkm_out",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES),
    expand("data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs_db_run_ncbi_cogs.done",assemblytype='singles_doublefiltered',hostcode=HOSTCODES),
    expand("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/CAT_taxonomy_imported.done",assemblytype='singles_doublefiltered',hostcode=HOSTCODES),
    expand("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE_db_imported-metabat2.done",assemblytype='singles_doublefiltered',hostcode=HOSTCODES),
    "analyses/assembly-hybrid_stats_and_taxonomy.tab",
    expand("data/bins_{assemblytype}/{hostcode}.BAT.names.txt",assemblytype='hybrid_doublefiltered',hostcode=HOSTS),
    expand("data/bins_{assemblytype}_checkm/{hostcode}/{hostcode}.checkm_out",assemblytype='hybrid_doublefiltered',hostcode=HOSTS),
    expand("data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs_db_run_ncbi_cogs.done",assemblytype='hybrid_doublefiltered',hostcode=HOSTS),
    expand("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/CAT_taxonomy_imported.done",assemblytype='hybrid_doublefiltered',hostcode=HOSTS),
    expand("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE_db_imported-metabat2.done",assemblytype='hybrid_doublefiltered',hostcode=HOSTS)

rule all_exported_bins:
  input:
    expand("data/curated_bins/{{collection}}/{hostcode}",hostcode=HOSTS),
    expand("data/curated_bins/{{collection}}/{hostcode}",hostcode=SINGLEHOSTS)
  output:
    touch("data/curated_bins/{collection}.exported")

rule all_exported_bins_refined_checkm_BAT:
  input:
    expand("data/curated_bins/{collection}/{hostcode}/{hostcode}.checkm_out.txt", collection='refined', hostcode=REFINED),
    expand("data/curated_bins/{collection}/{hostcode}/{hostcode}.BAT.names.txt" , collection='refined', hostcode=REFINED)

rule trim:
  input:
    expand("data/sequencing_genomic_trimmed/{hostcode}_R{PE}.fastq.gz",hostcode=HOSTCODES,PE=DIRECTIONS)

## Read QC stuff and rules for a data flow figure
rule fastqc_raw_data:
  input:
    "data/sequencing_genomic/{hostcode}_R{PE}.fastq.gz"
  output:
    directory("analyses/analyses_reads/{hostcode}_{PE}")
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule fastqc_trimmed_data:
  input:
    "data/sequencing_genomic_trimmed/{hostcode}_{PE}.fastq.gz"
  output:
    directory("analyses/analyses_reads/{hostcode}_{PE}")
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule fastqc_filtered_data:
  input:
    "data/sequencing_genomic_trimmed_filtered_corrected/{hostcode}/corrected/{hostcode}.{PE}.fastq.00.0_0.cor.fastq.gz"
  output:
    directory("analyses/analyses_reads_trimmed_filtered/{hostcode}_{PE}")
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule fastqc_doublefiltered_data:
  input:
    "data/sequencing_doublefiltered/{hostcode}/{hostcode}.{PE}.fastq.gz"
  output:
    directory("analyses/analyses_reads_doublefiltered/{hostcode}_{PE}")
  shell:
    "mkdir {output} 2> /dev/null && fastqc -o {output} {input}"

rule collect_reads_stats:
  input:
    expand("analyses/analyses_reads/{hostcode}_{PE}"                  ,hostcode=HOSTCODES,PE=DIRECTIONS),
    expand("analyses/analyses_reads_trimmed/{hostcode}_{PE}"                  ,hostcode=HOSTCODES,PE=DIRECTIONS),
    expand("analyses/analyses_reads_trimmed_filtered/{hostcode}_{PE}" ,hostcode=HOSTCODES,PE=DIRECTIONS),
    expand("analyses/analyses_reads_doublefiltered/{hostcode}_{PE}"   ,hostcode=HOSTCODES,PE=DIRECTIONS)
#    expand("data/assembly_{assemblytype}/{hostcode}/{assemblyfile}.fasta",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES,assemblyfile='contigs')
  output:
    "analyses/reads_stats.tab"
  shell:
    "bash ./scripts/collect_reads_stats.bash > {output}"

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
#    db=temp("references/CAT_prepare_20190108/2019-01-08_CAT_database"),
    tf=temp("references/CAT_prepare_20190108/2019-01-08_taxonomy"),
    nr=temp("references/CAT_prepare_20190108/2019-01-08_CAT_database/2019-01-08.nr.gz")
  shell:
    "cd ./references && wget -qO - http://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20190108.tar.gz | tar -xz"

rule CAT_customise:
  input:
    db="references/CAT_prepare_20190108/2019-01-08_CAT_database/2019-01-08.nr.gz",
    tf="references/CAT_prepare_20190108/2019-01-08_taxonomy",
    custom_proteins="references/host_genome/host_proteins.fasta"
  output:
    nr="references/CAT_customised_20190108/CAT_database_customised/2019-1-08.nr.gz",
    db=directory("references/CAT_customised_20190108/CAT_database_customised"),
    tf=directory("references/CAT_customised_20190108/taxonomy_customised"),
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
    expand("data/sequencing_genomic/{{hostcode}}_R{PE}.fastq.gz", PE=DIRECTIONS)
  params:
    "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36"
  output:
    p1=temp("data/sequencing_genomic_trimmed/{hostcode}_R1.fastq.gz"),
    p2=temp("data/sequencing_genomic_trimmed/{hostcode}_R2.fastq.gz")
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
    s1=expand("data/sequencing_genomic_trimmed/{{hostcode}}_R{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed/{{hostcode}}_R{PE}.fastq.gz",PE=2)
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
    s1=ancient(expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=1)),
    s2=ancient(expand("data/sequencing_genomic_trimmed_filtered_corrected/{{hostcode}}/corrected/{{hostcode}}.{PE}.fastq.00.0_0.cor.fastq.gz",PE=2))
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
    "cat {input} | grep -v Eukaryota | cut -f 1 | sort -n  > {output} 2> {log}"

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
    contigs=protected(expand("data/assembly_{assemblytype}/{{hostcode}}/contigs.fasta",assemblytype='singles_doublefiltered')),
    scaffolds=protected(expand("data/assembly_{assemblytype}/{{hostcode}}/scaffolds.fasta",assemblytype='singles_doublefiltered')),
    graph=protected(expand("data/assembly_{assemblytype}/{{hostcode}}/assembly_graph.fastg",assemblytype='singles_doublefiltered')),
    graph_scaffolds=protected(expand("data/assembly_{assemblytype}/{{hostcode}}/assembly_graph_with_scaffolds.gfa",assemblytype='singles_doublefiltered')),
    datasetyaml=protected(expand("data/assembly_{assemblytype}/{{hostcode}}/input_dataset.yaml",assemblytype='singles_doublefiltered')),
    paramfile=protected(expand("data/assembly_{assemblytype}/{{hostcode}}/params.txt",assemblytype='singles_doublefiltered'))
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
rule collect_assembly_stats_singles:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype=ASSEMBLYTYPES,hostcode=HOSTCODES,assemblyfile=ASSEMBLYFILES)
  params:
    assemblytype= lambda w : expand("{assemblytype}",assemblytype=ASSEMBLYTYPES),
    hostcode= lambda w : expand("{hostcode}",hostcode=HOSTCODES),
    assemblyfile= lambda w : expand("{assemblyfile}",assemblyfile=ASSEMBLYFILES)
  output:
    "analyses/assembly_stats_and_taxonomy.tab"
  threads: 12
  resources:
    mem_mb=1000
  shell:
    """
    scripts/make_assembly_stats_and_taxonomy.bash "{params.assemblytype}" "{params.hostcode}" "{params.assemblyfile}" {threads} {resources.mem_mb} {output}
    """

rule collect_assembly_stats_hybrid:
  input:
    expand("data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",assemblytype='hybrid_doublefiltered',hostcode=HOSTS,assemblyfile=ASSEMBLYFILES)
  params:
    assemblytype= lambda w : expand("{assemblytype}",assemblytype='hybrid_doublefiltered'),
    hostcode= lambda w : expand("{hostcode}",hostcode=HOSTS),
    assemblyfile= lambda w : expand("{assemblyfile}",assemblyfile=ASSEMBLYFILES)
  output:
    "analyses/assembly-hybrid_stats_and_taxonomy.tab"
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

rule CAT_filter_unclassified_and_eukaryotic_scaffolds_for_anvio:
  input:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_taxonomy.tab",assemblyfile='scaffolds')
  output:
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_minus-unclassified_minus-eukaryotic_filterlist.txt",assemblyfile='scaffolds')
  threads: 1
  log:
    expand("logs/CAT_assembly_{{assemblytype}}_{assemblyfile}filterlist_{{hostcode}}.stderr",assemblyfile='scaffolds')
  shell:
    "cat {input} | grep -v Eukaryota | grep -v unclassified | grep -v '#' | cut -f 1 | sort -n  > {output} 2> {log}"

rule filter_unclassfied_filter_eukaryotic:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}.fasta",
    filterlist=expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_minus-unclassified_minus-eukaryotic_filterlist.txt",assemblyfile='scaffolds')
  output:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}_minus-unclassified_minus-eukaryotic.fasta"
  log:
    stdout="logs/filter-final-fasta-eukaryotic-unclassfied_{assemblytype}_{hostcode}_{assemblyfile}.stdout",
    stderr="logs/filter-final-fasta-eukaryotic-unclassified_{assemblytype}_{hostcode}_{assemblyfile}.stderr"
  shell:
    "samtools faidx {input.scaffolds} -o {output.scaffolds} -r {input.filterlist} > {log.stdout} 2> {log.stderr}"

rule shorten_scaffold_names_anvi:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/{assemblyfile}_minus-unclassified_minus-eukaryotic.fasta"
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
    index={'index': expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb']) }
    if wildcards.assemblytype != 'hybrid_doublefiltered' :
        pathpe=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed_paired.R1.fastq.gz")
        pathse=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed.fastq.gz")
        if os.path.isfile(pathpe) ==  True :
            dict = {'reads' :  expand("data/sequencing_binning_signals/{binningsignal}.trimmed_paired.R{PE}.fastq.gz", PE=[1,2],binningsignal=wildcards.binningsignal) }
        elif os.path.isfile(pathse) == True :
            dict = {'reads' : expand("data/sequencing_binning_signals/{binningsignal}.trimmed.fastq.gz", binningsignal=wildcards.binningsignal) }
            return dict
        dict.update(index)
        return dict
    elif wildcards.assemblytype == 'hybrid_doublefiltered' :
        if len(list(filter(lambda x:wildcards.binningsignal in x, BINNINGSIGNALS))) > 0 :
            pathpe=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed_paired.R1.fastq.gz")
            pathse=("data/sequencing_binning_signals/" + wildcards.binningsignal + ".trimmed.fastq.gz")
            index={'index': expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb']) }
            if os.path.isfile(pathpe) ==  True :
                dict = {'reads' :  expand("data/sequencing_binning_signals/{binningsignal}.trimmed_paired.R{PE}.fastq.gz", PE=[1,2],binningsignal=wildcards.binningsignal) }
            elif os.path.isfile(pathse) == True :
                dict = {'reads' : expand("data/sequencing_binning_signals/{binningsignal}.trimmed.fastq.gz", binningsignal=wildcards.binningsignal) }
                return dict
            dict.update(index)
        elif len(list(filter(lambda x:wildcards.binningsignal in x, BINNINGSIGNALS))) == 0 :
            dict = { 'reads' : expand("data/sequencing_genomic_trimmed/{hostcode}_R{PE}.fastq.gz",PE=DIRECTIONS,hostcode=wildcards.binningsignal) }
        dict.update(index)
        return dict
    dict.update(index)
    return dict

rule backmap_bwa_mem:
  input:
    unpack(get_binning_reads),
    expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'])
  params:
    index=lambda w: expand("data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds",assemblytype=w.assemblytype,hostcode=w.hostcode)
  output:
    temp("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.bam")
  threads: 100
  log:
    stdout="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    samstderr="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/bwa_backmap_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  shell:
    "bwa mem -t {threads} {params.index} {input.reads} 2> {log.stderr} | samtools view -F 4 -@ {threads} -b -o {output}  2> {log.samstderr} > {log.stdout}"

# def get_source_binning_reads(wildcards):
#     if wildcards.assemblytype != 'hybrid_doublefiltered':
#         s1={'s1' : expand("data/sequencing_genomic_trimmed/{hostcode}_R{PE}.fastq.gz",PE=1,hostcode=wildcards.hostcode) }
#         s2={'s2' : expand("data/sequencing_genomic_trimmed/{hostcode}_R{PE}.fastq.gz",PE=2,hostcode=wildcards.hostcode) }
#         index={'index' : expand("data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'],hostcode=wildcards.hostcode,assemblytype=wildcards.assemblytype)}
#         input={}
#         input.update(s1)
#         input.update(s2)
#         input.update(index)
#         (input)
#     elif wildcards.assemblytype == 'hybrid_doublefiltered':
#         s1={'s1' : expand("data/sequencing_genomic_trimmed/{hostcode}_R{PE}.fastq.gz",PE=1,hostcode=wildcards.hostcode) }
#         s2={'s2' : expand("data/sequencing_genomic_trimmed/{hostcode}_R{PE}.fastq.gz",PE=2,hostcode=wildcards.hostcode) }
#         index={'index' : expand("data/assembly_{assemblytype}/{host}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'],host=wildcards.hostcode,assemblytype=wildcards.assemblytype)}
#         input={}
#         input.update(s1)
#         input.update(s2)
#         input.update(index)
#         return(input)
#     return(input)

ruleorder: BLASR_backmap_long_reads > backmap_bwa_mem_assemblysource

rule BLASR_backmap_long_reads:
  input:
    reads="data/sequencing_genomic-longreads_trimmed/Azfil_lab_longreads-selfcorrected_trimmed.fasta",
    scaffolds=expand("data/assembly_{{assemblytype}}/{{hostcode}}/{assemblyfile}_short_names.fasta",assemblyfile='scaffolds')
  output:
    temp("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+pacbio_reads.bam")
  threads: 100
  params:
    "--hitPolicy allbest --bam"
  shadow: "shallow"
  conda:
    "envs/blasr.yaml"
  log:
    stdout="logs/BLASR_filter_long_reads_{assemblytype}_{hostcode}.stdout",
    stderr="logs/BLASR_filter_long_reads_{assemblytype}_{hostcode}.stderr"
  shell:
    "blasr {input.reads} {input.scaffolds} {params} --nproc {threads} --out {output}  > {log.stdout} 2> {log.stderr}"

rule backmap_bwa_mem_assemblysource:
  input:
    #unpack(get_source_binning_reads)
    s1=expand("data/sequencing_genomic_trimmed/{{hostcode}}_R{PE}.fastq.gz",PE=1),
    s2=expand("data/sequencing_genomic_trimmed/{{hostcode}}_R{PE}.fastq.gz",PE=2),
    index=expand("data/assembly_{{assemblytype}}/{{hostcode}}/scaffolds_bwa_index/scaffolds.{ext}",ext=['bwt','pac','ann','sa','amb'])
  params:
    index=lambda w: expand("data/assembly_{assemblytype}/{hostcode}/scaffolds_bwa_index/scaffolds",assemblytype=w.assemblytype,hostcode=w.hostcode)
  output:
    temp("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{hostcode}.bam")
  threads: 100
  log:
    stdout="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}.stdout",
    samstderr="logs/bwa_backmap_samtools_{assemblytype}_{hostcode}.stdout",
    stderr="logs/bwa_backmap_{assemblytype}_{hostcode}.stderr"
  shell:
    "bwa mem -t {threads} {params.index} {input.s1} {input.s2} 2> {log.stderr} | samtools view -F 4 -@ {threads} -b -o {output}  2> {log.samstderr} > {log.stdout}"

rule backmap_samtools_sort:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.bam"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.sorted.bam"
  threads: 100
  resources:
    mem_mb=2000
  log:
    stdout="logs/bwa_backmap_samtools_sort_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/bwa_backmap_samtools_sort_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  shell:
    "samtools sort -@ {threads} -m {resources.mem_mb}M -o {output} {input} > {log.stdout} 2> {log.stderr}"

rule backmap_samtools_index_binningsignal:
  input:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.sorted.bam"
  output:
    "data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.sorted.bam.bai"
  log:
    stdout="logs/bwa_backmap_samtools_index_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/bwa_backmap_samtools_index_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  threads: 100
  shell:
    "samtools index -@ {threads} {input} > {log.stdout} 2> {log.stderr}"

def get_bams_for_binning(wildcards):
    if wildcards.assemblytype != 'hybrid_doublefiltered':
        input=expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{hostcode}.sorted.bam",assemblytype=wildcards.assemblytype,hostcode=wildcards.hostcode) +  expand("data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.sorted.bam",binningsignal=BINNINGSIGNALS,assemblytype=wildcards.assemblytype,hostcode=wildcards.hostcode)
        return(input)
    elif wildcards.assemblytype == 'hybrid_doublefiltered':
      HOST_LIBRARIES=list(filter(lambda x:wildcards.hostcode in x, HOSTCODES)) + BINNINGSIGNALS
      input=expand("data/assembly_{assemblytype}_binningsignals/{host}/{host}+{hostcode}.sorted.bam",assemblytype=wildcards.assemblytype,host=wildcards.hostcode,hostcode=HOST_LIBRARIES)
      return(input)

rule jgi_summarize_script:
  input:
    get_bams_for_binning
  output:
    "data/assembly_{assemblytype}/{hostcode}/{hostcode}_depthmatrix.tab"
  log:
    stdout="logs/jgi_summarize_script_{assemblytype}_{hostcode}.stdout",
    stderr="logs/jgi_summarize_script_{assemblytype}_{hostcode}.stdout"
  threads: 32
  shell:
    "jgi_summarize_bam_contig_depths --minContigLength 2500 --percentIdentity 80 --outputDepth {output} {input} > {log.stdout} 2> {log.stderr}"

checkpoint metabat2:
  input:
    scaffolds="data/assembly_{assemblytype}/{hostcode}/scaffolds_short_names.fasta",
    depthmatrix="data/assembly_{assemblytype}/{hostcode}/{hostcode}_depthmatrix.tab"
  output:
    bins=directory("data/bins_{assemblytype}/{hostcode,[A-Za-z0-9_]+}")
  params:
    prefix=lambda w: expand("data/bins_{assemblytype}/{hostcode}/{hostcode}_bin",assemblytype=w.assemblytype,hostcode=w.hostcode)
  threads: 72
  log:
    stdout="logs/metabat2_{assemblytype}_{hostcode}.stdout",
    stderr="logs/metabat2_{assemblytype}_{hostcode}.stdout"
  shell:
    "metabat2 -t {threads} -i {input.scaffolds} -a {input.depthmatrix} -o {params.prefix} > {log.stdout} 2> {log.stderr}"

rule checkm_set_data_folder:
  input:
    "references/checkm_data"
  output:
    touch("references/checkm_data_setroot.done")
  conda:
    "envs/checkm.yaml"
  shell:
    "checkm data setRoot {input}"

rule checkm:
  input:
    bins="data/bins_{assemblytype}/{hostcode}",
    set_root="references/checkm_data_setroot.done"
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
    """
    if [ -d {params.dir} ]
    then rm -rf {params.dir}
    fi
    checkm lineage_wf -t {threads} {params.options} {input.bins} {params.dir} -f {output.table} > {log.stdout} 2> {log.stderr}
    """
rule prodigal_get_ORFs_for_CAT_bins:
  input:
    bindir="data/bins_{assemblytype}/{hostcode}"
  output:
    a="data/bins_{assemblytype}/{hostcode}.BAT.concatenated.predicted_proteins.faa",
    o="data/bins_{assemblytype}/{hostcode}.BAT.concatenated.predicted_proteins.gff"
  params:
    "-p meta -g 11 -q -f gff"
  shell:
    "prodigal -i <(cat {input}/*.fa ) -a {output.a} -o {output.o} {params}"

rule CAT_bins:
  input:
#    a="data/bins_{assemblytype}/{hostcode}.BAT.concatenated.predicted_proteins.faa",
#    o="data/bins_{assemblytype}/{hostcode}.BAT.concatenated.predicted_proteins.gff",
    bindir="data/bins_{assemblytype}/{hostcode}",
    dmnd="references/CAT_customised_20190108/CAT_database_customised/2019-03-27.nr.dmnd",
    db="references/CAT_customised_20190108/CAT_database_customised",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
    "data/bins_{assemblytype}/{hostcode}.BAT.bin2classification.txt"
#  shadow: 'shallow'
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

rule anvi_setup_ncbi_cogs:
  output:
    dir=directory("references/anvi_ncbi_cogs")
  threads: 100
  log:
    stdout="logs/anvi-setup-cogs.stdout",
    stderr="logs/anvi-setup-cogs.stderr"
  conda:
    "envs/anvio.yaml"
  shell:
    "anvi-setup-ncbi-cogs -T {threads} --just-do-it --cog-data-dir {output.dir} > {log.stdout} 2> {log.stderr}"

rule anvi_run_ncbi_cogs:
  input:
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db",
    dir="references/anvi_ncbi_cogs"
  output:
    touch("data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs_db_run_ncbi_cogs.done")
  threads: 20
  params: "--sensitive"
  log:
    stdout="logs/anvi-run-cogs_{assemblytype}_{hostcode}.stdout",
    stderr="logs/anvi-run-cogs_{assemblytype}_{hostcode}.stderr"
  conda:
    "envs/anvio.yaml"
  shell:
    "anvi-run-ncbi-cogs {params} -c {input.db} -T {threads} --cog-data-dir {input.dir} > {log.stdout} 2> {log.stderr}"

rule prepare_anvi_import_cat_taxonomy:
  input:
    report="logs/anvi-script-reformat-fasta_{assemblytype}_{hostcode}_{assemblyfile}.report",
    taxonomy="data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy.tab",
    profile=ancient("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE.db")
  output:
    "data/assembly_{assemblytype}/{hostcode}/CAT_{hostcode}_{assemblyfile}_taxonomy_shortnames.tab"
  threads: 2
  conda:
    "envs/anvio.yaml"
  log:
    "logs/prepare-anvi-import-cat-taxonomy_{assemblytype}_{hostcode}_{assemblyfile}.stderr"
  shell:
    """
    echo "item_name\tcategorical_kingdom\tcategorical_phylum\tcategorical_class\tcategorical_order\tcategorical_family\tcategorical_genus\tcategorical_species" \
		> {output} 2> {log}
    join -1 2 -2 1						\
		<( sort -k2d {input.report} )			\
		<(tail -n +2 {input.taxonomy}			\
			| cut -f 1,7-				\
			| sed 's/: [01]\.[0-9][0-9]//g'		\
			| tr ' ' '_'				\
			| sort -k1d				)\
	| tr ' ' '\t'					\
	| cut -f 2-					\
	| sort -n					\
		 >> {output}.tmp 2>> {log}
    join -1 1 -2 1 --check-order -a 1						\
		<(anvi-get-split-coverages --list-splits -p {input.profile}	\
			2>> {log}						\
			| sed 's/_split/\t_split/g' | sort -k1d)		\
		<( sort -k1d {output}.tmp )					\
		| tr ' ' '\t'							\
		| sort -k1n,2n							\
		| sed 's/\t_split/_split/g'					\
			>> {output}  2>> {log}
    rm {output}.tmp 2>> {log}
    """

rule anvi_import_cat_taxonomy:
  input:
    taxonomy=expand("data/assembly_{{assemblytype}}/{{hostcode}}/CAT_{{hostcode}}_{assemblyfile}_taxonomy_shortnames.tab",assemblyfile='scaffolds'),
    profile=ancient("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE.db")
  output:
    touch("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/CAT_taxonomy_imported.done")
  threads: 1
  params:
    "--target-data-table items --just-do-it"
  conda:
    "envs/anvio.yaml"
  log:
    stdout="logs/anvi-import_cat_taxonomy_{assemblytype}_{hostcode}.stdout",
    stderr="logs/anvi-import_cat_taxonomy_{assemblytype}_{hostcode}.stderr"
  shell:
    "anvi-import-misc-data {params} -p {input.profile} {input.taxonomy}  > {log.stdout} 2> {log.stderr}"

rule anvi_profile_binningsignal:
  input:
    "data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs_db_run_hmms.done",
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db",
    bam="data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.sorted.bam",
    bai="data/assembly_{assemblytype}_binningsignals/{hostcode}/{hostcode}+{binningsignal}.sorted.bam.bai"
  output:
    profile="data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}+{binningsignal}/PROFILE.db"
  params:
    length="--min-contig-length 2500",
    name=lambda w: expand("-S assembly_{assemblytype}_sample_{hostcode}_binningsignal_{binningsignal}", assemblytype=w.assemblytype , hostcode=w.hostcode, binningsignal=w.binningsignal.replace('.','_')),
    path=lambda w: expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}+{binningsignal}", assemblytype=w.assemblytype , hostcode=w.hostcode, binningsignal=w.binningsignal)
  log:
    stdout="logs/anvi-profile_{assemblytype}_{hostcode}_{binningsignal}.stdout",
    stderr="logs/anvi-profile_{assemblytype}_{hostcode}_{binningsignal}.stderr",
    preperr="logs/anvi-profile_prep_{assemblytype}_{hostcode}_{binningsignal}.stderr"
  threads: 100
  conda:
    "envs/anvio.yaml"
  shell:
    """
    if [ -d {params.path} ]
    then rm -rf {params.path} > {log.preperr} 2>&1
    fi
    anvi-profile -c {input.db} -i {input.bam} -o {params.path} -T {threads} {params.length} {params.name} > {log.stdout} 2> {log.stderr}
    """

def get_anvi_merge_profiles(wildcards):
    if wildcards.assemblytype != 'hybrid_doublefiltered':
        source= expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}+{hostcode}/PROFILE.db",hostcode=wildcards.hostcode,assemblytype=wildcards.assemblytype)
        signal= expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}+{binningsignal}/PROFILE.db",binningsignal=BINNINGSIGNALS,hostcode=wildcards.hostcode,assemblytype=wildcards.assemblytype)
        signal.extend(source)
    elif wildcards.assemblytype == 'hybrid_doublefiltered':
        HOST=wildcards.hostcode
        HOST_LIBRARIES=list(filter(lambda x:wildcards.hostcode in x, HOSTCODES))
        signal=expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}+{library}/PROFILE.db",    library=HOST_LIBRARIES,    assemblytype=wildcards.assemblytype,hostcode=wildcards.hostcode)
        SPECIES=HOST.split('_',1)[0]
        longreadfiles=listdir("data/sequencing_genomic-longreads_trimmed/")
        longreadfiles_species=list(filter(lambda x:SPECIES in x, longreadfiles))
        if wildcards.hostcode != 'Azfil_wild':
            binningsignal = expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}+{binningsignal}/PROFILE.db",binningsignal=BINNINGSIGNALS,assemblytype=wildcards.assemblytype,hostcode=wildcards.hostcode)
            signal.extend(binningsignal)
        if len(longreadfiles_species) > 0 :
            pacbio = expand("data/assembly_{assemblytype}_binningsignals_anvio/{hostcode}+pacbio_reads/PROFILE.db",assemblytype=wildcards.assemblytype,hostcode=wildcards.hostcode)
            signal.extend(pacbio)
    dict={'signal' : signal }
    return(dict)

rule anvi_merge:
  input:
    unpack(get_anvi_merge_profiles),
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db"
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
    if [ -d {params.path} ]
    then rm -rf {params.path}
    fi
    anvi-merge -c {input.db} -o {params.path} {params.options} {params.name} {input.signal} > {log.stdout} 2> {log.stderr}
    """

def get_all_bins(wildcards):
    bins = checkpoints.metabat2.get(**wildcards).output[0]
    input= expand("data/bins_{assemblytype}/{hostcode}/{hostcode}_bin.{bin_nr}.fa",
                  assemblytype=wildcards.assemblytype,
                  hostcode=wildcards.hostcode,
                  bin_nr=glob_wildcards(os.path.join(bins,"{hostcode}_bin.{bin_nr}.fa")).bin_nr
                  )
    return input

rule prepare_anvi_import_metabat2:
  input:
    get_all_bins
  output:
    temp("data/bins_{assemblytype}/{hostcode}/{hostcode}_binlist.tab")
  log:
    "logs/prepare_anvi-import-metabat2-{assemblytype}-{hostcode}.stderr"
  threads: 3
  shell:
    """
    bins=( $(echo "{input}" | tr ' ' '\n' ) )
    for   f in ${{bins[@]}}
    do    nr=$(echo $f | rev | cut -f 2 -d '.' )
          cat $f | grep '>' | tr -d '>' | sed "s/$/\tbin$nr/g"  >> {output} 2> {log}
    done
    """

rule anvi_import_metabat2:
  input:
    db="data/assembly_{assemblytype}_anvio/{hostcode}/{hostcode}_contigs.db",
    profile=ancient("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE.db"),
    binlist="data/bins_{assemblytype}/{hostcode}/{hostcode}_binlist.tab"
  output:
    profile=touch("data/assembly_{assemblytype}_binningsignals_anvio/MERGED_{hostcode}/PROFILE_db_imported-metabat2.done")
  params:
    "-C 'metabat2' ",
    "--contigs-mode"
  log:
    stdout="logs/anvi-import-metabat2_{assemblytype}_{hostcode}.stdout",
    stderr="logs/anvi-import-metabat2_{assemblytype}_{hostcode}.stderr"
  conda:
    "envs/anvio.yaml"
  shell:
    "anvi-import-collection {input.binlist} -c {input.db} {params} -p {input.profile} > {log.stdout} 2> {log.stderr}"

def get_concatenate_hybrid_assemblies(wildcards):
    HOST=wildcards.host
    PE=wildcards.PE
    HOST_LIBRARIES=list(filter(lambda x:HOST in x, HOSTCODES))
    READS  = ['data/sequencing_doublefiltered/' + h + '/' + h +'.' + PE + '.fastq.gz' for h in HOST_LIBRARIES]
    return(READS)

rule concatenate_fastq_for_hybrid_assembly:
  input:
    get_concatenate_hybrid_assemblies
  output:
    temp("data/sequencing_doublefiltered_concatenated/{host}.{PE}.fastq.gz")
  shell:
    "cat {input} > {output}"

from os import listdir
def get_input_hybrid_assemblies(wildcards):
    HOST=wildcards.host
    SPECIES=HOST.split('_',1)[0]
    s1  = {'s1' : 'data/sequencing_doublefiltered_concatenated/' + HOST +'.' + '1' + '.fastq.gz' }
    s2  = {'s2' : 'data/sequencing_doublefiltered_concatenated/' + HOST +'.' + '2' + '.fastq.gz' }
    input = {}
    input.update(s1)
    input.update(s2)
    longreadfiles=listdir("data/sequencing_genomic-longreads_trimmed/")
    longreadfiles_species=list(filter(lambda x:SPECIES in x, longreadfiles))
    if len(longreadfiles_species) > 0 :
      longreadfiles_species_path= "data/sequencing_genomic-longreads_trimmed/" + longreadfiles_species[0]
      if os.path.isfile(longreadfiles_species_path) == True :
        pacbio = {'pacbio' : longreadfiles_species_path }
        input.update(pacbio)
        return(input)
      return(input)
    return input

rule SPADES_hybrid_assembly:
  input:
    unpack(get_input_hybrid_assemblies)
  output:
    contigs=protected(expand("data/assembly_{assemblytype}/{{host}}/contigs.fasta",assemblytype='hybrid_doublefiltered')),
    scaffolds=protected(expand("data/assembly_{assemblytype}/{{host}}/scaffolds.fasta",assemblytype='hybrid_doublefiltered')),
    graph=protected(expand("data/assembly_{assemblytype}/{{host}}/assembly_graph.fastg",assemblytype='hybrid_doublefiltered')),
    graph_scaffolds=protected(expand("data/assembly_{assemblytype}/{{host}}/assembly_graph_with_scaffolds.gfa",assemblytype='hybrid_doublefiltered')),
    datasetyaml=protected(expand("data/assembly_{assemblytype}/{{host}}/input_dataset.yaml",assemblytype='hybrid_doublefiltered')),
    paramfile=protected(expand("data/assembly_{assemblytype}/{{host}}/params.txt",assemblytype='hybrid_doublefiltered'))
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
  run:
    if len(input) >2  :
       shell("spades.py {params.options} -t {threads} -m {resources.mem_gb} -1 {input.s1} -2 {input.s2} --pacbio {input.pacbio} -o {params.basedir} > {log.stdout} 2> {log.stderr}")
    else :
       shell("spades.py {params.options} -t {threads} -m {resources.mem_gb} -1 {input.s1} -2 {input.s2} -o {params.basedir} > {log.stdout} 2> {log.stderr}")

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
    "pigz --best -p {threads} --keep {input}"

############## rules for analyses after the manual curation process.
def get_single_or_hibrid_anvi_profile_for_export(wildcards):
    HOSTCODE=wildcards.hostcode
    if HOSTCODE in HOSTS:
        input = {'profile' : 'data/assembly_hybrid_doublefiltered_binningsignals_anvio/MERGED_' + HOSTCODE + '/PROFILE.db'}
        contigsdb = {'contigdb' : 'data/assembly_hybrid_doublefiltered_anvio/' + HOSTCODE + '/' + HOSTCODE + '_contigs.db'}
        input.update(contigsdb)
    elif HOSTCODE in SINGLEHOSTS:
        input = {'profile' : 'data/assembly_singles_doublefiltered_binningsignals_anvio/MERGED_' + HOSTCODE + '/PROFILE.db' }
        contigsdb = {'contigdb' : 'data/assembly_singles_doublefiltered_anvio/' + HOSTCODE + '/' + HOSTCODE + '_contigs.db'}
        input.update(contigsdb)
    return(input)


rule extract_curated_bins_from_anvio:
  input:
    unpack(get_single_or_hibrid_anvi_profile_for_export)
  output:
    dir=directory("data/curated_bins/{collection}/{hostcode}"),
    sum="data/curated_bins/{collection}/{hostcode}/bins_summary.txt" 
  params:
    collection = lambda w:expand("-C {collection}",collection=w.collection)
  conda:
    "envs/anvio6.yaml"
  log:
    stdout="logs/anvi-export-bins-{collection}-{hostcode}.stdout",
    stderr="logs/anvi-export-bins-{collection}-{hostcode}.stderr"
  shell:
    "anvi-summarize -p {input.profile} -c {input.contigdb} {params.collection} -o {output.dir} > {log.stdout} 2> {log.stderr}"

rule collect_bins_in_folder:
  input:
    "data/curated_bins/{collection}/{hostcode}"
  output:
    directory("data/curated_bins/{collection}/{hostcode}.bin-fastas")
  shell:
    """
    if [ ! -d {output} ]
    then mkdir {output}
    fi
    cp --reflink=always {input}/bin_by_bin/*/*contigs.fa {output}
    """

rule checkm_curated_bins:
  input:
    folder="data/curated_bins/{collection}/{hostcode}",
    bins="data/curated_bins/{collection}/{hostcode}.bin-fastas",
    set_root="references/checkm_data_setroot.done"
  output:
    table="data/curated_bins/{collection}/{hostcode}/{hostcode}.checkm_out.txt"
  params:
    options="-x fa --pplacer_threads=12 --tab_table",
    dir=lambda w:expand("data/curated_bins/{collection}/{hostcode}/checkm/{hostcode}",collection=w.collection,hostcode=w.hostcode)
  threads: 72
  log:
    stdout="logs/checkm_curated_bins_{collection}_{hostcode}.stdout",
    stderr="logs/checkm_curated_bins_{collection}_{hostcode}.stderr"
  conda:
    "envs/checkm.yaml"
  shell:
    """
    if [ -d {params.dir} ]
    then rm -rf {params.dir}
    fi
    checkm lineage_wf -t {threads} {params.options} {input.bins} {params.dir} -f {output.table} > {log.stdout} 2> {log.stderr}
    """

rule BAT_curated_bins:
  input:
    bindir="data/curated_bins/{collection}/{hostcode}.bin-fastas",
    dmnd="references/CAT_customised_20190108/CAT_database_customised/2019-03-27.nr.dmnd",
    db=  "references/CAT_customised_20190108/CAT_database_customised",
    tf=  "references/CAT_customised_20190108/taxonomy_customised"
  output:
    "data/curated_bins/{collection}/{hostcode}/BAT/{hostcode}.BAT.bin2classification.txt"
  params:
    options= " -s '.fa' ",
    prefix=lambda w : expand( "data/curated_bins/{collection}/{hostcode}/BAT/{hostcode}.BAT" , collection=w.collection, hostcode=w.hostcode)
  threads: 12
  conda:
    "envs/cat.yaml"
  log:
    stdout="logs/BAT_{collection}_{hostcode}.stdout",
    stderr="logs/BAT_{collection}_{hostcode}.stderr"
  shell:
    "CAT bins -n {threads} -b {input.bindir} -d {input.db} -t {input.tf} {params.options} -o {params.prefix} > {log.stdout} 2> {log.stderr}"

rule BAT_add_names_curated_bins:
  input:
    i="data/curated_bins/{collection}/{hostcode}/BAT/{hostcode}.BAT.bin2classification.txt",
    tf="references/CAT_customised_20190108/taxonomy_customised"
  output:
    "data/curated_bins/{collection}/{hostcode}/{hostcode}.BAT.names.txt"
  params:
    "--only_official"
  log:
    stdout="logs/BAT_bins-{collection}_classification_taxonomy_{hostcode}.stdout",
    stderr="logs/BAT_bins-{collection}_classification_taxonomy_{hostcode}.stderr"
  threads: 1
  conda:
    "envs/cat.yaml"
  shell:
    "CAT add_names {params} -i {input.i} -t {input.tf} -o {output} > {log.stdout} 2> {log.stderr}"

rule collect_curated_bin_info:
  input:
    anvio="data/curated_bins/{collection}/{hostcode}/bins_summary.txt" ,
    checkm= "data/curated_bins/{collection}/{hostcode}/{hostcode}.checkm_out.txt",
    BAT=    "data/curated_bins/{collection}/{hostcode}/{hostcode}.BAT.names.txt"
  output:
    tmp = temp("data/curated_bins/{collection}/bin_info.{hostcode}.tmp"),
    tab =      "data/curated_bins/{collection}/bin_info.{hostcode}.tab"
  threads: 1
  log:
    stderr="logs/collect_curated_bin_info-{collection}_{hostcode}.stderr"
  shell:
    """
    join <(tail -n +2 {input.anvio} )                                        \
         <(sed -E 's/-contigs//g' {input.checkm} | tail -n +2 | tr ' ' "_")  \
         -j 1 -t "\t"                                                        \
         > {output.tmp} 2>  {log.stderr}
    paste <(head -n1 {input.anvio}   | tr ' ' '_') \
          <(head -n1 {input.checkm}  | tr ' ' '_' | sed "s/Bin_Id\t//g") \
          <(head -n1 {input.BAT}     | tr ' ' '_' | sed "s/\#_bin\t//g"  ) \
          > {output.tab}
    join {output.tmp}                                                         \
         <(sed -E 's/-contigs\.fa//g' {input.BAT} | tail -n +2 | sed  's/\:\ [01]\.[0-9][0-9]//g' | tr ' ' "_" ) \
         -j 1 -t "\t"                                                         \
          >> {output.tab} 2>> {log.stderr}
    """

rule merge_curated_bin_info:
  input:
    expand("data/curated_bins/{{collection}}/bin_info.{hostcode}.tab", hostcode=REFINED)
  output:
    "data/curated_bins/{collection}.bin_info.tab"
  shell:
    """
    # get array in bash:
    tabs=( {input} )
    # extract header and modify:
    head -n 1 ${{tabs[0]}} | sed -E "s/^/Sample\t/"> {output}
    # get content of all other files, add sample column, and discard header
    for t in ${{tabs[@]}}
    do  name=$(echo $t | rev | cut -f 1 -d '/' | rev | cut -f 2 -d '.')
        sed -E "s/^/$name\t/" $t | tail -n +2 >> {output}
    done

    """

rule collect_binning_stats:
  input:
    expand("data/curated_bins/{collection}.bin_info.tab", collection=['refined','metabat2','CONCOCT'])
  output:
    "data/curated_bins/summary.bin_info.tab"
  shell:
    """
    # get array in bash:
    tabs=( {input} )
    # extract header and modify:
    head -n 1 ${{tabs[0]}} | sed -E "s/^/Collection\t/"> {output}
    # get content of all other files, add sample column, and discard header
    for t in ${{tabs[@]}}
    do  name=$(echo $t | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
        sed -E "s/^/$name\t/" $t | tail -n +2 >> {output}
    done

    """

rule fastq_length:
  input:
    "{folder}/{file}.fastq.gz"
  output:
    "{folder}/{file}.fastq_length"
  log:
    stderr="logs/fastq_length_{folder}_{file}.stderr"
  shell:
    """
    echo "$(zcat {input} | wc -l) / 4" | bc > {output}
    """

rule collect_fastq_lengths:
  input:
    expand("data/sequencing_genomic/{hostcode}_R1.fastq_length",hostcode=HOSTCODES)
