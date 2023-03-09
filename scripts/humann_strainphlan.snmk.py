import pandas as pd
import glob

#SAMPLES
SAMPLES = list(pd.read_csv("df_for_thresholds.csv").iloc[:,0])

def get_file_names(wildcards):
    #basically makes the checkpoint run before this function; so that we have the output (sgbs) defined once the checkpoint is done.
    ck_output = checkpoints.generate_wildcard.get(**wildcards).output[0] 
    SGB, = glob_wildcards(os.path.join(ck_output, "{sgb}.txt"))
    return expand("db_markers/{sgb}.fna", sgb=SGB) + expand("strainphlan_output/{sgb}/{sgb}.StrainPhlAn4_concatenated.aln", sgb=SGB) + expand("strainphlan_output/{sgb}/Table_strainsharing_{sgb}.csv", sgb=SGB) + expand("strainphlan_output/{sgb}/Table_thresholds_{sgb}.csv", sgb=SGB)

def final_rule(wildcards):
    #same as the earlier function "get_file_names"
    #Just to make sure the final rule starts running at last and not straight after the checkpoint
    ck_output = checkpoints.generate_wildcard.get(**wildcards).output[0] 
    SGB, = glob_wildcards(os.path.join(ck_output, "{sgb}.txt"))
    return expand(rules.threshold_with_info.output.tabout, sgb = SGB) + expand(rules.threshold_with_info.output.strainsh, sgb = SGB)

rule all:
	input:
		expand("Humann_merged_gene_families_cpm_stratified.txt"),
		expand("Metaphlan_merged_abundance_table_SGB.txt"),
		expand("Consensus_markers/{sample}_mp4sam.pkl", sample=SAMPLES),
		expand("Clades/print_clades_only.tsv"),
		final_rule,
		expand("Table_strainsharing.txt")

rule fastp:
	input:
		fw_reads="/path/to/sample/{sample}_R1.fastq.gz",
		rv_reads="/path/to/sample/{sample}_R2.fastq.gz"
	output:
		fw_trimmed_reads=temp("{sample}_1_trim.fq.gz"),
		rv_trimmed_reads=temp("{sample}_2_trim.fq.gz"),
		stats_1="logs/{sample}/original_reads.out",
		fastp_out="logs/{sample}/fastp.out",
		json="logs/{sample}/filtering.json",
		html="logs/{sample}/filtering.html"
	threads: 16
	shell:
		"""
		seqkit stats -b {input.fw_reads} > {output.stats_1} && 
		fastp -i {input.fw_reads} -o {output.fw_trimmed_reads} -I {input.rv_reads} \
		-O {output.rv_trimmed_reads} --detect_adapter_for_pe -f 5 -r -W 4 -M 15 -l 70 \
		-h {output.html} -j {output.json} \
		2>&1 | tee {output.fastp_out}
		"""

rule human_mapping:
	input:
		fw_trimmed_reads=rules.fastp.output.fw_trimmed_reads,
		rv_trimmed_reads=rules.fastp.output.rv_trimmed_reads
	output:
		stats_2="logs/{sample}/trimmed_reads.out",
		bowtie_stats="logs/{sample}/bowtie.out",
		sam=temp("{sample}.sam"),
		bam=temp("{sample}.bam"),
		bam_unmapped=temp("{sample}_unmapped.bam"),
		bam_unm_sort=temp("{sample}_unmapped_sort.bam")
	threads: 16
	shell:
		"""
		seqkit stats -b {input.fw_trimmed_reads} > {output.stats_2} && 
		bowtie2 -x ~/bowtie2_indexes/hg19/hg19 \
		-1 {input.fw_trimmed_reads} -2 {input.rv_trimmed_reads} -p {threads} -S {output.sam} \
		--very-sensitive --dovetail --seed 128 2>&1 | tee {output.bowtie_stats} && 
		samtools view -b {output.sam} > {output.bam} && 
		samtools view -b -f 12 -F 256 {output.bam} > {output.bam_unmapped} && 
		sambamba sort -n -m 40G -o {output.bam_unm_sort} {output.bam_unmapped}
		"""

rule read_counting:
	input:
		bam_unm_sort=rules.human_mapping.output.bam_unm_sort
	output:
		fw_filtered=temp("{sample}_1_filt.fq.gz"),
		rv_filtered=temp("{sample}_2_filt.fq.gz")
	threads: 16
	shell:
		"""
		bedtools bamtofastq -i {input.bam_unm_sort} -fq {output.fw_filtered} -fq2 {output.rv_filtered}
		"""

rule subsampling:
	input:
		fw_filtered=rules.read_counting.output.fw_filtered,
		rv_filtered=rules.read_counting.output.rv_filtered
	output:
		fw_subsampled=temp("{sample}_1_subsample.fq.gz"),
		rv_subsampled=temp("{sample}_2_subsample.fq.gz"),
		stats_3="logs/{sample}/filtered_reads.out",
		concatenated=temp("{sample}.fq")
	threads: 16
	shell:
		"""
		seqkit stats -b {input.fw_filtered} > {output.stats_3} && 
		seqtk sample -s128 {input.fw_filtered} 20000000 > {output.fw_subsampled} && 
		seqtk sample -s128 {input.rv_filtered} 20000000 > {output.rv_subsampled} && 
		cat {output.fw_subsampled} {output.rv_subsampled} > {output.concatenated}
		"""

rule metaphlan:
	input:
		concatenated=rules.subsampling.output.concatenated
	output:
		stats_4="logs/{sample}/concatenated_reads.out",
		output_mp4=touch("metaphlan_output/{sample}_metaphlan/profiled_metagenome_{sample}.txt"),
		bt2out="metaphlan_output/{sample}_metaphlan/metagenome.bowtie2.bz2",
		sambz2="metaphlan_output/{sample}_metaphlan/{sample}_mp4sam.sam.bz2"
	threads: 16
	shell:
		"""
		seqkit stats -b {input.concatenated} > {output.stats_4} && 
		metaphlan {input.concatenated} --bowtie2out {output.bt2out} -s {output.sambz2} --nproc 5 \
		--unclassified_estimation -t rel_ab --input_type fastq -o {output.output_mp4} && 
		touch {output.output_mp4}
		"""

rule humann:
	input:
		fq=rules.subsampling.output.concatenated,
		output_mp4=rules.metaphlan.output.output_mp4
	output:
		out_dir=directory("humann_output/{sample}_humann"),
		path=touch("humann_output/{sample}_humann/{sample}_pathabundance.tsv"),
		gene=touch("humann_output/{sample}_humann/{sample}_genefamilies.tsv"),
		path_temp="humann_output/pathabundance/{sample}_pathabundance.tsv",
		gene_temp="humann_output/genefamilies/{sample}_genefamilies.tsv"
	threads: 64
	shell:
		"""
		humann --input {input.fq} --threads {threads} --memory-use maximum -v --search-mode uniref90 \
		--input-format fastq --output {output.out_dir} --o-log logs/{wildcards.sample}/humann3.log --taxonomic-profile {input.output_mp4} &&
		touch {output.path} && touch {output.gene} &&
		rm -r {output.out_dir}/{wildcards.sample}_humann_temp/ &&
		cat {output.path} > {output.path_temp} &&
		cat {output.gene} > {output.gene_temp}
		"""

rule humann_final_tables:
	input:
		path=expand(rules.humann.output.path, sample=SAMPLES),
		gene=expand(rules.humann.output.gene, sample=SAMPLES)
	output:
		mtab_path=temp("Humann_merged_pathway_abundance.txt"),
		mtab_gene=temp("Humann_merged_gene_families.txt"),
		mtab_path_cpm=temp("Humann_merged_path_abundance_cpm.txt"),
		mtab_gene_cpm=temp("Humann_merged_gene_families_cpm.txt"),
		tab_path=touch("Humann_merged_path_abundance_cpm_stratified.txt"),
		tab_gene=touch("Humann_merged_gene_families_cpm_stratified.txt")
	params:
		map_path="humann_output/pathabundance/",
		map_gene="humann_output/genefamilies/"
	threads: 16
	shell:
		"""
		humann_join_tables --input {params.map_path} --output {output.mtab_path} &&
		humann_join_tables --input {params.map_gene} --output {output.mtab_gene} &&
		humann_renorm_table --input {output.mtab_path} --units cpm --output {output.mtab_path} &&
		humann_renorm_table --input {output.mtab_gene} --units cpm --output {output.mtab_gene} &&
		cat {output.mtab_path} | awk '{{if(NR==1); gsub("_Abundance", ""); print; }}' > {output.mtab_path_cpm} &&
		cat {output.mtab_gene} | awk '{{if(NR==1); gsub("_Abundance-RPKs", ""); print; }}' > {output.mtab_gene_cpm} &&
		humann_split_stratified_table --input {output.mtab_path_cpm} --output ./ &&
		humann_split_stratified_table --input {output.mtab_gene_cpm} --output ./ &&
		touch {output.tab_gene} && touch {output.tab_path}
		"""

rule mp4_table:
	input:
		mp4=expand(rules.metaphlan.output.output_mp4, sample=SAMPLES)
	output:
		bugs_joined="metaphlan_output/Bugs_list_joined.tsv",
		merged_tab="Metaphlan_merged_abundance_table_SGB.txt"
	threads: 16
	shell:
		"""
		merge_metaphlan_tables.py {input.mp4} > {output.bugs_joined} && 
		grep -E "(t__)|(^clade_name)|UNCLASSIFIED" {output.bugs_joined} | sed "s/profiled_metagenome_//g" | sed "s/^.*t__//g" > {output.merged_tab}
		"""

rule consensus_markers:
	input:
		sambz2=rules.metaphlan.output.sambz2
	output:
		pkl="Consensus_markers/{sample}_mp4sam.pkl"
	threads: 1
	params:
		cm="Consensus_markers"
	shell:
		"""
		sample2markers.py -i {input.sambz2} -o {params.cm} -n {threads}
		"""

rule get_sgblist:
	input:
		pkl=expand(rules.consensus_markers.output.pkl, sample=SAMPLES)
	output:
		ctsv=touch("Clades/print_clades_only.tsv")
	params:
		cladesdir=directory("Clades"),
	threads: 1
	shell:
		"""
		strainphlan -s {input.pkl}  --mutation_rates --sample_with_n_markers 20 --marker_in_n_samples 50 --sample_with_n_markers_after_filt 10 \
		--print_clades_only -o {params.cladesdir} && 
		touch {output.ctsv}
		"""

checkpoint generate_wildcard: 
	input:
		sgb_list=rules.get_sgblist.output.ctsv
	output:
		sgbsdir=directory("Clades/sgb_list")
	shell:
		"""
		mkdir {output.sgbsdir} && 
		cat {input.sgb_list} | cut -f1 | awk "NR>1" | awk '{{fname=$0; print > "{output.sgbsdir}""/"fname".txt"; close("{output.sgbsdir}""/"fname".txt")}}'
		"""

rule extract_consensus_markers:
	output:
		fna="db_markers/{sgb}.fna"
	params:
		mapmarkers="db_markers",
		sgbs=get_file_names
	threads: 6
	shell:
		"""
		extract_markers.py -c {wildcards.sgb} -o {params.mapmarkers}
		"""

rule strainphlan:
	input:
		dbm=rules.extract_consensus_markers.output.fna
	output:
		tre="strainphlan_output/{sgb}/RAxML_bestTree.{sgb}.StrainPhlAn4.tre",
		info="strainphlan_output/{sgb}/{sgb}.info",
		aln="strainphlan_output/{sgb}/{sgb}.StrainPhlAn4_concatenated.aln"
	params:
		strainp="strainphlan_output/{sgb}/",
		extra=lambda wildcards, input: expand("Consensus_markers/{sample}_mp4sam.pkl", sample=SAMPLES),
		sgbs=get_file_names
	threads: 2
	shell:
		"""
		strainphlan -s {params.extra} -m {input.dbm} -o {params.strainp} -n {threads} -c {wildcards.sgb} --mutation_rates \
		--sample_with_n_markers 20 --marker_in_n_samples 50 --sample_with_n_markers_after_filt 10
		"""

rule ngd:
	input:
		tre=rules.strainphlan.output.tre
	output:
		dist="strainphlan_output/{sgb}/{sgb}_nGD.tsv"	
	params:
		sgbs=get_file_names
	threads: 64
	shell:
		"""
		python scripts/pyphlan/tree_pairwisedists.py -n {input.tre} {output.dist}
		"""

rule threshold_with_info:
	input:
		df_thres="df_for_thresholds.csv",
		sgb_dat=rules.mp4_table.output.bugs_joined,
		info=rules.strainphlan.output.info,
		aln=rules.strainphlan.output.aln,
		dist=rules.ngd.output.dist
	output:
		ggp="nGD_plots/{sgb}_nGD_plot.svg",
		tabout="strainphlan_output/{sgb}/Table_thresholds_{sgb}.csv",
		strainsh="strainphlan_output/{sgb}/Table_strainsharing_{sgb}.csv"
	params:
		sgbs=get_file_names
	shell:
		"""
		Rscript scripts/strainphlan_calc_threshold.R {input.df_thres} {input.sgb_dat} {input.info} {input.aln} {input.dist} \
		{output.ggp} {output.tabout} {output.strainsh}
		"""

def get_wildcards(wildcards):
    #same as the earlier function "get_file_names"
    #Just to make sure the final rule starts running at last and not straight after the checkpoint
    ck_output = checkpoints.generate_wildcard.get(**wildcards).output[0] 
    SGB, = glob_wildcards(os.path.join(ck_output, "{sgb}.txt"))
    return SGB

rule concat_tables_all_SGBs:
	input:
		final_rule
	output:
		concat_s="Table_strainsharing.txt",
		concat_t="Table_thresholds.txt"
	params:
		thres=lambda wildcards, input: expand(rules.threshold_with_info.output.tabout, sgb = get_wildcards(wildcards)),
		strain=lambda wildcards, input: expand(rules.threshold_with_info.output.strainsh, sgb = get_wildcards(wildcards))
	threads: 64
	shell:
		"""
		python scripts/concat_strainsharing_thresholds.py -t {params.thres} -s {params.strain} \
		-o {output.concat_s} -p {output.concat_t}
		"""
