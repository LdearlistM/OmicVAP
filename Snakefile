configfile: "config.yaml"

SAMPLES=[
	"CCMD19168690ST-11-0",
	"CCMD31134579ST-11-0",
	"CCMD41521570ST-11-0",
	"CCMD45812507ST-11-0",
	"CCMD46727384ST-11-0",
	"CCMD49461418ST-11-0",
	"CCMD51154251ST-11-0",
	"CCMD52145360ST-11-0",
	"CCMD53522274ST-11-0",
	"CCMD59540613ST-11-0",
	"CCMD59583015ST-11-0",
	"CCMD65222621ST-11-0",
	"CCMD71242853ST-11-0",
	"CCMD76222476ST-11-0",
	"CCMD85481373ST-11-0",
	"CCMD89967135ST-11-0",
	"CCMD95431029ST-11-0",
	"CCMD98531134ST-11-0",
	"CCMD79349503ST-11-0",
	"CCMD15562448ST-11-0"
]

rule all:
	input:
		# QC
		expand("01-QC/{sample}/{sample}_paired_1.fastq", sample=SAMPLES),
		expand("01-QC/{sample}/{sample}_paired_2.fastq", sample=SAMPLES),
		# SNP
		expand("02-variant-calling/SNP/{sample}/{sample}.snp.tsv", sample=SAMPLES),
		"02-variant-calling/CNV/across_sample/snps/snps_summary.tsv",
		# INDEL
		expand("02-variant-calling/INDEL/{sample}/{sample}.vcf.indel.txt", sample=SAMPLES),
		# SV
		"02-variant-calling/SV/SGVFinder2/merge/dsgv.csv",
		"02-variant-calling/SV/SGVFinder2/merge/vsgv.csv",
		"02-variant-calling/SV/SGVFinder/merge/dsgv.csv",
		"02-variant-calling/SV/SGVFinder/merge/vsgv.csv",
		expand("02-variant-calling/SV/PhaseFinder/{sample}/{sample}.inversion.txt", sample=SAMPLES),
		# CNV
		"02-variant-calling/CNV/across_sample/genes/genes_summary.tsv"

rule quality_control:
	input:
		R1 = "/home/data/CYM/CRC/rawdata/{sample}_1.fastq.gz",
		R2 = "/home/data/CYM/CRC/rawdata/{sample}_2.fastq.gz"
	output:
		"01-QC/{sample}/{sample}_paired_1.fastq",
		"01-QC/{sample}/{sample}_paired_2.fastq"
	conda:
		"envs/kneaddata-latest.yaml"
	threads: 8
	params:
		output_prefix = "{sample}",
		kneaddata_db = config["kneaddata_db"],
		trimmomatic_path = config["trimmomatic_path"],
		trimmomatic_options = config["trimmomatic_options"]
	shell:
		"kneaddata -i1 {input.R1} -i2 {input.R2} -o 01-QC/{params.output_prefix}/ -db {params.kneaddata_db} "
		"-t {threads} --output-prefix {params.output_prefix} --trimmomatic {params.trimmomatic_path} "
		"--trimmomatic-options '{params.trimmomatic_options}' "
		"--bypass-trf --reorder --remove-intermediate-output && "
		"rm 01-QC/{wildcards.sample}/{wildcards.sample}_hg_38*.fastq && "
		"rm 01-QC/{wildcards.sample}/{wildcards.sample}_unmatched*.fastq"

rule midas_run_species:
	input:
		R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
		R2 = "01-QC/{sample}/{sample}_paired_2.fastq"
	output:
		"02-variant-calling/CNV/single_sample/{sample}/species/species_profile.tsv"
	conda:
		"envs/midasv3.yaml"
	log:
		"logs/02-variant-calling/CNV/run-species/{sample}.log"
	threads: 8
	params:
		midasdb_name = config["midasdb_name"],
		midasdb_dir = config["midasv3_db"]
	shell:
		"(midas run_species --sample_name {wildcards.sample} -1 {input.R1} -2 {input.R2} "
		"--midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
		"--num_cores {threads} 02-variant-calling/CNV/single_sample/) 2> {log}"
	
rule midas_merge_species:
	input:
		expand("02-variant-calling/CNV/single_sample/{sample}/species/species_profile.tsv", sample=SAMPLES)
	output:
		"02-variant-calling/CNV/across_sample/species/species_prevalence.tsv"
	conda:
		"envs/midasv3.yaml"
	log:
		"logs/02-variant-calling/CNV/merge-species/species.log"
	threads: 4
	shell:
		"(midas merge_species --samples_list samples_list.tsv --min_cov 2 "
		"02-variant-calling/CNV/across_sample/) 2> {log}"

rule midas_run_snps:
	input:
		R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
		R2 = "01-QC/{sample}/{sample}_paired_2.fastq",
		species = "02-variant-calling/CNV/across_sample/species/species_prevalence.tsv"
	output:
		"02-variant-calling/CNV/single_sample/{sample}/snps/snps_summary.tsv"
	conda:
		"envs/midasv3.yaml"
	log:
		"logs/02-variant-calling/SNP/MIDAS/run-snps/{sample}.log"
	threads: 8
	params:
		midasdb_name = config["midasdb_name"],
		midasdb_dir = config["midasv3_db"]
	shell:
		"(midas run_snps --sample_name {wildcards.sample} -1 {input.R1} -2 {input.R2} "
		"--midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
		"--num_cores {threads} --chunk_size 1000000 --ignore_ambiguous "
		"--select_by median_marker_coverage,unique_fraction_covered --select_threshold=2,0.5 "
		"02-variant-calling/CNV/single_sample/) 2> {log}"

rule midas_merge_snps:
	input:
		expand("02-variant-calling/CNV/single_sample/{sample}/snps/snps_summary.tsv", sample=SAMPLES)
	output:
		"02-variant-calling/CNV/across_sample/snps/snps_summary.tsv"
	conda:
		"envs/midasv3.yaml"
	log:
		"logs/02-variant-calling/CNV/merge-snps/snps.log"
	threads: 8
	params:
		midasdb_name = config["midasdb_name"],
		midasdb_dir = config["midasv3_db"]
	shell:
		"(midas merge_snps --samples_list samples_list.tsv "
		"--midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
		"--num_cores {threads} --chunk_size 100000 --genome_coverage 0.7 "
		"02-variant-calling/CNV/across_sample/) 2> {log}"

rule midas_build_bowtie2db:
	input:
		"02-variant-calling/CNV/across_sample/species/species_prevalence.tsv"
	output:
		"02-variant-calling/CNV/across_sample/bt2_indexes/pangenomes.species"
	conda:
		"envs/midasv3.yaml"
	log:
		"logs/02-variant-calling/CNV/build-bowtie2db/bt2-db-build-pangenomes.log"
	threads: 8
	params:
		midasdb_name = config["midasdb_name"],
		midasdb_dir = config["midasv3_db"]
	shell:
		"(midas build_bowtie2db --midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
		"--bt2_indexes_name pangenomes --bt2_indexes_dir 02-variant-calling/CNV/across_sample/bt2_indexes "
		"--species_profile {input} --num_cores {threads} --select_by sample_counts --select_threshold 1) 2> {log}"

rule midas_run_genes:
	input:
		R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
		R2 = "01-QC/{sample}/{sample}_paired_2.fastq",
		species = "02-variant-calling/CNV/across_sample/bt2_indexes/pangenomes.species"
	output:
		"02-variant-calling/CNV/single_sample/{sample}/genes/genes_summary.tsv"
	conda:
		"envs/midasv3.yaml"
	log:
		"logs/02-variant-calling/CNV/run-genes/{sample}.log"
	threads: 8
	params:
		midasdb_name = config["midasdb_name"],
		midasdb_dir = config["midasv3_db"]
	shell:
		"(midas run_genes --sample_name {wildcards.sample} -1 {input.R1} -2 {input.R2} "
		"--midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
		"--num_cores {threads} --select_threshold=-1 "
		"--prebuilt_bowtie2_indexes 02-variant-calling/CNV/across_sample/bt2_indexes/pangenomes "
		"--prebuilt_bowtie2_species {input.species} "
		"02-variant-calling/CNV/single_sample/) 2> {log}"

rule midas_merge_genes:
	input:
		expand("02-variant-calling/CNV/single_sample/{sample}/genes/genes_summary.tsv", sample=SAMPLES)
	output:
		"02-variant-calling/CNV/across_sample/genes/genes_summary.tsv"
	conda:
		"envs/midasv3.yaml"
	log:
		"logs/02-variant-calling/CNV/merge-genes/genes.log"
	threads: 8
	params:
		midasdb_name = config["midasdb_name"],
		midasdb_dir = config["midasv3_db"]
	shell:
		"(midas merge_genes --samples_list samples_list.tsv "
		"--midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
		"--num_cores {threads} --sample_counts 2 --cluster_level_in 99 --genome_depth 0.4 "
		"02-variant-calling/CNV/across_sample/) 2> {log}"

rule GT_Pro_genotype:
	input:
		R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
		R2 = "01-QC/{sample}/{sample}_paired_2.fastq"
	output:
		"02-variant-calling/SNP/{sample}/{sample}.tsv"
	conda:
		"envs/GT-Pro.yaml"
	log:
		"logs/02-variant-calling/SNP/{sample}.log"
	threads: 8
	params:
		reference = config["GT_Pro_db"],
		GT_Pro = "/home/data/CYM/software/GT-Pro/gt-pro/GT_Pro"
	shell:
		"({params.GT_Pro} genotype -d {params.reference} {input.R1} {input.R2} -t {threads} "
		"-o 02-variant-calling/SNP/{wildcards.sample}/{wildcards.sample}) 2> {log}"

rule GT_Pro_parse:
	input:
		"02-variant-calling/SNP/{sample}/{sample}.tsv"
	output:
		"02-variant-calling/SNP/{sample}/{sample}.snp.tsv"
	conda:
		"envs/GT-Pro.yaml"
	log:
		"logs/02-variant-calling/SNP/{sample}.snp.log"
	threads: 2
	params:
		GT_Pro = "/home/data/CYM/software/GT-Pro/gt-pro/GT_Pro",
		dict_path = config["GT_dict_path"]
	shell:
		"({params.GT_Pro} parse --dict {params.dict_path} --in {input} --out {output}) 2> {log}"

rule bwa_mem:
	input:
		R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
		R2 = "01-QC/{sample}/{sample}_paired_2.fastq",
		reference = config['QuickVariant_db']
	output:
		"02-variant-calling/INDEL/{sample}/{sample}.sam"
	conda:
		"envs/bwa.yaml"
	log:
		"logs/02-variant-calling/INDEL/{sample}.sam.log"
	threads: 4
	shell:
		"(bwa mem -t {threads} {input.reference} {input.R1} {input.R2} "
		"> {output}) 2> {log}"

rule QuickVariant_vcf:
	input:
		sam = "02-variant-calling/INDEL/{sample}/{sample}.sam",
		reference = config['QuickVariant_db']
	output:
		"02-variant-calling/INDEL/{sample}/{sample}.vcf"
	conda:
		"envs/QuickVariants.yaml"
	threads: 4
	params:
		quickvariant = "/home/data/CYM/pipeline/.snakemake/conda/6d0f440ac4144d9731adccf20a121355_/bin/quick-variants-1.1.0.jar"
	shell:
		"java -jar {params.quickvariant} --out-vcf {output} "
		"--reference {input.reference} --in-sam {input.sam} "
		"--vcf-exclude-non-mutations --num-threads {threads}"

rule QuickVariant_indel:
	input:
		"02-variant-calling/INDEL/{sample}/{sample}.vcf"
	output:
		"02-variant-calling/INDEL/{sample}/{sample}.vcf.indel.txt"
	conda:
		"envs/QuickVariants.yaml"
	threads: 2
	params:
		output_dir = "/home/data/CYM/pipeline/02-variant-calling/INDEL/{sample}"
	shell:
		"python /home/data/CYM/software/QuickVariants/QuickVariants_Indelfilter.py -i {params.output_dir} && "
		"rm 02-variant-calling/INDEL/{wildcards.sample}/{wildcards.sample}.vcf.snp"

rule SGVFinder2_ICRA:
	input:
		R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
		R2 = "01-QC/{sample}/{sample}_paired_2.fastq"
	output:
		"02-variant-calling/SV/SGVFinder2/ICRA/{sample}_paired.smp"
	singularity: "./images/sgvfinder2.sif"
	threads: 2
	params:
		database = config["SGVFinder2_db"],
		output_dir = "/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder2/ICRA"
	script:
		"scripts/ICRA.py"

rule SGVFinder2_merge:
	input:
		expand("02-variant-calling/SV/SGVFinder2/ICRA/{sample}_paired.smp", sample=SAMPLES)
	output:
		"02-variant-calling/SV/SGVFinder2/merge/dsgv.csv",
		"02-variant-calling/SV/SGVFinder2/merge/vsgv.csv"
	singularity: "./images/sgvfinder2.sif"
	threads: 2
	params:
		database = config["SGVFinder2_db"],
		input_dir = "/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder2/ICRA",
		output_dir = "/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder2/merge"
	script:
		"scripts/merge.py"

rule SGVFinder_ICRA:
    input:
        R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
        R2 = "01-QC/{sample}/{sample}_paired_2.fastq"
    output:
        "02-variant-calling/SV/SGVFinder/ICRA/{sample}_paired.jsdel"
    singularity: "./images/sgvfinder.sif"
    log:
        "logs/02-variant-calling/SV/SGVFinder/ICRA/ICRA_{sample}.log"
    threads: 8
    params:
        output_dir = "/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder/ICRA",
        prefix = "01-QC/{sample}/{sample}_paired"
    shell:
        """
        source activate sgvfinder && \\
        python /opt/SGVFinder/src/ICRA_cmd.py {params.output_dir} {params.prefix} --pe 2> {log}
        """

rule SGVFinder_PerFile:
    input:
        jsdel = "02-variant-calling/SV/SGVFinder/ICRA/{sample}_paired.jsdel"
    output:
        map_jsdel = "02-variant-calling/SV/SGVFinder/ICRA/{sample}_paired.map.jsdel"
    singularity: "./images/sgvfinder.sif"
    log:
        "logs/02-variant-calling/SV/SGVFinder/SGVF_PerFile/{sample}.log"
    threads: 8
    shell:
        """
        source activate sgvfinder && \\
        python /opt/SGVFinder/src/SGVF_PerFile_cmd.py {input.jsdel} {output.map_jsdel} 100 --x_coverage 0.01 --rate_param 10 2> {log}
        """

rule SGVFinder_merge:
    input:
        expand("02-variant-calling/SV/SGVFinder/ICRA/{sample}_paired.map.jsdel", sample=SAMPLES)
    output:
        dsgv = "02-variant-calling/SV/SGVFinder/merge/dsgv.csv",
        vsgv = "02-variant-calling/SV/SGVFinder/merge/vsgv.csv"
    singularity: "./images/sgvfinder.sif"
    log:
        "logs/02-variant-calling/SV/SGVFinder/SGVF_merge/merge.log"
    threads: 8
    params:
        input_glob_string = "/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder/ICRA/*.jsdel",
        browser_path = "/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder/browser"
    shell:
        """
        source activate sgvfinder && \\
        mkdir -p {params.browser_path} && \\
        python /opt/SGVFinder/src/SGVF_cmd.py "/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder/ICRA/*.map.jsdel" {output.dsgv} {output.vsgv} --min_samp_cutoff 2 --x_coverage 0.01 --rate_param 10 --browser_path {params.browser_path} --csv_output 2> {log}
        """

rule PhaseFinder:
	input:
		reference = config['PhaseFinder_db'],
		R1 = "01-QC/{sample}/{sample}_paired_1.fastq",
		R2 = "01-QC/{sample}/{sample}_paired_2.fastq"
	output:
		"02-variant-calling/SV/PhaseFinder/{sample}/{sample}.inversion.txt"
	conda:
		"envs/PhaseFinder.yaml"
	log:
		"logs/02-variant-calling/SV/PhaseFinder/{sample}.log"
	threads: 8
	shell:
		"""
		python scripts/PhaseFinder.py ratio \
		-i {input.reference} -1 {input.R1} -2 {input.R2} -p {threads} \
		-o 02-variant-calling/SV/PhaseFinder/{wildcards.sample}/{wildcards.sample} && \
		awk '($4 != 0 && $4 != \"NA\") || ($7 != 0 && $7 != \"NA\")' 02-variant-calling/SV/PhaseFinder/{wildcards.sample}/{wildcards.sample}.ratio.txt > \
		02-variant-calling/SV/PhaseFinder/{wildcards.sample}/{wildcards.sample}.inversion.txt 2> {log}
		"""