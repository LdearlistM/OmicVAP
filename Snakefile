import os
import pandas as pd

configfile: "config.yaml"

def get_results_dir():
    return config.get("results_dir", "results")

try:
    SAMPLES_DF = pd.read_csv(config["samples_tsv"], sep='\t').set_index("sample", drop=False)
except FileNotFoundError:
    raise FileNotFoundError(f"Sample sheet '{config['samples_tsv']}' not found. Please create it.")
SAMPLES = list(SAMPLES_DF.index)

def get_raw_read_path(wildcards):
    read_type_map = {
        "paired_1": "r1",
        "paired_2": "r2",
        "se": "se"
    }
    wildcard_read_type = wildcards.read_type
    try:
        column_name_in_df = read_type_map[wildcard_read_type]
    except KeyError:
        raise ValueError(f"Unknown read_type wildcard '{wildcard_read_type}' in rule link_raw_reads. Expected one of {list(read_type_map.keys())}.")
    try:
        return SAMPLES_DF.loc[wildcards.sample, column_name_in_df]
    except KeyError:
        raise FileNotFoundError(
            f"Could not find entry for sample '{wildcards.sample}' "
            f"with read type '{column_name_in_df}' in the sample sheet."
        )

rule all:
    input:
        os.path.join(get_results_dir(), "02-variant-calling", "SNP", "GT-Pro", "across-sample", "snp.tsv"),
        os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "snps", "snps_summary.tsv"),
        os.path.join(get_results_dir(), "logs","MIDASv3","snp_done.txt"),
        os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "across-sample", "indel_anno.tsv"),
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "dsgv_anno.tsv"),
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "vsgv_anno.tsv"),
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "across-sample", "inversion_anno.tsv"),
        os.path.join(get_results_dir(), "logs","MIDASv3","cnv_done.txt"),
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "merge.genes_copynum.tsv"),

rule link_raw_reads:
    input:
        get_raw_read_path
    output:
        linked_read = os.path.join(get_results_dir(), "linked_reads", "{sample}_{read_type,paired_1|paired_2|se}.fastq.gz")
    shell:
        """
        mkdir -p $(dirname "{output.linked_read}")
        if [[ "{input}" == *.gz ]]; then
            cp -s "$(realpath {input})" "{output.linked_read}"
        else
            gzip -c "{input}" > "{output.linked_read}"
        fi
        """

rule quality_control_pe:
    input:
        r1 = os.path.join(get_results_dir(), "linked_reads", "{sample}_paired_1.fastq.gz"),
        r2 = os.path.join(get_results_dir(), "linked_reads", "{sample}_paired_2.fastq.gz")
    output:
        p1 = temp(os.path.join(get_results_dir(), "01-QC", "{sample}", "{sample}_paired_1.fastq")),
        p2 = temp(os.path.join(get_results_dir(), "01-QC", "{sample}", "{sample}_paired_2.fastq"))
    conda: "envs/kneaddata-latest.yaml"
    log: os.path.join(get_results_dir(), "logs", "kneaddata", "{sample}_pe.log")
    threads: 8
    params:
        output_dir = lambda wildcards: os.path.join(get_results_dir(), "01-QC", wildcards.sample),
        kneaddata_db = config["kneaddata_db"],
        trimmomatic_path = config["trimmomatic_path"],
        trimmomatic_options = config.get("trimmomatic_options", "SLIDINGWINDOW:4:20 MINLEN:50 LEADING:3 TRAILING:3"),
        contaminant_prefix = config.get("contaminant_db_prefix", "hg_38")
    shell:
        """
        (
            kneaddata -i1 {input.r1} -i2 {input.r2} -o {params.output_dir} -db {params.kneaddata_db} \
                    -t {threads} --output-prefix {wildcards.sample} \
                    --trimmomatic {params.trimmomatic_path} \
                    --trimmomatic-options "{params.trimmomatic_options}" \
                    --bypass-trf --reorder --remove-intermediate-output
        ) > {log} 2>&1
        (
            rm -f {params.output_dir}/{wildcards.sample}_{params.contaminant_prefix}*.fastq || true ; \
            rm -f {params.output_dir}/{wildcards.sample}_unmatched*.fastq || true
        ) >> {log} 2>&1
        """

rule quality_control_se:
    input:
        se = os.path.join(get_results_dir(), "linked_reads", "{sample}_se.fastq.gz")
    output:
        s1 = temp(os.path.join(get_results_dir(), "01-QC", "{sample}", "{sample}_se.fastq"))
    conda: "envs/kneaddata-latest.yaml"
    log: os.path.join(get_results_dir(), "logs", "kneaddata", "{sample}_se.log")
    threads: 8
    params:
        output_dir = lambda wildcards: os.path.join(get_results_dir(), "01-QC", wildcards.sample),
        kneaddata_db = config["kneaddata_db"],
        trimmomatic_path = config["trimmomatic_path"],
        trimmomatic_options = config.get("trimmomatic_options", "SLIDINGWINDOW:4:20 MINLEN:50 LEADING:3 TRAILING:3"),
        contaminant_prefix = config.get("contaminant_db_prefix", "hg_38")
    shell:
        """
        (
            kneaddata --unpaired {input.se} -o {params.output_dir} -db {params.kneaddata_db} \
                    -t {threads} --output-prefix {wildcards.sample} \
                    --trimmomatic {params.trimmomatic_path} \
                    --trimmomatic-options "{params.trimmomatic_options}" \
                    --bypass-trf --reorder --remove-intermediate-output
        ) > {log} 2>&1
        (
            mv {params.output_dir}/{wildcards.sample}.fastq {output.s1} ; \
            rm -f {params.output_dir}/{wildcards.sample}_{params.contaminant_prefix}*.fastq || true
        ) >> {log} 2>&1
        """

def get_reads_for_analysis(wildcards):
    if config.get("skip_qc", False):
        if config["sequencing_type"] == "PE":
            return {
                "r1": os.path.join(get_results_dir(), "linked_reads", f"{wildcards.sample}_paired_1.fastq.gz"),
                "r2": os.path.join(get_results_dir(), "linked_reads", f"{wildcards.sample}_paired_2.fastq.gz")
            }
        else:
            return {
                "se": os.path.join(get_results_dir(), "linked_reads", f"{wildcards.sample}_se.fastq.gz")
            }
    else:
        if config["sequencing_type"] == "PE":
            return {
                "r1": os.path.join(get_results_dir(), "01-QC", wildcards.sample, f"{wildcards.sample}_paired_1.fastq"),
                "r2": os.path.join(get_results_dir(), "01-QC", wildcards.sample, f"{wildcards.sample}_paired_2.fastq")
            }
        else:
            return {
                "se": os.path.join(get_results_dir(), "01-QC", wildcards.sample, f"{wildcards.sample}_se.fastq")
            }

rule midas_run_species_snp:
    input:
        unpack(get_reads_for_analysis)
    output:
        profile=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample", "{sample}", "species", "species_profile.tsv"),
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "run-species", "{sample}.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "MIDAS", "single_sample", "{sample}.tsv")
    params:
        midasdb_name=config["midasdb_name"],
        midasdb_dir=config["midasv3_db"],
        output_dir=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample"),
        input_reads_params = lambda wildcards, input: f"-1 {input.r1} -2 {input.r2}" if config["sequencing_type"] == "PE" else f"-1 {input.se}"
    shell:
        """
        (
        midas run_species --sample_name {wildcards.sample} \
            {params.input_reads_params} \
            --midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} \
            --num_cores {threads} {params.output_dir}
        ) > {log} 2>&1
        """

rule generate_samples_list_snp:
    input:
        expand(os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample", "{sample}", "species", "species_profile.tsv"),sample=SAMPLES)
    output:
        samples_list=temp(os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "samples_list.tsv"))
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "MIDAS", "samples_list.tsv")
    params:
        midas_single_sample_base=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample"),
        header="sample_name\tmidas_outdir"
    threads: 1
    run:
        with open(output.samples_list, "w") as f_out:
            f_out.write(params.header + "\n")
            for sample_id in SAMPLES:
                midas_outdir = os.path.join(params.midas_single_sample_base)
                f_out.write(f"{sample_id}\t{midas_outdir}\n")

rule midas_merge_species_snp:
    input:
        samples_list=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "samples_list.tsv")
    output:
        os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "species", "species_prevalence.tsv")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "merge-species", "species.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "MIDAS", "across_sample", "merge_species.tsv")
    params:
        output_dir=lambda wildcards, output: os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample")
    shell:
        "(midas merge_species --samples_list {input.samples_list} --min_cov 2 "
        "{params.output_dir}) > {log} 2>&1"

rule midas_run_snps:
    input:
        unpack(get_reads_for_analysis),
        species=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "species", "species_prevalence.tsv")
    output:
        os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample", "{sample}", "snps", "snps_summary.tsv")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "run-snps", "{sample}.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "MIDAS", "single_sample", "{sample}_snps.tsv")
    params:
        midasdb_name=config["midasdb_name"],
        midasdb_dir=config["midasv3_db"],
        output_dir=lambda wildcards: os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample"),
        input_reads_params = lambda wildcards, input: f"-1 {input.r1} -2 {input.r2}" if config["sequencing_type"] == "PE" else f"-1 {input.se}"
    shell:
        """
        (
        midas run_snps --sample_name {wildcards.sample} \
            {params.input_reads_params} \
            --midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} \
            --num_cores {threads} --chunk_size 1000000 --ignore_ambiguous \
            --select_by median_marker_coverage,unique_fraction_covered --select_threshold=2,0.5 \
            {params.output_dir}
        ) > {log} 2>&1
        """

rule midas_merge_snps:
    input:
        expand(os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample", "{sample}", "snps", "snps_summary.tsv"),sample=SAMPLES),
        samples_list=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "samples_list.tsv")
    output:
        os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "snps", "snps_summary.tsv")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "merge-snps", "snps.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "MIDAS", "across_sample", "merge_snps.tsv")
    params:
        midasdb_name=config["midasdb_name"],
        midasdb_dir=config["midasv3_db"],
        output_dir=lambda wildcards, output: os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample"),
        dir_to_delete=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample")
    shell:
        "(midas merge_snps --samples_list {input.samples_list} "
        "--midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
        "--num_cores {threads} --chunk_size 100000 --genome_coverage 0.7 "
        "{params.output_dir}) > {log} 2>&1 "

rule finalize_snp_analysis:
    input:
        snps_summary=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "snps", "snps_summary.tsv")
    output:
        merged_freqs=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "merge.snps_freqs.tsv"),
        merged_info=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "merge.snps_info.tsv"),
        done_file=touch(os.path.join(get_results_dir(), "logs", "MIDASv3", "snp_done.txt"))
    params:
        snps_dir=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "snps"),
        dir_to_delete_1=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "single_sample"),
        dir_to_delete_2=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "MIDAS", "across_sample", "temp")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "finalize_snp.log")
    threads: 1
    shell:
        """
        (
            python -c '
import pandas as pd
import lz4.frame
import os
import glob
import sys

snps_dir = r"{params.snps_dir}"
output_freqs_file = r"{output.merged_freqs}"
output_info_file = r"{output.merged_info}"

try:
    # --- Part 1: Merge snps_freqs files ---
    freqs_lz4_files = glob.glob(os.path.join(snps_dir, "*", "*.snps_freqs.tsv.lz4"))

    if not freqs_lz4_files:
        print("No snps_freqs.tsv.lz4 files found to merge.", file=sys.stderr)
        with open(output_freqs_file, "w") as f:
            f.write("site_id\\n")
    else:
        freqs_dfs = []
        for lz4_file in freqs_lz4_files:
            try:
                with lz4.frame.open(lz4_file, "rt", encoding="utf-8") as f:
                    df = pd.read_csv(f, sep="\\t")
                    if not df.empty:
                        df.set_index(df.columns[0], inplace=True)
                        freqs_dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not process freqs file {{lz4_file}}. Error: {{e}}", file=sys.stderr)

        if freqs_dfs:
            final_freqs_df = pd.concat(freqs_dfs, axis=0, sort=False).groupby(level=0).sum().fillna(-1.0)
            final_freqs_df.reset_index(inplace=True)
            final_freqs_df.to_csv(output_freqs_file, sep="\\t", index=False)
            print(f"Successfully merged {{len(freqs_dfs)}} snps_freqs files.")
        else:
            print("Failed to read any snps_freqs files.", file=sys.stderr)
            with open(output_freqs_file, "w") as f:
                f.write("site_id\\n")

    # --- Part 2: Merge snps_info files ---
    info_lz4_files = glob.glob(os.path.join(snps_dir, "*", "*.snps_info.tsv.lz4"))
    
    if not info_lz4_files:
        print("No snps_info.tsv.lz4 files found to merge.", file=sys.stderr)
        with open(output_info_file, "w") as f:
            f.write("site_id\\tmajor_allele\\tminor_allele\\tsample_counts\\tsnp_type\\trc_A\\trc_C\\trc_G\\trc_T\\tsc_A\\tsc_C\\tsc_G\\tsc_T\\tlocus_type\\tgene_id\\tsite_type\\tamino_acids\\n")
    else:
        info_dfs_generator = (pd.read_csv(lz4.frame.open(f, "rt", encoding="utf-8"), sep="\\t") for f in info_lz4_files)
        
        final_info_df = pd.concat(info_dfs_generator, ignore_index=True)
        
        final_info_df.drop_duplicates(inplace=True)
        
        final_info_df.to_csv(output_info_file, sep="\\t", index=False)
        print(f"Successfully merged {{len(info_lz4_files)}} snps_info files.")

except Exception as e:
    print(f"An unexpected error occurred: {{e}}", file=sys.stderr)
    sys.exit(1)
' > {log} 2>&1
        ) && rm -rf {params.dir_to_delete_1} {params.dir_to_delete_2}
        """




rule midas_run_species_cnv:
    input:
        unpack(get_reads_for_analysis)
    output:
        profile=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample", "{sample}", "species", "species_profile.tsv")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "run-species", "{sample}.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "CNV", "single_sample", "{sample}.tsv")
    params:
        midasdb_name=config["midasdb_name"],
        midasdb_dir=config["midasv3_db"],
        output_dir=lambda wildcards: os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample"),
        input_reads_params = lambda wildcards, input: f"-1 {input.r1} -2 {input.r2}" if config["sequencing_type"] == "PE" else f"-1 {input.se}"
    shell:
        """
        (
        midas run_species --sample_name {wildcards.sample} \
            {params.input_reads_params} \
            --midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} \
            --num_cores {threads} {params.output_dir}
        ) > {log} 2>&1
        """

rule generate_samples_list_cnv:
    input:
        expand(os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample", "{sample}", "species", "species_profile.tsv"), sample=SAMPLES)
    output:
        samples_list=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "samples_list.tsv"))
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "CNV", "samples_list.tsv")
    params:
        midas_single_sample_base=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample"),
        header="sample_name\tmidas_outdir"
    threads: 1
    run:
        with open(output.samples_list, "w") as f_out:
            f_out.write(params.header + "\n")
            for sample_id in SAMPLES:
                midas_outdir = os.path.join(params.midas_single_sample_base)
                f_out.write(f"{sample_id}\t{midas_outdir}\n")

rule midas_merge_species_cnv:
    input:
        samples_list=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "samples_list.tsv")
    output:
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "species", "species_prevalence.tsv")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "merge-species", "species.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "CNV", "across_sample", "merge_species.tsv")
    params:
        output_dir=lambda wildcards, output: os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample")
    shell:
        "(midas merge_species --samples_list {input.samples_list} --min_cov 2 "
        "{params.output_dir}) > {log} 2>&1"

rule midas_build_bowtie2db:
    input:
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "species", "species_prevalence.tsv")
    output:
        index_marker=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "bt2_indexes", "pangenomes.species")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "build-bowtie2db", "bt2-db-build-pangenomes.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "CNV", "across_sample", "build_bowtie2db.tsv")
    params:
        midasdb_name=config["midasdb_name"],
        midasdb_dir=config["midasv3_db"],
        bt2_indexes_dir=lambda wildcards, output: os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "bt2_indexes")
    shell:
        "(midas build_bowtie2db --midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
        "--bt2_indexes_name pangenomes --bt2_indexes_dir {params.bt2_indexes_dir} "
        "--species_profile {input} --num_cores {threads} --select_by sample_counts --select_threshold 1) > {log} 2>&1"

rule midas_run_genes:
    input:
        unpack(get_reads_for_analysis),
        species=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "bt2_indexes", "pangenomes.species")
    output:
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample", "{sample}", "genes", "genes_summary.tsv")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "run-genes", "{sample}.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "CNV", "single_sample", "{sample}_genes.tsv")
    params:
        midasdb_name=config["midasdb_name"],
        midasdb_dir=config["midasv3_db"],
        output_dir=lambda wildcards: os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample"),
        prebuilt_indexes_base=lambda wildcards: os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "bt2_indexes", "pangenomes"),
        input_reads_params = lambda wildcards, input: f"-1 {input.r1} -2 {input.r2}" if config["sequencing_type"] == "PE" else f"-1 {input.se}"
    shell:
        """
        (
        midas run_genes --sample_name {wildcards.sample} \
            {params.input_reads_params} \
            --midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} \
            --num_cores {threads} --select_threshold=-1 \
            --prebuilt_bowtie2_indexes {params.prebuilt_indexes_base} \
            --prebuilt_bowtie2_species {input.species} \
            {params.output_dir}
        ) > {log} 2>&1
        """

rule midas_merge_genes:
    input:
        expand(os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample", "{sample}", "genes", "genes_summary.tsv"), sample=SAMPLES),
        samples_list=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "samples_list.tsv")
    output:
        os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "genes", "genes_summary.tsv")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "merge-genes", "genes.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "CNV", "across_sample", "merge_genes.tsv")
    params:
        midasdb_name=config["midasdb_name"],
        midasdb_dir=config["midasv3_db"],
        output_dir=lambda wildcards, output: os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample"),
        dir_to_delete_1=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample"),
        dir_to_delete_2=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "bt2_indexes")
    shell:
        "(midas merge_genes --samples_list {input.samples_list} "
        "--midasdb_name {params.midasdb_name} --midasdb_dir {params.midasdb_dir} "
        "--num_cores {threads} --sample_counts 2 --cluster_level_in 99 --genome_depth 0.4 "
        "{params.output_dir}) > {log} 2>&1 "

rule finalize_cnv_analysis:
    input:
        genes_summary=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "genes", "genes_summary.tsv")
    output:
        merged_tsv=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "merge.genes_copynum.tsv"),
        done_file=touch(os.path.join(get_results_dir(), "logs", "MIDASv3", "cnv_done.txt"))
    params:
        genes_dir=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "genes"),
        dir_to_delete_1=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "single_sample"),
        dir_to_delete_2=os.path.join(get_results_dir(), "02-variant-calling", "SV", "CNV", "across_sample", "bt2_indexes")
    conda:
        "envs/midasv3.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "MIDASv3", "finalize_cnv.log")
    threads: 1
    shell:
        """
        (
            python -c '
import pandas as pd
import lz4.frame
import os
import glob
import sys

genes_dir = r"{params.genes_dir}"
output_file = r"{output.merged_tsv}"

try:
    lz4_files = glob.glob(os.path.join(genes_dir, "*", "*.genes_copynum.tsv.lz4"))

    if not lz4_files:
        print("No .lz4 files found to merge.", file=sys.stderr)
        with open(output_file, "w") as f:
            f.write("cluster_id\\n")
        sys.exit(0)

    all_dfs = []
    for lz4_file in lz4_files:
        try:
            with lz4.frame.open(lz4_file, "rt", encoding="utf-8") as f:
                df = pd.read_csv(f, sep="\\t")
                if not df.empty:
                    df.set_index(df.columns[0], inplace=True)
                    all_dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not process file {{lz4_file}}. Error: {{e}}", file=sys.stderr)

    if not all_dfs:
        print("Failed to read any .lz4 files.", file=sys.stderr)
        with open(output_file, "w") as f:
            f.write("cluster_id\\n")
        sys.exit(0)

    final_df = pd.concat(all_dfs, axis=0, sort=False).groupby(level=0).sum().fillna(0.000000)
    final_df.reset_index(inplace=True)
    final_df.to_csv(output_file, sep="\\t", index=False)

    print(f"Successfully merged {{len(all_dfs)}} files into the output file.")

except Exception as e:
    print(f"An unexpected error occurred: {{e}}", file=sys.stderr)
    sys.exit(1)
' > {log} 2>&1
        ) && rm -rf {params.dir_to_delete_1} {params.dir_to_delete_2}
        """



def get_gtpro_input_params(wildcards, input):
    is_pe = config["sequencing_type"] == "PE"
    is_compressed = config.get("skip_qc", False)

    if is_pe:
        r1 = f"<(gunzip -c {input.r1}) " if is_compressed else f"{input.r1} "
        r2 = f"<(gunzip -c {input.r2})" if is_compressed else f"{input.r2}"
        return r1 + r2
    else:
        return f"<(gunzip -c {input.se})" if is_compressed else f"{input.se}"

rule GT_Pro_genotype:
    input:
        unpack(get_reads_for_analysis),
    output:
        genotype_tsv=temp(os.path.join(get_results_dir(), "02-variant-calling", "SNP", "GT-Pro", "single-sample", "{sample}.tsv"))
    conda:
        "envs/GT-Pro.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "GT-Pro", "single-sample","{sample}.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "GT-Pro", "single-sample", "{sample}.tsv")
    params:
        GT_Pro="/pipeline/software/GT-Pro/gt-pro/GT_Pro",
        output_prefix=lambda wildcards, output: os.path.splitext(output.genotype_tsv)[0],
        reference=config["GT_Pro_db"],
        input_reads_str = get_gtpro_input_params
    shell:
        """
        (
        mkdir -p $(dirname {params.output_prefix}) &&
        {params.GT_Pro} genotype -d {params.reference} \
        {params.input_reads_str} \
        -t {threads} \
        -o {params.output_prefix}
        ) > {log} 2>&1
        """

rule GT_Pro_parse:
    input:
        genotype_tsv=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "GT-Pro", "single-sample", "{sample}.tsv")
    output:
        snp_tsv=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "GT-Pro", "single-sample", "{sample}.snp.tsv")
    conda:
        "envs/GT-Pro.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "GT-Pro", "single-sample","{sample}.snp.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "GT-Pro", "single-sample", "{sample}.snp.tsv")
    params:
        GT_Pro="/pipeline/software/GT-Pro/gt-pro/GT_Pro", 
        dict_path=config["GT_dict_path"]
    shell:
        "({params.GT_Pro} parse --dict {params.dict_path} --in {input.genotype_tsv} --out {output.snp_tsv}) > {log} 2>&1" 

rule GT_Pro_merge:
    input:
        snp_files = expand(os.path.join(get_results_dir(), "02-variant-calling", "SNP", "GT-Pro", "single-sample", "{sample}.snp.tsv"), sample=SAMPLES),
    output:
        merged_snp=os.path.join(get_results_dir(), "02-variant-calling", "SNP", "GT-Pro", "across-sample", "snp.tsv")
    log:
        os.path.join(get_results_dir(), "logs", "GT-Pro", "across-sample","merge_snp.log")
    conda:
        "envs/QuickVariants.yaml" 
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SNP", "GT-Pro", "across-sample", "merge_snp.tsv")
    script:
        "scripts/merge_snp.py" 

rule bwa_mem_pe:
    input:
        r1=lambda wildcards: get_reads_for_analysis(wildcards).get('r1'),
        r2=lambda wildcards: get_reads_for_analysis(wildcards).get('r2'),
        reference=config['QuickVariant_db']
    output:
        sam=temp(os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.sam"))
    conda:
        "envs/bwa.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "QuickVariant", "single-sample","{sample}_pe.sam.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "INDEL", "single-sample", "{sample}_pe.sam.tsv")
    shell:
        "(bwa mem -t {threads} {input.reference} {input.r1} {input.r2} "
        "> {output.sam}) > {log} 2>&1"

rule bwa_mem_se:
    input:
        se=lambda wildcards: get_reads_for_analysis(wildcards).get('se'),
        reference=config['QuickVariant_db']
    output:
        sam=temp(os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.sam"))
    conda:
        "envs/bwa.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "QuickVariant", "single-sample","{sample}_se.sam.log")
    threads: 8
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "INDEL", "single-sample", "{sample}_se.sam.tsv")
    shell:
        "(bwa mem -t {threads} {input.reference} {input.se} "
        "> {output.sam}) > {log} 2>&1"

rule QuickVariant_vcf:
    input:
        sam=os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.sam"),
        reference=config['QuickVariant_db']
    output:
        vcf=temp(os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.vcf"))
    conda:
        "envs/QuickVariants.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "QuickVariant", "single-sample","{sample}.vcf.log")
    threads: 4
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "INDEL", "single-sample", "{sample}.vcf.tsv")
    params:
        quickvariant_jar="/home/data/CYM/pipeline/.snakemake/conda/6d0f440ac4144d9731adccf20a121355_/bin/quick-variants-1.1.0.jar"
    shell:
        "(java -jar {params.quickvariant_jar} --out-vcf {output.vcf} "
        "--reference {input.reference} --in-sam {input.sam} "
        "--vcf-exclude-non-mutations --num-threads {threads}) > {log} 2>&1"

rule QuickVariant_indel:
    input:
        vcf=expand(os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.vcf"), sample=SAMPLES)
    output:
        indel_tsv=expand(os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.vcf.indel.txt"), sample=SAMPLES)
    conda:
        "envs/QuickVariants.yaml"
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "INDEL", "single-sample", "indelfilter.tsv")
    log:
        os.path.join(get_results_dir(), "logs", "QuickVariant", "single-sample","indelfilter.log")
    params:
        output_dir=os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample"),
        indel_script="/pipeline/software/QuickVariants/QuickVariants_Indelfilter.py" 
    shell:
        "(python {params.indel_script} -i {params.output_dir}) > {log} 2>&1"

rule QuickVariant_mv:
    input:
        os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.vcf.indel.txt")
    output:
        indel_tsv=os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.indel.tsv")
    conda:
        "envs/QuickVariants.yaml"
    threads: 1
    params:
        snp_to_remove=lambda wildcards, output: os.path.join(os.path.dirname(output.indel_tsv), f"{wildcards.sample}.vcf.snp"),
    shell:
        "mv {input} {output.indel_tsv} && "
        "rm -f {params.snp_to_remove}"

rule QuickVariant_merge:
    input:
        data_files=expand(os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "single-sample", "{sample}.indel.tsv"), sample=SAMPLES),
    output:
        merged_indel=os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "across-sample", "indel.tsv")
    log:
        os.path.join(get_results_dir(), "logs", "QuickVariant", "across-sample","merge_indels.log")
    conda:
        "envs/QuickVariants.yaml"
    params: 
        sample_ids=SAMPLES
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "INDEL", "across-sample", "merge_indels.tsv")
    script:
        "scripts/merge_indel.py"

rule QuickVariant_anno:
    input:
        merged_indel=os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "across-sample", "indel.tsv"),
        annotation_file=config["annotation_file"]
    output:
        anno_indel=os.path.join(get_results_dir(), "02-variant-calling", "INDEL", "across-sample", "indel_anno.tsv")
    log:
        os.path.join(get_results_dir(), "logs", "QuickVariant", "across-sample","anno_indels.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "INDEL", "across-sample", "anno_indels.tsv")
    script:
        "scripts/anno_indel.py"

rule SGVFinder2_ICRA_se:
    input:
        se = lambda wildcards: get_reads_for_analysis(wildcards).get('se'),
        db_index_file = config["SGVFinder2_db"] + ".1.bt2"
    output:
        smp_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}.smp")),
        pmp_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_se.pmp")),
        jspi_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_se.jspi")),
        jsdel_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_se.jsdel")),
        bam_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_se.bam"))
    conda: "envs/SGVFinder2.yaml"
    log: os.path.join(get_results_dir(), "logs", "dSV-vSV", "ICRA", "{sample}_se.log")
    threads: 8
    benchmark: os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_se.smp.tsv")
    params:
        fq1 = lambda wildcards, input: input.se,
        fq2 = "none",  
        output_dir=lambda wildcards, output: os.path.dirname(output.smp_file)
    script:
        "scripts/ICRA.py"

rule SGVFinder2_ICRA_pe:
    input:
        r1 = lambda wildcards: get_reads_for_analysis(wildcards).get('r1'),
        r2 = lambda wildcards: get_reads_for_analysis(wildcards).get('r2'),
        db_index_file = config["SGVFinder2_db"] + ".1.bt2"
    output:
        smp_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}.smp")),
        pmp_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_paired.pmp")),
        jspi_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_paired.jspi")),
        jsdel_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_paired.jsdel")),
        bam_file=temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_paired.bam"))
    conda: "envs/SGVFinder2.yaml"
    log: os.path.join(get_results_dir(), "logs", "dSV-vSV", "ICRA", "{sample}_pe.log")
    threads: 8
    benchmark: os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}_pe.smp.tsv")
    params:
        fq1 = lambda wildcards, input: input.r1,
        fq2 = lambda wildcards, input: input.r2,
        output_dir=lambda wildcards, output: os.path.dirname(output.smp_file)
    script:
        "scripts/ICRA.py"

rule SGVFinder2_merge:
    input:
        expand(os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA", "{sample}.smp"), sample=SAMPLES)
    output:
        dsgv=os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "dsgv.csv"),
        vsgv=os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "vsgv.csv")
    conda:
        "envs/SGVFinder2.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "dSV-vSV", "across-sample","merge.log")
    threads: 2
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "dSV-vSV", "across-sample", "merge.tsv")
    params:
        database=config["SGVFinder2_db"],
        input_dir=lambda wildcards, input: os.path.commonpath([os.path.dirname(f) for f in input]) if input else os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "ICRA"),
        output_dir=lambda wildcards, output: os.path.dirname(output.dsgv),
        frames_path=os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "single-species")
    script:
        "scripts/merge_sv.py"

rule SGVFinder2_anno_dsgv:
    input:
        sgv_matrix=os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "dsgv.csv"),
        annotation_file=config["annotation_file"]
    output:
        annotated_sgv=os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "dsgv_anno.tsv")
    log:
        dsgv_log=os.path.join(get_results_dir(), "logs", "dSV-vSV", "across-sample","anno_dsgv.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "dSV-vSV", "across-sample", "anno_dsgv.tsv")
    script:
        "scripts/anno_sv.py"

rule SGVFinder2_anno_vsgv:
    input:
        sgv_matrix=os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "vsgv.csv"),
        annotation_file=config["annotation_file"]
    output:
        annotated_sgv=os.path.join(get_results_dir(), "02-variant-calling", "SV", "dSV-vSV", "across-sample", "vsgv_anno.tsv")
    log:
        vsgv_log=os.path.join(get_results_dir(), "logs", "dSV-vSV", "across-sample","anno_vsgv.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "dSV-vSV", "across-sample", "anno_vsgv.tsv")
    script:
        "scripts/anno_sv.py"

rule decompress_for_phasefinder:
    input:
        r1 = os.path.join(get_results_dir(), "linked_reads", "{sample}_paired_1.fastq.gz"),
        r2 = os.path.join(get_results_dir(), "linked_reads", "{sample}_paired_2.fastq.gz")
    output:
        r1 = temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "temp_unzipped", "{sample}_paired_1.fastq")),
        r2 = temp(os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "temp_unzipped", "{sample}_paired_2.fastq"))
    log:
        os.path.join(get_results_dir(), "logs", "decompress", "{sample}.log")
    shell:
        """
        (
        gunzip -c {input.r1} > {output.r1}
        gunzip -c {input.r2} > {output.r2}
        ) > {log} 2>&1
        """

def get_reads_for_phasefinder(wildcards):
    if config.get("skip_qc", False):
        return {
            "r1": os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "temp_unzipped", f"{wildcards.sample}_paired_1.fastq"),
            "r2": os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "temp_unzipped", f"{wildcards.sample}_paired_2.fastq")
        }
    else:
        return {
            "r1": os.path.join(get_results_dir(), "01-QC", wildcards.sample, f"{wildcards.sample}_paired_1.fastq"),
            "r2": os.path.join(get_results_dir(), "01-QC", wildcards.sample, f"{wildcards.sample}_paired_2.fastq")
        }

rule PhaseFinder:
    input:
        unpack(get_reads_for_phasefinder)
    output:
        inv_tsv=os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "single-sample", "{sample}.inv.tsv")
    params:
        PhaseFinder_db=config['PhaseFinder_db'],
        output_prefix=lambda wildcards: os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "single-sample", wildcards.sample),
        ratio_file=lambda wildcards: f"{os.path.join(get_results_dir(), '02-variant-calling', 'SV', 'Inversion', 'single-sample', wildcards.sample)}.ratio.txt",
    log:
        os.path.join(get_results_dir(), "logs", "Inversion", "single-sample", "{sample}.log")
    conda: "envs/PhaseFinder.yaml"
    shell:
        """
        (
            python {workflow.basedir}/scripts/PhaseFinder.py ratio \
            -i {params.PhaseFinder_db} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} -o {params.output_prefix}
        ) && \
        (
            awk '($4 != 0 && $4 != "NA") || ($7 != 0 && $7 != "NA")' {params.ratio_file} > {output.inv_tsv}
        ) && \
        (
            rm -f {params.ratio_file}
        )
        """

rule PhaseFinder_merge:
    input:
        data_files=expand(os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "single-sample", "{sample}.inv.tsv"), sample=SAMPLES),
    output:
        merged_inv=os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "across-sample", "inversion.tsv")
    conda:
        "envs/PhaseFinder.yaml"
    log:
        os.path.join(get_results_dir(), "logs", "Inversion", "across-sample","merge_inv.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "Inversion", "across-sample", "merge_inv.tsv")
    script:
        "scripts/merge_inversion.py"

rule PhaseFinder_anno:
    input:
        merged_inversion=os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "across-sample", "inversion.tsv"),
        annotation_file=config["annotation_file"]
    output:
        anno_inversion=os.path.join(get_results_dir(), "02-variant-calling", "SV", "Inversion", "across-sample", "inversion_anno.tsv")
    log:
        os.path.join(get_results_dir(), "logs", "Inversion", "across-sample","anno_inv.log")
    threads: 1
    benchmark:
        os.path.join(get_results_dir(), "benchmarks", "02-variant-calling", "SV", "Inversion", "across-sample", "anno_inversions.tsv")
    script:
        "scripts/anno_inversion.py"
