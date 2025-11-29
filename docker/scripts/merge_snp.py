# scripts/merge_filter_snp_ratios_v2.py
import pandas as pd
import math
import os
import sys
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(snakemake.log[0]),
        logging.StreamHandler(sys.stdout) # Also print to stdout
    ]
)

# --- Input and Output ---
snp_files = snakemake.input.snp_files
output_file = snakemake.output[0]
FIRST_COLUMN_NAME = "SiteID"

logging.info(f"Starting SNP binarization and merge.")
logging.info(f"Input SNP files: {len(snp_files)}")
logging.info(f"Output file: {output_file}")

# --- Data Processing ---
all_sample_series = []
sample_ids_processed = []

if not snp_files:
    logging.warning("No input SNP files found. Creating an empty output file with headers.")
    with open(output_file, 'w') as f_out:
        f_out.write(FIRST_COLUMN_NAME + "\n")
    sys.exit(0)

for file_path in snp_files:
    try:
        filename = os.path.basename(file_path)
        sample_id = filename.replace(".snp.tsv", "")
        logging.info(f"Processing sample: {sample_id} from file: {file_path}")

        try:
            df = pd.read_csv(
                file_path,
                sep='\t',
                usecols=['species', 'global_pos', 'ref_count', 'alt_count'],
                dtype={'species': str, 'global_pos': str, 'ref_count': str, 'alt_count': str},
                low_memory=False
            )
        except ValueError as ve:
            logging.error(f"Skipping file {file_path} due to pd.read_csv error: {ve}")
            sample_ids_processed.append(sample_id)
            continue

        if df.empty:
            logging.warning(f"Skipping empty file: {file_path}")
            sample_ids_processed.append(sample_id)
            continue

        df['ref_count'] = pd.to_numeric(df['ref_count'], errors='coerce')
        df['alt_count'] = pd.to_numeric(df['alt_count'], errors='coerce')
        df_filtered = df[df['alt_count'] > 0].copy()

        if df_filtered.empty:
            logging.info(f"No SNPs with alt_count > 0 in sample: {sample_id}")
            sample_ids_processed.append(sample_id)
            continue

        df_filtered['snp_id'] = df_filtered['species'] + ':' + df_filtered['global_pos']
        sample_series = pd.Series(1, index=df_filtered['snp_id'], name=sample_id)

        if sample_series.index.has_duplicates:
            logging.warning(f"Duplicate snp_ids in {file_path}. Taking first occurrence.")
            sample_series = sample_series[~sample_series.index.duplicated(keep='first')]

        all_sample_series.append(sample_series)
        sample_ids_processed.append(sample_id)
        logging.info(f"Sample {sample_id} contributed {len(sample_series)} SNPs.")

    except Exception as e:
        logging.error(f"Error processing file {file_path}: {e}", exc_info=True)
        sample_id_error = os.path.basename(file_path).replace(".snp.tsv", "")
        if sample_id_error not in sample_ids_processed:
            sample_ids_processed.append(sample_id_error)
        continue

logging.info(f"Finished reading input files.")

# --- Merge Series into DataFrame ---
processed_sample_columns = sorted(list(set(sample_ids_processed)))

if not all_sample_series:
    logging.warning("No SNPs with alt alleles found. Creating empty output file.")
    header_cols = [FIRST_COLUMN_NAME] + processed_sample_columns
    with open(output_file, 'w') as f_out:
        f_out.write("\t".join(header_cols) + "\n")
    sys.exit(0)

merged_df = pd.concat(all_sample_series, axis=1)
merged_df = merged_df.reindex(columns=processed_sample_columns)
merged_df.fillna(0, inplace=True)
merged_df = merged_df.astype(int)
merged_df.index.name = FIRST_COLUMN_NAME
merged_df.sort_index(inplace=True)

logging.info(f"Initial matrix shape ({FIRST_COLUMN_NAME} x Samples): {merged_df.shape}")

# --- Filter SNPs: present in at least max(25% samples, 2 samples) ---
num_samples_total = merged_df.shape[1]

if num_samples_total > 1:
    threshold = max(math.ceil(num_samples_total * 0.25), 2)
    logging.info(f"Filtering SNPs present in >= {threshold} samples (out of {num_samples_total})")

    presence_counts = merged_df.sum(axis=1)
    snps_to_keep = presence_counts >= threshold
    final_df = merged_df[snps_to_keep]

    logging.info(f"After filtering: {final_df.shape[0]} SNPs retained (before: {merged_df.shape[0]})")
else:
    logging.info("Skipping filtering due to only one sample.")
    final_df = merged_df

# --- Save Result ---
final_df.reset_index(inplace=True)
final_df.to_csv(output_file, sep='\t', index=False)
logging.info(f"Output written to: {output_file}")
