# scripts/merge_filter_snp_ratios.py
import pandas as pd
import math
import os
import sys
import logging
import numpy as np # Import numpy for np.nan

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
snp_files = snakemake.input
output_file = snakemake.output[0]

logging.info(f"Starting SNP ratio merge and filter (v2: filter based on non-zero ratios).")
logging.info(f"Input files: {len(snp_files)}")
logging.info(f"Output file: {output_file}")

# --- Data Processing ---
all_snp_data = {} # Dictionary to hold data: {snp_id: {sample_id: snp_ratio}}

if not snp_files:
    logging.warning("No input SNP files found. Creating an empty output file.")
    with open(output_file, 'w') as f_out:
        f_out.write("sampleID\n") # Write header for empty file
    sys.exit(0) # Exit successfully

sample_ids = []

for file_path in snp_files:
    try:
        # Extract sampleID from file path
        # sample_id = os.path.basename(os.path.dirname(file_path))
        # 获取文件名，例如 "sample1.snp.tsv"
        filename = os.path.basename(file_path)
        sample_id = filename.replace(".snp.tsv", "")
        sample_ids.append(sample_id)
        logging.info(f"Processing sample: {sample_id} from file: {file_path}")

        df = pd.read_csv(file_path, sep='\t', low_memory=False)

        # Check for empty or invalid files
        if df.empty or not all(col in df.columns for col in ['species', 'global_pos', 'ref_count', 'alt_count']):
             logging.warning(f"Skipping empty or malformed file: {file_path}")
             continue

        # Ensure counts are numeric, coerce errors to NaN
        df['ref_count'] = pd.to_numeric(df['ref_count'], errors='coerce')
        df['alt_count'] = pd.to_numeric(df['alt_count'], errors='coerce')

        # Drop rows where counts could not be converted to numeric OR where both are NaN
        df.dropna(subset=['ref_count', 'alt_count'], how='all', inplace=True)
        # Fill potential remaining NaNs in counts with 0 before calculation
        df['ref_count'] = df['ref_count'].fillna(0)
        df['alt_count'] = df['alt_count'].fillna(0)

        # Calculate total count and snp_ratio
        df['total_count'] = df['ref_count'] + df['alt_count']

        # Calculate snp_ratio. IMPORTANT: If total_count is 0, ratio is NaN. If alt_count is 0, ratio is 0.
        # We will replace 0 ratios with NaN later before merging, per the new requirement.
        df['snp_ratio'] = df['alt_count'] / df['total_count']

        # --- NEW LOGIC: Treat ratio=0 as missing/not-a-SNP for this sample ---
        # Replace 0 with NaN. Keep original NaNs (from 0/0 division) as NaN.
        df['snp_ratio'] = df['snp_ratio'].replace(0, np.nan)

        # Create unique SNP identifier
        df['snp_id'] = df['species'].astype(str) + ':' + df['global_pos'].astype(str)

        # Populate the dictionary, only adding non-NaN ratios
        for _, row in df.iterrows():
            snp_id = row['snp_id']
            snp_ratio = row['snp_ratio']
            # Only store if the ratio is valid (not NaN, meaning alt_count > 0)
            if pd.notna(snp_ratio):
                if snp_id not in all_snp_data:
                    all_snp_data[snp_id] = {}
                all_snp_data[snp_id][sample_id] = snp_ratio

    except Exception as e:
        logging.error(f"Error processing file {file_path}: {e}", exc_info=True)
        continue

logging.info(f"Finished reading individual files. Total unique SNP positions with non-zero ratios found across all samples: {len(all_snp_data)}")

# --- Convert to DataFrame ---
if not all_snp_data:
     logging.warning("No valid SNP data (ratio > 0) found across all files. Creating an empty output file.")
     with open(output_file, 'w') as f_out:
        f_out.write("sampleID\n") # Write header for empty file
     sys.exit(0)

# Create DataFrame. Values are non-zero ratios. Missing entries will be NaN.
merged_df = pd.DataFrame.from_dict(all_snp_data, orient='index')
merged_df = merged_df.reindex(columns=sorted(list(set(sample_ids)))). T # Transpose

logging.info(f"Initial merged matrix shape (Samples x SNPs): {merged_df.shape}")
logging.info(f"Matrix contains non-zero SNP ratios. Missing values or ratios=0 are represented as NaN.")

# --- Filtering based on non-zero ratio presence ---
num_samples = merged_df.shape[0]
logging.info(f"Total number of samples: {num_samples}")

filtered_df = merged_df # Start with the merged df

if num_samples == 1:
    logging.info("Filtering for single sample: Removing columns (SNPs) that are all NaN (which includes original zeros).")
    # For a single sample, any column that has data (i.e., not NaN) should be kept.
    # Since we converted 0s to NaNs earlier, columns with data inherently have ratio > 0.
    filtered_df = merged_df.dropna(axis=1, how='all')
    logging.info(f"Single sample matrix shape after removing non-SNP columns: {filtered_df.shape}")

elif num_samples > 1:
    # Calculate threshold: 50% of samples, rounded down, minimum 1
    threshold = math.floor(num_samples * 0.5)
    threshold = max(1, threshold) # Ensure threshold is at least 1
    logging.info(f"Filtering threshold: SNP must have ratio > 0 in at least {threshold} samples.")

    # Count non-NA occurrences for each SNP (column). This now counts only samples with ratio > 0.
    snp_non_zero_counts = merged_df.notna().sum(axis=0) # axis=0 sums vertically

    # Select columns (SNPs) meeting the threshold
    snps_to_keep = snp_non_zero_counts[snp_non_zero_counts > threshold].index
    filtered_df = merged_df[snps_to_keep]
    logging.info(f"Filtered matrix shape (Samples x SNPs): {filtered_df.shape}")
    logging.info(f"Number of SNPs removed by filtering: {merged_df.shape[1] - filtered_df.shape[1]}")

else: # num_samples == 0 case (should have been caught earlier, but good practice)
    logging.info("Filtering skipped: 0 samples found.")
    # filtered_df is already empty or correctly handled

# --- Final Step: Fill remaining NaN values with 0 ---
# This applies to SNPs that passed the filter but were missing (NaN) in some specific samples.
if not filtered_df.empty:
    logging.info("Replacing remaining NA values (representing missing data for retained SNPs) with 0 for the final output.")
    filtered_df = filtered_df.fillna(0)
else:
     logging.info("Filtered DataFrame is empty. No NA values to fill.")

# --- Save Output ---
logging.info(f"Saving final matrix to {output_file}")
# Ensure index label is correct
filtered_df.to_csv(output_file, sep='\t', index=True, index_label='sampleID')

logging.info(f"Successfully wrote filtered SNP ratio matrix to {output_file}")