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
# Assumes Snakefile rule defines named inputs:
# input:
#   snp_files = expand(...),
#   taxonomy_file = "/home/data/CYM/pipeline/database/gt-pro/species_taxonomy_ext.tsv"
snp_files = snakemake.input.snp_files
taxonomy_file_path = snakemake.input.taxonomy_file
output_file = snakemake.output[0]
FIRST_COLUMN_NAME = "SiteID" # Desired name for the site/SNP ID column

logging.info(f"Starting SNP binarization, merge, and taxonomy annotation.")
logging.info(f"Input SNP files: {len(snp_files)}")
logging.info(f"Taxonomy file: {taxonomy_file_path}")
logging.info(f"Output file: {output_file}")

# --- Load Taxonomy Data ---
taxonomy_map = {}
try:
    logging.info(f"Loading taxonomy map from {taxonomy_file_path}")
    # GUT_GENOME000004    100002    d__Bacteria;p__Firmicutes_A;...
    # We need to map column 1 (0-indexed) to column 2 (0-indexed)
    tax_df = pd.read_csv(
        taxonomy_file_path,
        sep='\t',
        header=None, # No header in the provided example
        usecols=[1, 2], # Use the second and third columns
        dtype={1: str, 2: str} # Read both as strings
    )
    # Create a dictionary from the second column (as key) and third column (as value)
    taxonomy_map = pd.Series(tax_df[2].values, index=tax_df[1]).to_dict()
    logging.info(f"Successfully loaded {len(taxonomy_map)} taxonomy entries.")
    if not taxonomy_map:
        logging.warning("Taxonomy map is empty. Taxonomy column will likely be 'Unknown'.")
except FileNotFoundError:
    logging.error(f"Taxonomy file not found: {taxonomy_file_path}. Cannot add Taxonomy column.")
    sys.exit(1)
except Exception as e:
    logging.error(f"Error loading taxonomy file {taxonomy_file_path}: {e}", exc_info=True)
    # Depending on strictness, you might want to sys.exit(1) here too.
    # For now, we'll allow proceeding, but Taxonomy column will be mostly 'Unknown'.
    logging.warning("Proceeding without full taxonomy information.")


# --- Data Processing ---
all_sample_series = [] # List to hold Series for each sample

if not snp_files:
    logging.warning("No input SNP files found. Creating an empty output file with headers.")
    # Define expected headers for an empty file
    empty_headers = [FIRST_COLUMN_NAME, "Taxonomy"]
    # If snakemake.params.samples is available and contains expected sample IDs, add them
    # For simplicity, if no actual processing, we might not know sample IDs yet.
    # If sample_ids_processed was populated somehow, use it. Here, it would be empty.
    with open(output_file, 'w') as f_out:
        f_out.write("\t".join(empty_headers) + "\n")
    sys.exit(0)

sample_ids_processed = [] # To keep track of processed sample IDs for final column order

for file_path in snp_files:
    try:
        filename = os.path.basename(file_path)
        # Adjust sample_id extraction if naming is different, e.g. using a regex or params
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
            logging.error(f"Skipping file {file_path} due to missing columns or other pd.read_csv error: {ve}")
            sample_ids_processed.append(sample_id) # Record attempt
            continue

        if df.empty:
             logging.warning(f"Skipping empty file: {file_path}")
             sample_ids_processed.append(sample_id) # Record attempt
             continue

        df['ref_count'] = pd.to_numeric(df['ref_count'], errors='coerce')
        df['alt_count'] = pd.to_numeric(df['alt_count'], errors='coerce')
        df_filtered_alleles = df[df['alt_count'] > 0].copy()

        if df_filtered_alleles.empty:
            logging.info(f"No SNPs with alt_count > 0 found in sample: {sample_id}")
            sample_ids_processed.append(sample_id)
            continue

        df_filtered_alleles['snp_id'] = df_filtered_alleles['species'] + ':' + df_filtered_alleles['global_pos']
        sample_series = pd.Series(1, index=df_filtered_alleles['snp_id'], name=sample_id)
        if sample_series.index.has_duplicates:
            logging.warning(f"Duplicate snp_ids found in {file_path} for sample {sample_id}. Taking first occurrence.")
            sample_series = sample_series[~sample_series.index.duplicated(keep='first')]

        all_sample_series.append(sample_series)
        sample_ids_processed.append(sample_id)
        logging.info(f"Sample {sample_id} contributed {len(sample_series)} SNPs with alt alleles.")

    except Exception as e:
        logging.error(f"Error processing file {file_path}: {e}", exc_info=True)
        filename = os.path.basename(file_path) # Ensure filename is defined
        sample_id_error = filename.replace(".snp.tsv", "") # Ensure sample_id is defined for error case
        if sample_id_error not in sample_ids_processed:
           sample_ids_processed.append(sample_id_error)
        continue

logging.info(f"Finished reading individual files.")

# --- Convert to DataFrame ---
# Define processed_sample_columns early, sorted for consistency
processed_sample_columns = sorted(list(set(sample_ids_processed)))

if not all_sample_series:
     logging.warning("No SNPs with alternative alleles found across all files. Creating an empty output file with headers.")
     # Define expected headers for an empty file
     header_cols = [FIRST_COLUMN_NAME, "Taxonomy"]
     if processed_sample_columns: # Add sample columns if any were processed/attempted
         header_cols.extend(processed_sample_columns)
     with open(output_file, 'w') as f_out:
        f_out.write("\t".join(header_cols) + "\n")
     sys.exit(0)

logging.info("Concatenating sample series into a matrix...")
merged_df = pd.concat(all_sample_series, axis=1)
merged_df = merged_df.reindex(columns=processed_sample_columns) # Ensure all samples are columns
merged_df.fillna(0, inplace=True)
merged_df = merged_df.astype(int)
merged_df.index.name = FIRST_COLUMN_NAME
merged_df.sort_index(inplace=True) # Sort rows (SNP IDs)

logging.info(f"Initial binarized matrix shape ({FIRST_COLUMN_NAME} x Samples): {merged_df.shape}")

# --- Filtering: Keep SNPs present based on 50% rule with a minimum of 2 ---
num_samples_total = merged_df.shape[1]

if num_samples_total > 1:  # 只在样本数大于1时进行筛选
    # 计算50%样本数，并向上取整
    threshold_from_percentage = math.ceil(num_samples_total * 0.25)
    
    # 要求的最小出现样本数
    minimum_presence_count = 2
    
    # 最终的筛选阈值取上述两者中的较大值
    final_threshold = max(threshold_from_percentage, minimum_presence_count)
    
    logging.info(f"Applying filter: SNP must be present in at least {final_threshold} samples. "
                 f"(Calculated from max(ceil({num_samples_total}*0.5)={threshold_from_percentage}, min_required={minimum_presence_count})). "
                 f"Total samples: {num_samples_total}.")

    # 计算每个SNP（行）在多少个样本中出现（值为1）
    snp_presence_counts = merged_df.sum(axis=1) # axis=1 表示按行求和

    # 筛选出满足条件的SNP
    snps_to_keep_filter = snp_presence_counts >= final_threshold
    final_df_for_taxonomy = merged_df[snps_to_keep_filter]

    logging.info(f"Matrix shape before presence filtering ({FIRST_COLUMN_NAME} x Samples): {merged_df.shape}")
    logging.info(f"Matrix shape after presence filtering: {final_df_for_taxonomy.shape}")
    if final_df_for_taxonomy.shape[0] < merged_df.shape[0]:
        logging.info(f"Number of SNPs removed by presence filtering: {merged_df.shape[0] - final_df_for_taxonomy.shape[0]}")
    else:
        logging.info("No SNPs were removed by presence filtering.")
else:
    logging.info(f"Skipping presence filtering: Number of samples ({num_samples_total}) is not greater than 1.")
    final_df_for_taxonomy = merged_df # 不进行筛选

# --- Add Taxonomy Column and Prepare for Save ---
# `final_df_for_taxonomy` has SiteID as index. We need to make it a column to insert Taxonomy.
if not final_df_for_taxonomy.empty:
    df_to_save = final_df_for_taxonomy.reset_index() # SiteID becomes the first column

    # Extract species ID from SiteID (e.g., '100002' from '100002:1016131')
    # This species_id is used to lookup in taxonomy_map
    # The keys in taxonomy_map are strings like '100002'
    species_ids_for_map = df_to_save[FIRST_COLUMN_NAME].str.split(':', expand=True)[0]
    
    # Map to get taxonomy strings; fill missing with "Unknown"
    taxonomy_col_values = species_ids_for_map.map(taxonomy_map).fillna("Unknown")
    
    # Insert 'Taxonomy' column as the second column (at index 1)
    df_to_save.insert(1, 'Taxonomy', taxonomy_col_values)
    
    logging.info(f"Taxonomy column added. Final matrix shape for saving: {df_to_save.shape}")
else:
    # If final_df_for_taxonomy is empty after filtering
    logging.warning("DataFrame is empty after filtering. Creating file with headers only.")
    # Create an empty DataFrame with the correct columns for header writing
    # Use processed_sample_columns for sample names
    header_list = [FIRST_COLUMN_NAME, 'Taxonomy'] + processed_sample_columns
    df_to_save = pd.DataFrame(columns=header_list)


# --- Save Output ---
logging.info(f"Saving final matrix to {output_file}")
if df_to_save.empty and not ([FIRST_COLUMN_NAME, 'Taxonomy'] + processed_sample_columns): # Check if header list is also effectively empty
    # This case should ideally be caught earlier, but as a safeguard:
    logging.warning("Final DataFrame is empty and no column headers could be determined. Writing a minimal header.")
    with open(output_file, 'w') as f_out:
        f_out.write(f"{FIRST_COLUMN_NAME}\tTaxonomy\n")
elif df_to_save.empty: # Has columns defined, but no data rows
    logging.warning("Final DataFrame has headers but is empty (no data rows). Writing headers only.")
    df_to_save.to_csv(output_file, sep='\t', index=False, header=True)
else:
    # df_to_save has SiteID as a regular column, and Taxonomy inserted.
    # index=False because SiteID is already a column.
    df_to_save.to_csv(output_file, sep='\t', index=False, header=True)

logging.info(f"Successfully wrote final matrix to {output_file}")