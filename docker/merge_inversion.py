# scripts/merge_inversions.py
import sys
import os
import pandas as pd
import logging

# --- Configuration Section ---
input_files = snakemake.input.data_files  # List of input files
output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Define the desired name for the first column (variant IDs)
FIRST_COLUMN_NAME = "Site"  # Or "Locus", "InversionID", etc.

# --- Configure Logging ---
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# --- Script Start ---
logging.info("Starting inversion merging script.")
logging.info(f"Input files count: {len(input_files)}")
logging.info(f"Output file: {output_file}")
logging.info(f"Log file: {log_file}")

# --- Helper Function ---
def parse_ratio(value_str):
    if value_str == "NA" or value_str == "":
        return 0.0
    try:
        return float(value_str)
    except ValueError:
        logging.warning(f"Could not convert value '{value_str}' to float. Treating as 0.0.")
        return 0.0

# --- Main Processing Logic ---
all_data = []

try:
    for f_path in input_files:
        if not os.path.exists(f_path):
            logging.warning(f"Input file not found, skipping: {f_path}")
            continue
        if os.path.getsize(f_path) == 0:
            logging.warning(f"Input file is empty, skipping: {f_path}")
            continue

        filename = os.path.basename(f_path)
        sample_id = filename.replace(".inv.tsv", "")  # Adjust if sample ID extraction is more complex
        logging.info(f"Processing: {sample_id} from {f_path}")

        try:
            temp_df = pd.read_csv(f_path, sep='\t', na_values=["NA", ""], keep_default_na=True)

            if temp_df.empty:
                logging.warning(f"File {f_path} for sample {sample_id} resulted in an empty DataFrame after read_csv.")
                continue

            required_cols = ['ID', 'Pe_ratio', 'Span_ratio']
            if not all(col in temp_df.columns for col in required_cols):
                logging.error(f"File {f_path} is missing one or more required columns: {required_cols}. Found: {temp_df.columns.tolist()}")
                continue

            temp_df['Pe_ratio'] = temp_df['Pe_ratio'].fillna(0.0)
            temp_df['Span_ratio'] = temp_df['Span_ratio'].fillna(0.0)
            temp_df['inv_ratio'] = temp_df[['Pe_ratio', 'Span_ratio']].max(axis=1)

            temp_df_processed = temp_df[['ID', 'inv_ratio']].copy()
            temp_df_processed['sampleID'] = sample_id

            all_data.append(temp_df_processed)

        except Exception as e:
            logging.error(f"Error processing file {f_path}: {e}", exc_info=True)

    if not all_data:
        logging.error("No data was successfully processed from any input file. Creating empty output.")
        with open(output_file, 'w') as outfile:
            outfile.write(f"{FIRST_COLUMN_NAME}\n")
    else:
        logging.info("Concatenating all processed dataframes...")
        df = pd.concat(all_data, ignore_index=True)

        if df.empty:
            logging.error("Concatenated DataFrame is empty. Creating empty output.")
            with open(output_file, 'w') as outfile:
                outfile.write(f"{FIRST_COLUMN_NAME}\n")
            sys.exit(0)

        logging.info("Pivoting DataFrame...")
        try:
            df_pivot = df.pivot_table(index='sampleID', columns='ID', values='inv_ratio', fill_value=0)
        except Exception as e:
            logging.error(f"Error pivoting data: {e}", exc_info=True)
            logging.error("\nDataFrame head before pivot attempt:\n%s", df.head())
            sys.exit(1)

        logging.info("Sorting DataFrame columns (Variant IDs)...")
        df_pivot = df_pivot.reindex(sorted(df_pivot.columns), axis=1)

        logging.info("Binarizing data (values > 0 become 1)...")
        df_binary = (df_pivot > 0).astype(int)

        df_transposed = df_binary.T
        df_transposed.index.name = FIRST_COLUMN_NAME

        df_final = df_transposed

        logging.info(f"Writing output to {output_file}...")
        if df_final.empty:
            with open(output_file, 'w') as outfile:
                outfile.write(f"{FIRST_COLUMN_NAME}\n")  # Write only header
        else:
            df_final.to_csv(output_file, sep='\t', index=True, header=True)

    logging.info("Script finished successfully.")

except Exception as e:
    logging.exception("An unexpected error occurred during script execution.")
    sys.exit(1)
