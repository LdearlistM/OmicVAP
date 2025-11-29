# scripts/merge_inversions.py
import sys
import os
import pandas as pd
import logging

# --- Configuration Section ---
input_files = snakemake.input.data_files # List of input files
output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Define the path to your mapping file (adjust if passed via Snakemake)
# SPECIES_MAP_FILE = "/home/data/CYM/pipeline/database/contig_to_genus_species.tsv"
SPECIES_MAP_FILE = snakemake.input.species_map

# Define the desired name for the first column (variant IDs)
FIRST_COLUMN_NAME = "Site" # Or "Locus", "InversionID", etc.

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

# --- Load Species Mapping ---
try:
    logging.info(f"Loading species mapping from: {SPECIES_MAP_FILE}")
    species_map_df = pd.read_csv(SPECIES_MAP_FILE, sep='\t', header=None, names=['ContigID', 'Species'])
    contig_to_species_dict = pd.Series(species_map_df.Species.values, index=species_map_df.ContigID).to_dict()
    logging.info(f"Successfully loaded {len(contig_to_species_dict)} species mappings.")
except FileNotFoundError:
    logging.error(f"Species mapping file not found: {SPECIES_MAP_FILE}")
    sys.exit(1) # Exit if critical mapping file is missing
except Exception as e:
    logging.error(f"Error loading or processing species mapping file {SPECIES_MAP_FILE}: {e}", exc_info=True)
    sys.exit(1)
# --- End Load Species Mapping ---

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
        sample_id = filename.replace(".inv.tsv", "") # Adjust if sample ID extraction is more complex
        logging.info(f"Processing: {sample_id} from {f_path}")

        try:
            # Use pandas to read the file for robustness and ease
            temp_df = pd.read_csv(f_path, sep='\t', na_values=["NA", ""], keep_default_na=True)
            
            if temp_df.empty:
                logging.warning(f"File {f_path} for sample {sample_id} resulted in an empty DataFrame after read_csv.")
                continue

            # Ensure required columns exist
            required_cols = ['ID', 'Pe_ratio', 'Span_ratio']
            if not all(col in temp_df.columns for col in required_cols):
                logging.error(f"File {f_path} is missing one or more required columns: {required_cols}. Found: {temp_df.columns.tolist()}")
                continue
            
            # Fill NA values in ratio columns with 0 after reading
            temp_df['Pe_ratio'] = temp_df['Pe_ratio'].fillna(0.0)
            temp_df['Span_ratio'] = temp_df['Span_ratio'].fillna(0.0)

            # Calculate inv_ratio (vectorized)
            temp_df['inv_ratio'] = temp_df[['Pe_ratio', 'Span_ratio']].max(axis=1)
            
            # Select and rename columns for a cleaner structure to append
            temp_df_processed = temp_df[['ID', 'inv_ratio']].copy()
            temp_df_processed['sampleID'] = sample_id
            
            all_data.append(temp_df_processed)

        except Exception as e:
            logging.error(f"Error processing file {f_path}: {e}", exc_info=True)

    if not all_data:
        logging.error("No data was successfully processed from any input file. Creating empty output.")
        # Create an empty file with just the headers
        with open(output_file, 'w') as outfile:
            outfile.write(f"{FIRST_COLUMN_NAME}\tSpecies\n") # Minimal output with new headers
    else:
        logging.info("Concatenating all processed dataframes...")
        df = pd.concat(all_data, ignore_index=True)

        if df.empty:
            logging.error("Concatenated DataFrame is empty. Creating empty output.")
            with open(output_file, 'w') as outfile:
                outfile.write(f"{FIRST_COLUMN_NAME}\tSpecies\n")
            sys.exit(0)


        logging.info("Pivoting DataFrame...")
        try:
            # Pivot: samples as index, variant IDs as columns, inv_ratio as values
            df_pivot = df.pivot_table(index='sampleID', columns='ID', values='inv_ratio', fill_value=0)
            # If there can be duplicate (sampleID, ID) pairs and you want to sum/mean them:
            # df_pivot = df.pivot_table(index='sampleID', columns='ID', values='inv_ratio', fill_value=0, aggfunc='max') # or 'sum', 'mean'
        except Exception as e:
             logging.error(f"Error pivoting data: {e}", exc_info=True)
             logging.error("\nDataFrame head before pivot attempt:\n%s", df.head())
             sys.exit(1)

        logging.info("Sorting DataFrame columns (Variant IDs)...")
        df_pivot = df_pivot.reindex(sorted(df_pivot.columns), axis=1)
        # Sorting rows (sampleIDs) is already handled by pivot or can be done: df_pivot = df_pivot.sort_index()

        # --- Binarize Data ---
        logging.info("Binarizing data (values > 0 become 1)...")
        df_binary = (df_pivot > 0).astype(int)

        # --- Transpose Data ---
        # Now Variant IDs become rows (index), and Sample IDs become columns
        df_transposed = df_binary.T
        df_transposed.index.name = FIRST_COLUMN_NAME # Rename the index (which contains variant IDs)

        # --- Add Species Column ---
        if not df_transposed.empty:
            # Extract ContigID from the 'VariantID' index (e.g., NC_009614.1 from NC_009614.1:348986-...)
            # Assuming the ContigID is the part before the first colon in the VariantID
            df_transposed['ContigID_temp'] = df_transposed.index.str.split(':').str[0]
            
            df_transposed['Species'] = df_transposed['ContigID_temp'].map(
                lambda x: contig_to_species_dict.get(x, "Unknown_Species")
            )
            df_transposed.drop(columns=['ContigID_temp'], inplace=True)
            
            # Reorder columns to make 'Species' the first data column (second in CSV)
            current_cols = [col for col in df_transposed.columns if col != 'Species']
            df_final = df_transposed[['Species'] + current_cols]
            logging.info("Species column added and DataFrame reordered.")
        else:
            # If df_transposed is empty, create a structure for header writing
            df_final = pd.DataFrame(columns=['Species'])
            df_final.index.name = FIRST_COLUMN_NAME


        logging.info(f"Writing output to {output_file}...")
        if df_final.empty and df_final.columns.tolist() == ['Species']: # Special case for no data but species column exists
             with open(output_file, 'w') as outfile:
                outfile.write(f"{FIRST_COLUMN_NAME}\tSpecies\n") # Write headers
        else:
            df_final.to_csv(output_file, sep='\t', index=True, header=True)
            # index=True writes the index (now named FIRST_COLUMN_NAME)
            # header=True writes the column names (Species, SampleID1, SampleID2...)

    logging.info("Script finished successfully.")

except Exception as e:
    logging.exception("An unexpected error occurred during script execution.")
    sys.exit(1)