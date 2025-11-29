# scripts/merge_inversions.py
import sys
import os
import pandas as pd
import logging  # Import the logging library

# --- Configuration Section ---
# Access Snakemake variables (available when run via 'script:')
input_files = snakemake.input
output_file = snakemake.output[0]
log_file = snakemake.log[0]

# --- Configure Logging ---
# Setup logging to write to the file specified in the Snakemake rule
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,  # Set the logging level (e.g., INFO, DEBUG, WARNING)
    format='%(asctime)s [%(levelname)s] %(message)s', # Define log message format
    datefmt='%Y-%m-%d %H:%M:%S' # Define date format
)
# Optional: To also see logs on console during development/debugging, uncomment next line
# logging.getLogger().addHandler(logging.StreamHandler(sys.stderr))

# --- Script Start ---
logging.info("Starting inversion merging script.")
logging.info(f"Input files count: {len(input_files)}")
logging.info(f"Output file: {output_file}")
logging.info(f"Log file: {log_file}")

# --- Helper Function ---
def parse_ratio(value_str):
    """Safely converts 'NA' or empty strings to 0.0, otherwise converts to float."""
    if value_str == "NA" or value_str == "":
        return 0.0
    try:
        return float(value_str)
    except ValueError:
        # Use logging.warning for warnings, no need for file=sys.stderr
        logging.warning(f"Could not convert value '{value_str}' to float. Treating as 0.0.")
        return 0.0

# --- Main Processing Logic ---
all_data = [] # List to store dictionaries, one per row read

try: # Wrap main logic in try...except for better error logging
    for f_path in input_files:
        if not os.path.exists(f_path):
            logging.warning(f"Input file not found, skipping: {f_path}")
            continue
        if os.path.getsize(f_path) == 0:
            logging.warning(f"Input file is empty, skipping: {f_path}")
            continue

        filename = os.path.basename(f_path)
        # Handle potential complex sample names if needed, this is basic:
        sample_id = filename.replace(".inv.tsv", "")
        logging.info(f"Processing: {sample_id} from {f_path}")

        try:
            with open(f_path, 'r') as infile:
                header = next(infile).strip().split('\t') # Read and skip header
                expected_cols = 7 # Based on the header description

                for i, line in enumerate(infile):
                    fields = line.strip().split('\t')
                    if len(fields) < expected_cols:
                        logging.warning(f"Skipping malformed line {i+2} in {f_path}. Expected >= {expected_cols} fields, got {len(fields)}. Line: '{line.strip()}'")
                        continue

                    variant_id = fields[0]
                    pe_ratio = parse_ratio(fields[3])   # Column 4 (index 3)
                    span_ratio = parse_ratio(fields[6]) # Column 7 (index 6)

                    inv_ratio = max(pe_ratio, span_ratio)

                    all_data.append({
                        'sampleID': sample_id,
                        'ID': variant_id,
                        'inv_ratio': inv_ratio
                    })

        except Exception as e:
            # Log errors related to specific file processing
            logging.error(f"Error processing file {f_path}: {e}", exc_info=True) # exc_info=True adds traceback

    if not all_data:
        logging.error("No data was successfully processed from any input file. Creating empty output.")
        # Create an empty file or a file with just the header
        with open(output_file, 'w') as outfile:
            outfile.write("sampleID\n") # Minimal output
        # No need to exit here if you want the script to finish "successfully" according to snakemake
        # sys.exit(0) # If you want to exit early

    else: # Only proceed if data exists
        logging.info("Creating DataFrame...")
        df = pd.DataFrame(all_data)

        logging.info("Pivoting DataFrame...")
        try:
            df_pivot = df.pivot_table(index='sampleID', columns='ID', values='inv_ratio', fill_value=0, aggfunc='first')
        except Exception as e:
             # Log error during pivoting and exit
             logging.error(f"Error pivoting data: {e}", exc_info=True)
             logging.error("\nDataFrame head before pivot attempt:\n%s", df.head())
             logging.error("\nDataFrame info before pivot attempt:\n%s", df.info())
             sys.exit(1) # Exit with error status

        logging.info("Sorting DataFrame columns and rows...")
        # Sort columns (IDs) alphabetically for consistent output
        df_pivot = df_pivot.reindex(sorted(df_pivot.columns), axis=1)
        # Sort rows (sampleIDs) alphabetically for consistent output
        df_pivot = df_pivot.sort_index()

        logging.info(f"Writing output to {output_file}...")
        df_pivot.T.to_csv(output_file, sep='\t', index=True, header=True, index_label='sampleID')

    logging.info("Script finished successfully.")

except Exception as e:
    # Catch any unexpected errors in the main script logic
    logging.exception("An unexpected error occurred during script execution.") # Logs traceback
    sys.exit(1) # Exit with a non-zero status code to indicate failure