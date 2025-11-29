import pandas as pd
import sys
import logging
import numpy as np

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Define the path to your mapping file
# In a real Snakemake workflow, you might pass this as a parameter: snakemake.input.species_map
# SPECIES_MAP_FILE = "/home/data/CYM/pipeline/contig_to_genus_species.tsv"
SPECIES_MAP_FILE = snakemake.input.species_map
# input_files = snakemake.input.data_files

try:
    input_files = snakemake.input.data_files # If species_map is an input, it would be snakemake.input.data_files
    output_file = snakemake.output[0]
    sample_ids = snakemake.params.sample_ids

    logging.info(f"Starting indel ratio merging for {len(input_files)} samples.")
    logging.info(f"Sample IDs: {sample_ids}")

    # --- Load Species Mapping ---
    try:
        logging.info(f"Loading species mapping from: {SPECIES_MAP_FILE}")
        species_map_df = pd.read_csv(SPECIES_MAP_FILE, sep='\t', header=None, names=['ContigID', 'Species'])
        # Create a dictionary for faster lookup: ContigID -> Species
        contig_to_species_dict = pd.Series(species_map_df.Species.values, index=species_map_df.ContigID).to_dict()
        logging.info(f"Successfully loaded {len(contig_to_species_dict)} species mappings.")
    except FileNotFoundError:
        logging.error(f"Species mapping file not found: {SPECIES_MAP_FILE}")
        raise
    except Exception as e:
        logging.error(f"Error loading or processing species mapping file {SPECIES_MAP_FILE}: {e}")
        raise
    # --- End Load Species Mapping ---

    if len(input_files) != len(sample_ids):
        raise ValueError(f"Mismatch between number of input files ({len(input_files)}) and sample IDs ({len(sample_ids)})")

    all_sample_series = []

    for i, file_path in enumerate(input_files):
        sample_id = sample_ids[i]
        logging.info(f"Processing file for sample {sample_id}: {file_path}")
        try:
            temp_df = pd.read_csv(
                file_path,
                sep='\t',
                usecols=["CHR", "POS", "Depth", "Total_depth"],
                dtype={"CHR": str, "POS": str, "Depth": float, "Total_depth": float},
                comment='#'
            )

            if temp_df.empty:
                logging.warning(f"File {file_path} for sample {sample_id} is empty or contains no data rows.")
                sample_series = pd.Series(dtype=float, name=sample_id)
                all_sample_series.append(sample_series)
                continue

            temp_df['site'] = temp_df['CHR'] + ':' + temp_df['POS']
            temp_df['indel_ratio'] = 0.0
            mask_valid_total_depth = temp_df['Total_depth'] > 0
            temp_df.loc[mask_valid_total_depth, 'indel_ratio'] = \
                temp_df['Depth'][mask_valid_total_depth] / temp_df['Total_depth'][mask_valid_total_depth]
            sample_series = pd.Series(temp_df['indel_ratio'].values, index=temp_df['site'], name=sample_id)
            all_sample_series.append(sample_series)
            logging.info(f"Successfully processed {len(sample_series)} sites for sample {sample_id}")

        except FileNotFoundError:
            logging.error(f"Input file not found: {file_path}")
            raise
        except ValueError as e:
            logging.error(f"Error processing file {file_path} (likely column issue or data type): {e}")
            raise
        except Exception as e:
            logging.error(f"An unexpected error occurred processing file {file_path}: {e}")
            raise

    logging.info("Finished reading all input files. Creating combined DataFrame.")

    if not all_sample_series:
        logging.warning("No data processed from input files. Writing empty output with only Site and Species headers.")
        # Adjust header for empty output
        with open(output_file, 'w') as f:
            f.write("Site\tSpecies\n")
        sys.exit(0)

    df_combined = pd.concat(all_sample_series, axis=1)
    df_combined.fillna(0, inplace=True)
    df_combined.index.name = "Site"

    logging.info(f"Initial combined DataFrame shape (Sites x Samples): {df_combined.shape}")

    # --- Add Species Column ---
    if not df_combined.empty:
        # Extract ContigID from the 'Site' index (format: CHR:POS)
        # The index is df_combined.index, which is a pandas Index object
        df_combined['ContigID_temp'] = df_combined.index.str.split(':').str[0]
        # Map ContigID to Species using the loaded dictionary
        # .get(key, default_value) is useful to handle contigs not in the map
        df_combined['Species'] = df_combined['ContigID_temp'].map(lambda x: contig_to_species_dict.get(x, "Unknown_Species"))
        # Drop the temporary ContigID column
        df_combined.drop(columns=['ContigID_temp'], inplace=True)
        # Reorder columns to make 'Species' the first data column (it will be second in CSV after 'Site' index)
        # Get current columns, put 'Species' first, then the rest
        current_cols = [col for col in df_combined.columns if col != 'Species']
        df_combined = df_combined[['Species'] + current_cols]
        logging.info("Species column added and DataFrame reordered.")
    # --- End Add Species Column ---


    if df_combined.empty or df_combined.shape[0] == 0:
        logging.warning("Combined DataFrame is empty or has no sites. Writing header-only output.")
        header_line = df_combined.index.name if df_combined.index.name else "Site"
        header_line += "\tSpecies" # Add Species to header
        if not df_combined.columns.drop('Species', errors='ignore').empty: # Sample columns excluding 'Species'
             header_line += "\t" + "\t".join(df_combined.columns.drop('Species', errors='ignore'))
        with open(output_file, 'w') as f:
            f.write(header_line + "\n")
    else:
        # Filtering and binarization should happen AFTER species column is added and reordered
        # The site_counts logic remains the same as it operates on numerical columns
        # Ensure species column is not included in numerical operations like sum or >0
        numerical_cols = [col for col in df_combined.columns if col != 'Species']

        if not numerical_cols: # Only Site and Species columns, no sample data
            logging.warning("DataFrame has Site and Species but no sample data columns after processing. Writing current state.")
            df_filtered_final = df_combined # Keep all rows if no numerical data to filter on
        else:
            site_counts = (df_combined[numerical_cols] > 0).sum(axis=1)
            sites_to_keep = site_counts[site_counts >= 2].index
            df_filtered = df_combined.loc[sites_to_keep]
            logging.info(f"Filtered DataFrame shape (Sites x (Species + Samples)): {df_filtered.shape}")

            if df_filtered.empty:
                logging.warning("DataFrame is empty after filtering. Writing header-only output.")
                header_line = df_filtered.index.name if df_filtered.index.name else "Site"
                header_line += "\tSpecies"
                if not df_filtered.columns.drop('Species', errors='ignore').empty:
                     header_line += "\t" + "\t".join(df_filtered.columns.drop('Species', errors='ignore'))
                with open(output_file, 'w') as f:
                    f.write(header_line + "\n")
                sys.exit(0) # Exit if empty after filtering

            # Binarize only the numerical (sample) columns
            df_filtered_final = df_filtered.copy() # Create a copy to avoid SettingWithCopyWarning
            df_filtered_final[numerical_cols] = (df_filtered[numerical_cols] > 0).astype(int)


        # df_filtered_final now has:
        # Index: Site
        # Column 0: Species
        # Column 1...N: Sample_IDs (binarized)
        df_filtered_final.to_csv(output_file, sep='\t', index=True)
        # index=True writes the 'Site' index. 'Species' is the first data column.

    logging.info(f"Successfully wrote results with Species column to {output_file}")

except Exception as e:
    logging.exception("An error occurred during script execution.")
    sys.exit(1)