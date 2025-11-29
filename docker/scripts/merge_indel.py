import pandas as pd
import sys
import logging
import numpy as np

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

try:
    input_files = snakemake.input.data_files
    output_file = snakemake.output[0]
    sample_ids = snakemake.params.sample_ids

    logging.info(f"Starting indel ratio merging for {len(input_files)} samples.")
    logging.info(f"Sample IDs: {sample_ids}")

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
                # ***** MODIFICATION START *****
                usecols=["CHR", "POS", "REF", "INDEL", "Depth", "Total_depth"], # Added REF and INDEL
                dtype={
                    "CHR": str,
                    "POS": str,
                    "REF": str, # Read REF as string
                    "INDEL": str, # Read INDEL as string
                    "Depth": float,
                    "Total_depth": float
                },
                # ***** MODIFICATION END *****
                comment='#'
            )

            if temp_df.empty:
                logging.warning(f"File {file_path} for sample {sample_id} is empty or contains no data rows.")
                sample_series = pd.Series(dtype=float, name=sample_id) # Will be converted to int later
                all_sample_series.append(sample_series)
                continue

            # ***** MODIFICATION START: Determine INDEL Type and create new site ID *****
            # Ensure REF and INDEL columns are treated as strings and handle NaN as empty strings for logic
            temp_df['REF'] = temp_df['REF'].fillna('').astype(str)
            temp_df['INDEL'] = temp_df['INDEL'].fillna('').astype(str)

            # Determine type: -1 for deletion, 1 for insertion
            # Conditions:
            # Deletion: REF has content, INDEL is empty
            # Insertion: REF is empty, INDEL has content
            conditions = [
                (temp_df['REF'].str.len() > 0) & (temp_df['INDEL'].str.len() == 0), # Deletion
                (temp_df['REF'].str.len() == 0) & (temp_df['INDEL'].str.len() > 0)  # Insertion
            ]
            choices = ['-1', '1'] # Type as string initially
            
            temp_df['indel_type_str'] = np.select(conditions, choices, default='0') # Default '0' for ambiguous cases or non-INDELs if any

            # Filter out rows where type could not be determined (default '0')
            # Or handle them as errors if they are unexpected
            # For now, let's assume all rows are valid INDELs as per your description
            # If a row results in '0', it means it didn't fit the deletion/insertion criteria
            # You might want to log or filter these:
            ambiguous_indels = temp_df[temp_df['indel_type_str'] == '0']
            if not ambiguous_indels.empty:
                logging.warning(f"Sample {sample_id}, file {file_path} has {len(ambiguous_indels)} rows with ambiguous INDEL type (REF/INDEL columns might be both empty or both filled). These will be excluded or marked with type '0'. Original data:\n{ambiguous_indels[['CHR', 'POS', 'REF', 'INDEL']].head()}")
                # Option: filter them out if they are not expected
                # temp_df = temp_df[temp_df['indel_type_str'] != '0']
                # If you keep them with type '0', they will have site like CHR:POS:0

            temp_df['site'] = temp_df['CHR'] + ':' + temp_df['POS'] + ':' + temp_df['indel_type_str']
            # ***** MODIFICATION END *****

            # The indel_ratio calculation logic might not be directly used if the goal is a 0/1 matrix,
            # but it doesn't hurt to keep it if it's part of an intermediate step or for other purposes.
            # If you only need presence/absence, this part could be simplified.
            temp_df['indel_presence'] = 0 # Initialize as 0 (absent)
            # Assuming 'Depth' > 0 means the INDEL is present for this specific CHR:POS:Type
            # Or, if your input file already implies presence for every row, then it's always 1.
            # Let's assume every row in the input file for a sample means that INDEL (CHR:POS:Type) is present in that sample.
            # The filtering later (df_filtered[numerical_cols] > 0).astype(int) will handle the 0/1.
            # So, we just need to associate a value (e.g., 1, or the ratio if you were using it) with the site.
            # For a 0/1 matrix, simply having the site in the series implies presence.
            # The value in the series doesn't strictly matter before the final ( > 0).astype(int) conversion,
            # as long as it's > 0 for present sites.
            # We can use 1.0 to indicate presence.
            
            sample_series = pd.Series(1.0, index=temp_df['site'], name=sample_id)
            # If a site is duplicated (e.g. CHR:POS has both an insertion and deletion reported as separate rows,
            # but after adding :Type they become distinct), this direct Series creation is fine.
            # If CHR:POS:Type could still be duplicated within a single sample file (which shouldn't happen if CHR:POS:REF:INDEL is unique),
            # you might need a groupby:
            # sample_series = temp_df.groupby('site').apply(lambda x: 1.0 if (x['Depth'] > 0).any() else 0.0)
            # But given your input format, each row is a distinct INDEL call.

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
        logging.warning("No data processed from input files. Writing empty output with only Site header.")
        with open(output_file, 'w') as f:
            f.write("Site\n") # Default if no sample IDs
            if sample_ids: # Add sample IDs if available, even if no data
                 f.write("Site\t" + "\t".join(sample_ids) + "\n")
        sys.exit(0)

    df_combined = pd.concat(all_sample_series, axis=1)
    df_combined.fillna(0, inplace=True) # Sites not present in a sample get 0
    df_combined.index.name = "Site"

    logging.info(f"Initial combined DataFrame shape (Sites x Samples): {df_combined.shape}")

    # The rest of the script for filtering and converting to 0/1 should work as before.
    # The `df_filtered_final[numerical_cols] = (df_filtered[numerical_cols] > 0).astype(int)`
    # will convert any positive value (like our 1.0 for presence) to 1, and NaN/0 to 0.

    if df_combined.empty or df_combined.shape[0] == 0:
        logging.warning("Combined DataFrame is empty or has no sites. Writing header-only output.")
        header_line = df_combined.index.name if df_combined.index.name else "Site"
        if not df_combined.columns.empty: # Use df_combined.columns, not sample_ids directly here
            header_line += "\t" + "\t".join(df_combined.columns)
        else: # If no columns (e.g. no samples processed successfully)
             header_line += "\t" + "\t".join(sample_ids) # Fallback to original sample_ids if df_combined.columns is empty
        with open(output_file, 'w') as f:
            f.write(header_line + "\n")
    else:
        numerical_cols = list(df_combined.columns)

        if not numerical_cols:
            logging.warning("DataFrame has Site but no sample data columns. Writing current state.")
            # This case means df_combined has an index (sites) but no sample columns,
            # which shouldn't happen if all_sample_series had named Series.
            # If it does, we just write the sites.
            df_combined.to_csv(output_file, sep='\t', index=True, header=False, columns=[]) # Write only index
        else:
            # Filter sites present in at least 2 samples (original logic)
            # (df_combined[numerical_cols] > 0) will be True for our 1.0 values
            site_counts = (df_combined[numerical_cols] > 0).sum(axis=1)
            sites_to_keep = site_counts[site_counts >= 2].index # Keep sites present in at least 2 samples
            
            # If you want to keep sites present in at least 1 sample, change to:
            # sites_to_keep = site_counts[site_counts >= 1].index

            df_filtered = df_combined.loc[sites_to_keep]
            logging.info(f"Filtered DataFrame shape (Sites x Samples) after keeping sites in >=2 samples: {df_filtered.shape}")

            if df_filtered.empty:
                logging.warning("DataFrame is empty after filtering (no site present in >=2 samples). Writing header-only output.")
                header_line = df_filtered.index.name if df_filtered.index.name else "Site"
                if not df_filtered.columns.empty:
                    header_line += "\t" + "\t".join(df_filtered.columns)
                else: # Fallback if df_filtered.columns is empty
                    header_line += "\t" + "\t".join(sample_ids)

                with open(output_file, 'w') as f:
                    f.write(header_line + "\n")
                # sys.exit(0) # Decide if you want to exit or write an empty data frame with headers
                # Writing an empty dataframe with headers might be better for downstream tools
                # Create an empty DataFrame with correct columns to write headers
                empty_df_for_header = pd.DataFrame(columns=df_combined.columns)
                empty_df_for_header.index.name = "Site"
                empty_df_for_header.to_csv(output_file, sep='\t', index=True)

            else:
                # Convert to 0/1 matrix
                df_filtered_final = df_filtered.copy()
                df_filtered_final[numerical_cols] = (df_filtered[numerical_cols] > 0).astype(int)
                df_filtered_final.to_csv(output_file, sep='\t', index=True)

    logging.info(f"Successfully wrote results to {output_file}")

except Exception as e:
    logging.exception("An error occurred during script execution.") # This will log the full traceback
    sys.exit(1)
