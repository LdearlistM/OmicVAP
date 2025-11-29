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
        logging.warning("No data processed from input files. Writing empty output with only Site header.")
        with open(output_file, 'w') as f:
            f.write("Site\n")
        sys.exit(0)

    df_combined = pd.concat(all_sample_series, axis=1)
    df_combined.fillna(0, inplace=True)
    df_combined.index.name = "Site"

    logging.info(f"Initial combined DataFrame shape (Sites x Samples): {df_combined.shape}")

    if df_combined.empty or df_combined.shape[0] == 0:
        logging.warning("Combined DataFrame is empty or has no sites. Writing header-only output.")
        header_line = df_combined.index.name if df_combined.index.name else "Site"
        if not df_combined.columns.empty:
            header_line += "\t" + "\t".join(df_combined.columns)
        with open(output_file, 'w') as f:
            f.write(header_line + "\n")
    else:
        numerical_cols = list(df_combined.columns)

        if not numerical_cols:
            logging.warning("DataFrame has Site but no sample data columns. Writing current state.")
            df_filtered_final = df_combined
        else:
            site_counts = (df_combined[numerical_cols] > 0).sum(axis=1)
            sites_to_keep = site_counts[site_counts >= 2].index
            df_filtered = df_combined.loc[sites_to_keep]
            logging.info(f"Filtered DataFrame shape (Sites x Samples): {df_filtered.shape}")

            if df_filtered.empty:
                logging.warning("DataFrame is empty after filtering. Writing header-only output.")
                header_line = df_filtered.index.name if df_filtered.index.name else "Site"
                if not df_filtered.columns.empty:
                    header_line += "\t" + "\t".join(df_filtered.columns)
                with open(output_file, 'w') as f:
                    f.write(header_line + "\n")
                sys.exit(0)

            df_filtered_final = df_filtered.copy()
            df_filtered_final[numerical_cols] = (df_filtered[numerical_cols] > 0).astype(int)

        df_filtered_final.to_csv(output_file, sep='\t', index=True)

    logging.info(f"Successfully wrote results to {output_file}")

except Exception as e:
    logging.exception("An error occurred during script execution.")
    sys.exit(1)
