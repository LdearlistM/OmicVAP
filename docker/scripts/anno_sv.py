import pandas as pd
import pyranges as pr
import sys
import logging
import re # For parsing SGV cluster IDs

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)]
)
log = logging.getLogger(__name__)

ALL_ANNOT_COLUMNS = [
    'genome_id', 'genome_name', 'accession', 'feature_type', 'patric_id',
    'start', 'end', 'strand', 'gene', 'product', 'global_start', 'global_end'
]
FINAL_OUTPUT_ORDER = ['SGV_Cluster_ID'] + ALL_ANNOT_COLUMNS

PR_CHROM_COL = "Chromosome"
PR_START_COL = "Start"
PR_END_COL = "End"
PR_STRAND_COL = "Strand"

ORIG_ANNOT_GENOME_ID_COL = "genome_id"
ORIG_ANNOT_GENOME_NAME_COL = "genome_name"
ORIG_ANNOT_ACCESSION_COL = "accession"
ORIG_ANNOT_FEATURE_TYPE_COL = "feature_type"
ORIG_ANNOT_PATRIC_ID_COL = "patric_id"
ORIG_ANNOT_LOCAL_START_COL = "start"
ORIG_ANNOT_LOCAL_END_COL = "end"
ORIG_ANNOT_STRAND_COL = "strand"
ORIG_ANNOT_GENE_COL = "gene"
ORIG_ANNOT_PRODUCT_COL = "product"
ORIG_ANNOT_GLOBAL_START_COL = "global_start"
ORIG_ANNOT_GLOBAL_END_COL = "global_end"

def parse_sgv_cluster_id(sgv_cluster_id_str):
    """
    Parses an SGV Cluster ID string like "genome_id:start1_end1;start2_end2;..."
    Coordinates are in kb and 1-based inclusive.
    Returns the genome_id and a list of (start_bp_0based, end_bp_0based_exclusive) tuples for PyRanges.
    """
    try:
        sgv_cluster_id_str = str(sgv_cluster_id_str)
        genome_id, coord_groups_str = sgv_cluster_id_str.split(':', 1)
        parsed_intervals = []
        if not coord_groups_str.strip():
            return genome_id, []

        coord_groups = coord_groups_str.split(';')
        for group in coord_groups:
            if not group.strip(): continue
            start_kb_str, end_kb_str = group.split('_')
            start_kb = int(start_kb_str)
            end_kb = int(end_kb_str)

            # Coordinates are 1-based inclusive in kb, e.g., 421_423 means [421000, 423000]
            # For PyRanges [Start, End) (0-based start, 0-based exclusive end):
            start_bp_0based = (start_kb * 1000) - 1
            end_bp_0based_exclusive = (end_kb * 1000)

            if start_bp_0based >= end_bp_0based_exclusive:
                log.warning(f"SGV interval {group} in '{sgv_cluster_id_str}' results in invalid 0-based range [{start_bp_0based}, {end_bp_0based_exclusive}). Skipping interval.")
                continue
            parsed_intervals.append((start_bp_0based, end_bp_0based_exclusive))
        return genome_id, parsed_intervals
    except ValueError as e:
        log.warning(f"Could not parse SGV_Cluster_ID due to value error (e.g., non-integer coords): '{sgv_cluster_id_str}'. Error: {e}. Returning None for genome_id.")
        return None, []
    except Exception as e:
        log.warning(f"Unexpected error parsing SGV_Cluster_ID: '{sgv_cluster_id_str}'. Error: {e}. Returning None for genome_id.")
        return None, []

def annotate_sgv_clusters(annotation_fp, sgv_matrix_fp, output_fp):
    log.info("Starting SGV cluster annotation.")
    log.info(f"Annotation file: {annotation_fp}")
    log.info(f"SGV matrix file: {sgv_matrix_fp}")
    log.info(f"Output file: {output_fp}")

    log.info(f"Reading annotation file: {annotation_fp}")
    try:
        annot_df_dtypes = {col: str for col in ALL_ANNOT_COLUMNS if col not in [ORIG_ANNOT_LOCAL_START_COL, ORIG_ANNOT_LOCAL_END_COL, ORIG_ANNOT_GLOBAL_START_COL, ORIG_ANNOT_GLOBAL_END_COL]}
        annot_df_original = pd.read_csv(annotation_fp, sep='\t', low_memory=False, dtype=annot_df_dtypes)
        log.info(f"Annotation file loaded. Shape: {annot_df_original.shape}")
    except Exception as e:
        log.error(f"Error reading annotation file {annotation_fp}: {e}")
        raise

    missing_annot_cols = [col for col in ALL_ANNOT_COLUMNS if col not in annot_df_original.columns]
    if missing_annot_cols:
        log.error(f"Annotation file '{annotation_fp}' is missing required columns: {missing_annot_cols}")
        raise ValueError(f"Annotation file missing columns: {missing_annot_cols}")

    try:
        genomeid_to_genomename_map = annot_df_original.drop_duplicates(subset=[ORIG_ANNOT_GENOME_ID_COL]) \
                                                      [[ORIG_ANNOT_GENOME_ID_COL, ORIG_ANNOT_GENOME_NAME_COL]] \
                                                      .set_index(ORIG_ANNOT_GENOME_ID_COL)[ORIG_ANNOT_GENOME_NAME_COL].to_dict()
        log.info("Genome_id to genome_name map created.")
    except KeyError as e:
        log.error(f"Error creating genome_id to name map: Missing column '{e}'. Needed: '{ORIG_ANNOT_GENOME_ID_COL}', '{ORIG_ANNOT_GENOME_NAME_COL}'.")
        raise

    annot_df_for_pyranges = annot_df_original.copy()
    annot_df_for_pyranges.rename(columns={
        ORIG_ANNOT_GENOME_ID_COL: PR_CHROM_COL,
        ORIG_ANNOT_STRAND_COL: PR_STRAND_COL
    }, inplace=True)

    if PR_STRAND_COL not in annot_df_for_pyranges.columns:
        log.warning(f"'{PR_STRAND_COL}' (from '{ORIG_ANNOT_STRAND_COL}') column not found. Intervals will be unstranded.")
        annot_df_for_pyranges[PR_STRAND_COL] = "."

    try:
        annot_df_for_pyranges[PR_START_COL] = annot_df_original[ORIG_ANNOT_GLOBAL_START_COL].astype(int) - 1
        annot_df_for_pyranges[PR_END_COL] = annot_df_original[ORIG_ANNOT_GLOBAL_END_COL].astype(int)
    except KeyError as e:
        log.error(f"Missing original global start/end columns ('{ORIG_ANNOT_GLOBAL_START_COL}' or '{ORIG_ANNOT_GLOBAL_END_COL}') for 0-based conversion: {e}.")
        raise
    except ValueError as e:
        log.error(f"Could not convert global start/end columns to integer for 0-based conversion: {e}")
        raise

    try:
        gr_annotations = pr.PyRanges(annot_df_for_pyranges)
        log.info("Annotation data (using global coords) converted to PyRanges object.")
    except Exception as e:
        log.error(f"Error creating PyRanges object from annotation_df: {e}")
        log.error(f"Columns in df passed to PyRanges for annotations: {annot_df_for_pyranges.columns.tolist()}")
        raise

    log.info(f"Reading SGV matrix file: {sgv_matrix_fp}")
    try:
        header_df = pd.read_csv(sgv_matrix_fp, sep=',', nrows=0)
        all_column_names = header_df.columns.tolist()

        if not all_column_names:
            log.error(f"SGV matrix file '{sgv_matrix_fp}' has no header or is empty.")
            raise ValueError("Empty or no header in SGV matrix file.")

        first_col_lower = str(all_column_names[0]).lower() # Ensure first col name is string for .lower()
        # Common sample ID column names, or pandas default for unnamed index col
        sample_id_indicators = ["sample", "samples", "id", "unnamed: 0", ""]
        
        if first_col_lower in sample_id_indicators and len(all_column_names) > 1 :
            sgv_cluster_ids = all_column_names[1:]
            log.info(f"First column ('{all_column_names[0]}') assumed as sample/index. Using subsequent {len(sgv_cluster_ids)} columns as SGV IDs.")
        else:
            # Try to detect if file has only one column (and that column contains SGV IDs)
            # This check might be fragile if a file genuinely has only one SGV ID in the header.
            # A more robust way might be to check if the first column name looks like an SGV ID.
            df_check_structure = pd.read_csv(sgv_matrix_fp, sep=',', header=None, nrows=1)
            if df_check_structure.shape[1] == 1 and ":" in str(df_check_structure.iloc[0,0]): # Check if looks like an SGV ID
                 log.info("SGV matrix file appears to have SGV IDs in the first column (no header recognized as sample column).")
                 sgv_cluster_ids = pd.read_csv(sgv_matrix_fp, sep=',', header=None, usecols=[0])[0].astype(str).tolist()
            else: # Default to assuming all header columns are SGV IDs if no clear sample column
                log.info(f"Assuming all {len(all_column_names)} columns in the header are SGV Cluster IDs (first is '{all_column_names[0]}').")
                sgv_cluster_ids = all_column_names

        if not sgv_cluster_ids:
            log.error(f"No SGV Cluster IDs could be extracted from {sgv_matrix_fp}")
            raise ValueError("No SGV Cluster IDs extracted.")
        log.info(f"Loaded {len(sgv_cluster_ids)} SGV Cluster IDs.")

    except Exception as e:
        log.error(f"Error reading SGV matrix file header {sgv_matrix_fp}: {e}")
        raise

    sgv_intervals_for_pr = []
    for sgv_id_str in sgv_cluster_ids:
        genome_id, intervals = parse_sgv_cluster_id(sgv_id_str)
        if genome_id is None or not intervals:
            continue # Warning already logged in parse_sgv_cluster_id
        for start_0b, end_0b in intervals:
            sgv_intervals_for_pr.append({
                PR_CHROM_COL: genome_id,
                PR_START_COL: start_0b,
                PR_END_COL: end_0b,
                'SGV_Cluster_ID_orig': sgv_id_str
            })

    if not sgv_intervals_for_pr:
        log.warning("No SGV intervals could be parsed. Output will be empty.")
        final_empty_df = pd.DataFrame(columns=FINAL_OUTPUT_ORDER)
        final_empty_df.to_csv(output_fp, sep='\t', index=False, header=True)
        log.info("Empty annotated SGV file created as no SGVs were parsed.")
        return

    sgv_pyranges_df = pd.DataFrame(sgv_intervals_for_pr)

    try:
        gr_sgv_intervals = pr.PyRanges(sgv_pyranges_df)
        log.info("SGV interval data converted to PyRanges object.")
    except Exception as e:
        log.error(f"Error creating PyRanges object from SGV intervals_df: {e}")
        log.error(f"Columns in df passed to PyRanges for SGVs: {sgv_pyranges_df.columns.tolist()}")
        raise

    log.info("Performing interval-based join to annotate SGV intervals...")
    annotated_sgv_gr = gr_sgv_intervals.join(gr_annotations, how='left', suffix='_annot', apply_strand_suffix=False)
    annotated_sgv_result_df = annotated_sgv_gr.df
    log.info(f"SGV Interval join complete. Result shape: {annotated_sgv_result_df.shape}")

    log.info("Formatting final output for SGVs...")
    final_output_list = []

    for _, row in annotated_sgv_result_df.iterrows():
        output_row = {'SGV_Cluster_ID': row['SGV_Cluster_ID_orig']}
        
        sgv_genome_id = row[PR_CHROM_COL] # This is the genome_id of the SGV interval (from left table)
        output_row[ORIG_ANNOT_GENOME_ID_COL] = sgv_genome_id
        output_row[ORIG_ANNOT_GENOME_NAME_COL] = genomeid_to_genomename_map.get(sgv_genome_id, pd.NA)

        # If patric_id from annotation (right table) is notNA, it's a match
        feature_matched = (ORIG_ANNOT_PATRIC_ID_COL in row and pd.notna(row[ORIG_ANNOT_PATRIC_ID_COL])) or \
                          (PR_START_COL + "_annot" in row and pd.notna(row[PR_START_COL + "_annot"]))


        if feature_matched:
            # Populate ALL_ANNOT_COLUMNS from the matched annotation feature part of the row
            for col_name in ALL_ANNOT_COLUMNS:
                if col_name in [ORIG_ANNOT_GENOME_ID_COL, ORIG_ANNOT_GENOME_NAME_COL]:
                    continue

                # Get values from the annotation part of the joined row
                # Non-core PyRanges columns from gr_annotations should retain original names
                # Core PyRanges columns from gr_annotations will have '_annot' suffix
                if col_name == ORIG_ANNOT_ACCESSION_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_ACCESSION_COL, pd.NA)
                elif col_name == ORIG_ANNOT_LOCAL_START_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_LOCAL_START_COL, pd.NA)
                elif col_name == ORIG_ANNOT_LOCAL_END_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_LOCAL_END_COL, pd.NA)
                elif col_name == ORIG_ANNOT_STRAND_COL:
                    # Strand from annotation (gr_annotations) will have suffix because gr_sgv_intervals
                    # might not have a strand, or if it did, it would conflict.
                    # Since we set apply_strand_suffix=False and gr_sgv_intervals has no strand,
                    # the strand from gr_annotations should be in a column named PR_STRAND_COL.
                    output_row[col_name] = row.get(PR_STRAND_COL, pd.NA) # This is the Strand from gr_annotations
                elif col_name == ORIG_ANNOT_GLOBAL_START_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_GLOBAL_START_COL, pd.NA)
                elif col_name == ORIG_ANNOT_GLOBAL_END_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_GLOBAL_END_COL, pd.NA)
                elif col_name == ORIG_ANNOT_FEATURE_TYPE_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_FEATURE_TYPE_COL, pd.NA)
                elif col_name == ORIG_ANNOT_PATRIC_ID_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_PATRIC_ID_COL, pd.NA)
                elif col_name == ORIG_ANNOT_GENE_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_GENE_COL, pd.NA)
                elif col_name == ORIG_ANNOT_PRODUCT_COL:
                    output_row[col_name] = row.get(ORIG_ANNOT_PRODUCT_COL, pd.NA)
                else: # Should not happen if ALL_ANNOT_COLUMNS is exhaustive
                    output_row[col_name] = pd.NA
        else:
            # No feature matched this SGV interval.
            # SGV_Cluster_ID, genome_id, genome_name are already set.
            # Other annotation cols are NA.
            for col_name in ALL_ANNOT_COLUMNS:
                if col_name not in output_row: # Fill remaining annotation-specific columns with NA
                    output_row[col_name] = pd.NA
        final_output_list.append(output_row)

    final_output_df = pd.DataFrame(final_output_list)

    for col in FINAL_OUTPUT_ORDER:
        if col not in final_output_df.columns:
            log.warning(f"Final output column '{col}' was not generated. Adding it as NA.")
            final_output_df[col] = pd.NA
    final_output_df = final_output_df[FINAL_OUTPUT_ORDER]

    log.info(f"Saving annotated SGVs to: {output_fp}")
    try:
        final_output_df.to_csv(output_fp, sep='\t', index=False, header=True)
        log.info("SGV Annotation complete. Output file saved successfully.")
    except Exception as e:
        log.error(f"Error saving final SGV file {output_fp}: {e}")
        raise

if __name__ == '__main__':
    try:
        annotation_input_file = snakemake.input.annotation_file
        sgv_matrix_input_file = snakemake.input.sgv_matrix
        output_file = snakemake.output.annotated_sgv

        annotate_sgv_clusters(annotation_input_file, sgv_matrix_input_file, output_file)

    except NameError:
        log.info("Running script in standalone mode (not via Snakemake). Using hardcoded file paths for testing.")
        TEST_ANNOTATION_FILE_SGV = '/home/data/CYM/pipeline/database/43_species_features.tsv'
        TEST_SGV_MATRIX_FILE_DSGV = '/home/data/CYM/pipeline/00-20-16core/02-variant-calling/SV/SGVFinder2/across-sample/vsgv.csv'
        TEST_OUTPUT_FILE_SGV = '/home/data/CYM/pipeline/00-20-16core/02-variant-calling/SV/SGVFinder2/across-sample/vsgv_annotated_test.vFinal.tsv'
        
        annotate_sgv_clusters(TEST_ANNOTATION_FILE_SGV, TEST_SGV_MATRIX_FILE_DSGV, TEST_OUTPUT_FILE_SGV)
