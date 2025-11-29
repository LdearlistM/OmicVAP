import pandas as pd
import pyranges as pr
import sys
import logging

# --- 0. 设置日志 ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)]
)
log = logging.getLogger(__name__)

# --- 1. 定义文件路径和列名常量 (将从 Snakemake 对象获取) ---
# ANNOTATION_FILE = snakemake.input.annotation_file
# VARIANT_MATRIX_FILE = snakemake.input.merged_inversion # 假设 snakemake input 名
# OUTPUT_ANNOTATED_VARIANTS_FILE = snakemake.output.anno_inversion # 假设 snakemake output 名

ALL_ANNOT_COLUMNS = [
    'genome_id', 'genome_name', 'accession', 'feature_type', 'patric_id',
    'start', 'end', 'strand', 'gene', 'product', 'global_start', 'global_end'
]
FINAL_OUTPUT_ORDER = ['Site'] + ALL_ANNOT_COLUMNS

COL_CHROMOSOME = "Chromosome"
COL_START_0B = "Start"
COL_END_0B = "End"
COL_STRAND = "Strand"
COL_START_1B_FROM_ANNOT = "Start_1based_from_annot"
COL_END_1B_FROM_ANNOT = "End_1based_from_annot"

def annotate_inversions(annotation_fp, inversion_matrix_fp, output_fp):
    """
    Annotates DNA inversions from a matrix file using a feature annotation file.
    """
    log.info(f"Starting inversion annotation.")
    log.info(f"Annotation file: {annotation_fp}")
    log.info(f"Inversion matrix file: {inversion_matrix_fp}")
    log.info(f"Output file: {output_fp}")

    # --- 2. 读取注释文件 ---
    log.info(f"Reading annotation file: {annotation_fp}")
    try:
        annot_df_dtypes = {col: str for col in ['genome_id', 'accession', 'patric_id', 'gene', 'strand', 'feature_type', 'product']}
        annot_df_original = pd.read_csv(annotation_fp, sep='\t', low_memory=False, dtype=annot_df_dtypes)
        log.info(f"Annotation file loaded. Shape: {annot_df_original.shape}")
    except FileNotFoundError:
        log.error(f"Annotation file not found at {annotation_fp}")
        raise
    except Exception as e:
        log.error(f"Error reading annotation file {annotation_fp}: {e}")
        raise

    missing_annot_cols = [col for col in ALL_ANNOT_COLUMNS if col not in annot_df_original.columns]
    if missing_annot_cols:
        log.error(f"Annotation file '{annotation_fp}' is missing required columns: {missing_annot_cols}")
        raise ValueError(f"Annotation file missing columns: {missing_annot_cols}")

    # --- 2a. 创建 accession 到 genome_id 和 genome_name 的映射 ---
    try:
        accession_to_genome_info_map = annot_df_original.drop_duplicates(subset=['accession']) \
                                                        [['accession', 'genome_id', 'genome_name']] \
                                                        .set_index('accession').to_dict('index')
        log.info("Accession to genome_id/genome_name map created.")
    except KeyError as e:
        log.error(f"Error creating accession map: Missing column '{e}'. Needed: 'accession', 'genome_id', 'genome_name'.")
        raise

    # --- 2b. 准备注释文件用于 PyRanges ---
    annot_df_for_pyranges = annot_df_original.copy()
    annot_df_for_pyranges.rename(columns={'accession': COL_CHROMOSOME,
                                          'start': COL_START_1B_FROM_ANNOT,
                                          'end': COL_END_1B_FROM_ANNOT},
                                 inplace=True)
    if 'strand' in annot_df_for_pyranges.columns and COL_STRAND not in annot_df_for_pyranges.columns:
        annot_df_for_pyranges.rename(columns={'strand': COL_STRAND}, inplace=True)
    elif COL_STRAND not in annot_df_for_pyranges.columns and 'strand' not in annot_df_for_pyranges.columns:
        log.warning(f"Neither '{COL_STRAND}' nor 'strand' column found in annotation. Intervals will be unstranded.")
        annot_df_for_pyranges[COL_STRAND] = "."

    try:
        annot_df_for_pyranges[COL_START_0B] = annot_df_for_pyranges[COL_START_1B_FROM_ANNOT].astype(int) - 1
        annot_df_for_pyranges[COL_END_0B] = annot_df_for_pyranges[COL_END_1B_FROM_ANNOT].astype(int)
    except KeyError as e:
        log.error(f"Missing renamed 1-based start/end columns for 0-based conversion: {e}.")
        raise
    except ValueError as e:
        log.error(f"Could not convert 1-based start/end columns to integer: {e}")
        raise

    try:
        gr_annotations = pr.PyRanges(annot_df_for_pyranges)
        log.info("Annotation data converted to PyRanges object for interval joining.")
    except Exception as e:
        log.error(f"Error creating PyRanges object from annotation_df: {e}")
        raise

    # --- 3. 读取并处理 Inversion 文件 ---
    log.info(f"Reading inversion matrix file: {inversion_matrix_fp}")
    try:
        inv_matrix_df = pd.read_csv(inversion_matrix_fp, sep='\t', dtype={'Site': str})
        log.info(f"Inversion matrix file loaded. Shape: {inv_matrix_df.shape}")
    except FileNotFoundError:
        log.error(f"Inversion matrix file not found at {inversion_matrix_fp}")
        raise
    except Exception as e:
        log.error(f"Error reading inversion matrix file {inversion_matrix_fp}: {e}")
        raise

    if 'Site' not in inv_matrix_df.columns:
        log.error(f"'Site' column not found in inversion matrix file: {inversion_matrix_fp}")
        raise ValueError("'Site' column missing in inversion matrix file.")

    parsed_inversions = []
    for site_str in inv_matrix_df['Site']:
        try:
            parts = site_str.split(':')
            chrom = parts[0]
            coords_str = parts[1]
            coords = [int(c) for c in coords_str.split('-')]
            if len(coords) != 4:
                log.warning(f"Site string '{site_str}' does not have 4 coordinate parts. Skipping.")
                continue
            
            # Inversion region is between the two inverted repeats
            # Your definition: END1 (coords[1]) and START2 (coords[2]) are the boundaries (1-based inclusive)
            # For PyRanges (0-based, half-open [Start, End) ):
            # Start = END1 - 1
            # End = START2
            # However, if you mean the region *between* END1 and START2, then:
            # Start = END1 (which is coords[1], 0-based would be coords[1])
            # End = START2 -1 (which is coords[2]-1, 0-based would be coords[2]-1)
            # Let's stick to your definition: "起始终止位点选取：182954和183130"
            # This means for site ADLE01000001:182930-182954-183130-183154,
            # inversion_start_1based = 182954 (coords[1])
            # inversion_end_1based   = 183130 (coords[2])
            
            inversion_start_1based = coords[1] # This is END1
            inversion_end_1based = coords[2]   # This is START2

            # Ensure start <= end for the inversion region itself
            if inversion_start_1based > inversion_end_1based:
                log.warning(f"Parsed inversion start ({inversion_start_1based}) is greater than end ({inversion_end_1based}) for site '{site_str}'. Skipping.")
                continue

            parsed_inversions.append({
                COL_CHROMOSOME: chrom,
                COL_START_0B: inversion_start_1based - 1, # Convert 1-based inclusive start to 0-based
                COL_END_0B: inversion_end_1based,       # Convert 1-based inclusive end to 0-based exclusive
                'Site': site_str,
                'Site_Chromosome_parsed': chrom
            })
        except (ValueError, IndexError) as e:
            log.warning(f"Could not parse site string '{site_str}': {e}. Skipping.")
            continue

    if not parsed_inversions:
        log.warning("No inversions could be parsed from the matrix file. Output might be empty.")
        variants_for_pr_df = pd.DataFrame(columns=[COL_CHROMOSOME, COL_START_0B, COL_END_0B, 'Site', 'Site_Chromosome_parsed'])
    else:
        variants_for_pr_df = pd.DataFrame(parsed_inversions)


    try:
        gr_variants = pr.PyRanges(variants_for_pr_df) # 'variants' here refers to inversions
        log.info("Inversion data converted to PyRanges object.")
    except Exception as e:
        log.error(f"Error creating PyRanges object from inversions_df: {e}")
        raise

    # --- 4. 使用 pyranges.join() 进行注释 ---
    log.info("Performing interval-based join to annotate inversions...")
    annotated_gr = gr_variants.join(gr_annotations, how='left', suffix='_annot', apply_strand_suffix=False)
    annotated_result_df = annotated_gr.df
    log.info(f"Interval join complete. Result shape: {annotated_result_df.shape}")

    # --- 5. 整理输出结果 ---
    log.info("Formatting final output...")
    final_output_list = []

    for _, row in annotated_result_df.iterrows():
        output_row = {'Site': row['Site']}
        site_chromosome_parsed = row['Site_Chromosome_parsed']

        genome_info = accession_to_genome_info_map.get(site_chromosome_parsed)
        if genome_info:
            output_row['genome_id'] = genome_info.get('genome_id', pd.NA)
            output_row['genome_name'] = genome_info.get('genome_name', pd.NA)
        else:
            output_row['genome_id'] = pd.NA
            output_row['genome_name'] = pd.NA

        feature_matched = ('patric_id' in row and pd.notna(row['patric_id'])) or \
                          (COL_START_0B + "_annot" in row and pd.notna(row[COL_START_0B + "_annot"]))

        if feature_matched:
            for col_name in ALL_ANNOT_COLUMNS:
                if col_name in ['genome_id', 'genome_name']:
                    continue
                if col_name == 'accession':
                    output_row[col_name] = row.get(COL_CHROMOSOME + "_annot", row.get(COL_CHROMOSOME, pd.NA)) # Prefer annot's chrom if available
                elif col_name == 'start':
                    output_row[col_name] = row.get(COL_START_1B_FROM_ANNOT, pd.NA)
                elif col_name == 'end':
                    output_row[col_name] = row.get(COL_END_1B_FROM_ANNOT, pd.NA)
                elif col_name == 'strand':
                    output_row[col_name] = row.get(COL_STRAND, pd.NA) # From previous debugging, join result has 'Strand'
                else:
                    output_row[col_name] = row.get(col_name, pd.NA)
        else:
            output_row['accession'] = site_chromosome_parsed
            for col_name in ALL_ANNOT_COLUMNS:
                if col_name not in output_row:
                    output_row[col_name] = pd.NA

        final_output_list.append(output_row)

    final_output_df = pd.DataFrame(final_output_list)

    for col in FINAL_OUTPUT_ORDER:
        if col not in final_output_df.columns:
            log.warning(f"Final output column '{col}' was not generated. Adding it as NA.")
            final_output_df[col] = pd.NA
    final_output_df = final_output_df[FINAL_OUTPUT_ORDER]

    # --- 6. 保存结果 ---
    log.info(f"Saving annotated inversions to: {output_fp}")
    try:
        final_output_df.to_csv(output_fp, sep='\t', index=False, header=True)
        log.info("Annotation complete. Output file saved successfully.")
    except Exception as e:
        log.error(f"Error saving final file {output_fp}: {e}")
        raise

if __name__ == '__main__':
    try:
        annotation_input_file = snakemake.input.annotation_file
        # Assuming your Snakemake rule input for inversions is named 'merged_inversion'
        variant_matrix_input_file = snakemake.input.merged_inversion
        output_file = snakemake.output.anno_inversion

        annotate_inversions(annotation_input_file, variant_matrix_input_file, output_file)

    except NameError:
        log.info("Running script in standalone mode (not via Snakemake). Using hardcoded file paths for testing.")
        # Define hardcoded paths for testing if not run by Snakemake
        # These should match the constants at the top if you want to test with them
        TEST_ANNOTATION_FILE = '/home/data/CYM/pipeline/database/43_species_features.tsv'
        TEST_INVERSION_MATRIX_FILE = '/home/data/CYM/pipeline/00-20-16core/02-variant-calling/SV/PhaseFinder/across-sample/inversion.tsv' # Example path
        TEST_OUTPUT_FILE = '/home/data/CYM/pipeline/00-20-16core/02-variant-calling/SV/PhaseFinder/across-sample/inversion_anno.tsv'
        annotate_inversions(TEST_ANNOTATION_FILE, TEST_INVERSION_MATRIX_FILE, TEST_OUTPUT_FILE)