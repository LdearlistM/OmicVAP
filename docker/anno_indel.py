import pandas as pd
import pyranges as pr
import sys # 用于从 Snakemake 获取参数
import logging

# --- 0. 设置日志 ---
# Snakemake 会将日志重定向到指定的文件
# 如果直接运行脚本，日志会输出到标准错误
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)] # 或者 sys.stdout 如果你更喜欢
)
log = logging.getLogger(__name__)


# --- 1. 定义文件路径和列名常量 (将从 Snakemake 对象获取) ---
# 这些将在脚本的 main 部分从 snakemake.input 和 snakemake.output 获取
# ANNOTATION_FILE = snakemake.input.annotation_file
# VARIANT_MATRIX_FILE = snakemake.input.merged_indel
# OUTPUT_ANNOTATED_VARIANTS_FILE = snakemake.output.anno_indel

# 注释文件中的所有期望列（也是最终输出中来自注释的部分）
# 这些列名应与 ANNOTATION_FILE 中的列名完全一致
ALL_ANNOT_COLUMNS = [
    'genome_id', 'genome_name', 'accession', 'feature_type', 'patric_id',
    'start', 'end', 'strand', 'gene', 'product', 'global_start', 'global_end'
]
# 最终输出文件的列顺序
FINAL_OUTPUT_ORDER = ['Site'] + ALL_ANNOT_COLUMNS

# PyRanges 使用的列名常量
COL_CHROMOSOME = "Chromosome"
COL_START_0B = "Start" # 0-based start for PyRanges
COL_END_0B = "End"     # 0-based end (exclusive) for PyRanges
COL_STRAND = "Strand"

# 注释文件中原始的1-based坐标列名 (在准备 PyRanges DataFrame 时创建)
COL_START_1B_FROM_ANNOT = "Start_1based_from_annot" # 标识来自注释的原始1-based start
COL_END_1B_FROM_ANNOT = "End_1based_from_annot"   # 标识来自注释的原始1-based end


def annotate_variants(annotation_fp, variant_matrix_fp, output_fp):
    """
    Annotates variants from a matrix file using a feature annotation file.

    Args:
        annotation_fp (str): Path to the feature annotation TSV file.
        variant_matrix_fp (str): Path to the variant matrix TSV file.
        output_fp (str): Path to save the annotated variants TSV file.
    """
    log.info(f"Starting variant annotation.")
    log.info(f"Annotation file: {annotation_fp}")
    log.info(f"Variant matrix file: {variant_matrix_fp}")
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

    # 检查注释文件是否包含所有期望的列
    missing_annot_cols = [col for col in ALL_ANNOT_COLUMNS if col not in annot_df_original.columns]
    if missing_annot_cols:
        log.error(f"Annotation file '{annotation_fp}' is missing required columns: {missing_annot_cols}")
        # Snakemake 会捕获异常并失败
        raise ValueError(f"Annotation file missing columns: {missing_annot_cols}")


    # --- 2a. 创建 accession 到 genome_id 和 genome_name 的映射 ---
    try:
        accession_to_genome_info_map = annot_df_original.drop_duplicates(subset=['accession']) \
                                                        [['accession', 'genome_id', 'genome_name']] \
                                                        .set_index('accession').to_dict('index')
        log.info("Accession to genome_id/genome_name map created.")
    except KeyError as e:
        log.error(f"Error creating accession map: Missing column '{e}' in annotation file. Needed: 'accession', 'genome_id', 'genome_name'.")
        raise

    # --- 2b. 准备注释文件用于 PyRanges ---
    annot_df_for_pyranges = annot_df_original.copy()
    annot_df_for_pyranges.rename(columns={'accession': COL_CHROMOSOME,
                                          'start': COL_START_1B_FROM_ANNOT, # 使用新常量名
                                          'end': COL_END_1B_FROM_ANNOT},   # 使用新常量名
                                 inplace=True)

    if 'strand' in annot_df_for_pyranges.columns and COL_STRAND not in annot_df_for_pyranges.columns:
        annot_df_for_pyranges.rename(columns={'strand': COL_STRAND}, inplace=True)
    elif COL_STRAND not in annot_df_for_pyranges.columns and 'strand' not in annot_df_for_pyranges.columns:
        log.warning(f"Neither '{COL_STRAND}' nor 'strand' column found in annotation. Intervals will be unstranded.")
        annot_df_for_pyranges[COL_STRAND] = "." # 添加一个默认的无链列

    try:
        annot_df_for_pyranges[COL_START_0B] = annot_df_for_pyranges[COL_START_1B_FROM_ANNOT].astype(int) - 1
        annot_df_for_pyranges[COL_END_0B] = annot_df_for_pyranges[COL_END_1B_FROM_ANNOT].astype(int)
    except KeyError as e:
        log.error(f"Missing renamed 1-based start/end columns for 0-based conversion: {e}.")
        raise
    except ValueError as e:
        log.error(f"Could not convert 1-based start/end columns to integer for 0-based conversion: {e}")
        raise

    try:
        gr_annotations = pr.PyRanges(annot_df_for_pyranges)
        log.info("Annotation data converted to PyRanges object for interval joining.")
    except Exception as e:
        log.error(f"Error creating PyRanges object from annotation_df: {e}")
        raise

    # --- 3. 读取并处理变异文件 ---
    log.info(f"Reading variant matrix file: {variant_matrix_fp}")
    try:
        var_matrix_df = pd.read_csv(variant_matrix_fp, sep='\t', dtype={'Site': str})
        log.info(f"Variant matrix file loaded. Shape: {var_matrix_df.shape}")
    except FileNotFoundError:
        log.error(f"Variant matrix file not found at {variant_matrix_fp}")
        raise
    except Exception as e:
        log.error(f"Error reading variant matrix file {variant_matrix_fp}: {e}")
        raise

    if 'Site' not in var_matrix_df.columns:
        log.error(f"'Site' column not found in variant matrix file: {variant_matrix_fp}")
        raise ValueError("'Site' column missing in variant matrix file.")

    parsed_variants = []
    for site_str in var_matrix_df['Site']:
        try:
            # chrom, pos_str = site_str.split(':', 1)
            # pos = int(pos_str)
            parts = site_str.split(':')
            if len(parts) == 3:
                chrom, pos_str, indel_type_str = parts[0], parts[1], parts[2]
                pos = int(pos_str)
                # indel_type = int(indel_type_str) # Type is already -1 or 1 (or 0) as string from previous script
                                                 # We don't strictly need to convert it to int here unless used for logic
            elif len(parts) == 2: # Fallback or handling for old CHR:POS format if mixed data is possible
                log.warning(f"Site string '{site_str}' appears to be in old CHR:POS format. Assuming no type info.")
                chrom, pos_str = parts[0], parts[1]
                pos = int(pos_str)
                # indel_type_str = "0" # Or some other default if type is missing
            else:
                raise ValueError(f"Site string '{site_str}' is not in expected CHR:POS:Type or CHR:POS format.")
            
            parsed_variants.append({
                COL_CHROMOSOME: chrom,
                COL_START_0B: pos - 1,
                COL_END_0B: pos,
                'Site': site_str,
                'Site_Chromosome_parsed': chrom
            })
        except ValueError:
            log.warning(f"Could not parse site string: '{site_str}'. Skipping.")
            continue

    if not parsed_variants:
        log.error("No variants could be parsed from the variant matrix file. Output will be empty or incomplete.")
        # 创建一个空的 DataFrame 以便后续步骤能运行，但结果会是空的
        variants_for_pr_df = pd.DataFrame(columns=[COL_CHROMOSOME, COL_START_0B, COL_END_0B, 'Site', 'Site_Chromosome_parsed'])
    else:
        variants_for_pr_df = pd.DataFrame(parsed_variants)


    try:
        gr_variants = pr.PyRanges(variants_for_pr_df)
        log.info("Variant data converted to PyRanges object.")
    except Exception as e:
        log.error(f"Error creating PyRanges object from variants_df: {e}")
        raise

    # --- 4. 使用 pyranges.join() 进行基于区间的特征注释 ---
    log.info("Performing interval-based join to annotate features...")
    # apply_strand_suffix=False 避免 PyRanges 自动添加 _strand 后缀，
    # 因为我们希望从右表（gr_annotations）获取 Strand 信息。
    # 同时这也抑制了之前看到的提示信息。
    annotated_gr = gr_variants.join(gr_annotations, how='left', suffix='_annot', apply_strand_suffix=False)
    annotated_result_df = annotated_gr.df
    log.info(f"Interval join complete. Result shape: {annotated_result_df.shape}")

    # --- 5. 整理输出结果 ---
    log.info("Formatting final output...")
    final_output_list = []
    joined_cols = annotated_result_df.columns.tolist() # 获取 join 后实际的列名

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

        # 判断是否匹配到特征。'patric_id' 是一个只在注释中存在的明确特征。
        # 或者检查 PyRanges join 后通常会为右表 Start 列添加的后缀。
        feature_matched = ('patric_id' in row and pd.notna(row['patric_id'])) or \
                          (COL_START_0B + "_annot" in row and pd.notna(row[COL_START_0B + "_annot"]))

        if feature_matched:
            for col_name in ALL_ANNOT_COLUMNS:
                if col_name in ['genome_id', 'genome_name']:
                    continue # 已经从 map 回填

                # 获取来自注释特征的值
                # PyRanges核心列在join后会带后缀 '_annot' (来自右表 gr_annotations)
                # 非核心列（如 patric_id, gene, product 等）会直接保留原名，因为 gr_variants 中没有这些列
                if col_name == 'accession':
                    # 注释特征的染色体/accession。gr_annotations 的 Chromosome 列会变成 Chromosome_annot
                    output_row[col_name] = row.get(COL_CHROMOSOME, pd.NA)
                elif col_name == 'start': # 我们要的是原始的 1-based start from annotation
                    output_row[col_name] = row.get(COL_START_1B_FROM_ANNOT, pd.NA)
                elif col_name == 'end':   # 我们要的是原始的 1-based end from annotation
                    output_row[col_name] = row.get(COL_END_1B_FROM_ANNOT, pd.NA)
                elif col_name == 'strand': # 注释特征的 strand
                    # gr_annotations 的 Strand 列在 join 后会是 Strand_annot
                    # (因为gr_variants默认没有Strand,所以右表的Strand列会被保留并加后缀)
                    # 如果apply_strand_suffix=False, 右表的Strand会直接成为结果中的Strand列
                    # 从之前的调试看，它直接是 'Strand'，所以这里我们尝试 'Strand'
                    output_row[col_name] = row.get(COL_STRAND, pd.NA)
                else:
                    output_row[col_name] = row.get(col_name, pd.NA)
        else: # 未匹配到具体特征区间 (基因间区)
            output_row['accession'] = site_chromosome_parsed
            for col_name in ALL_ANNOT_COLUMNS:
                if col_name not in output_row:
                    output_row[col_name] = pd.NA

        final_output_list.append(output_row)

    final_output_df = pd.DataFrame(final_output_list)

    # 重新排序并确保所有期望的列都存在
    for col in FINAL_OUTPUT_ORDER:
        if col not in final_output_df.columns:
            log.warning(f"Final output column '{col}' was not generated. Adding it as NA.")
            final_output_df[col] = pd.NA
    final_output_df = final_output_df[FINAL_OUTPUT_ORDER]

    # --- 6. 保存结果 ---
    log.info(f"Saving annotated variants to: {output_fp}")
    try:
        final_output_df.to_csv(output_fp, sep='\t', index=False, header=True)
        log.info("Annotation complete. Output file saved successfully.")
    except Exception as e:
        log.error(f"Error saving final file {output_fp}: {e}")
        raise

if __name__ == '__main__':
    # This part is for Snakemake integration
    # It expects snakemake object to be defined when run by Snakemake
    try:
        # Input files from Snakemake
        annotation_input_file = snakemake.input.annotation_file
        variant_matrix_input_file = snakemake.input.merged_indel
        # Output file from Snakemake
        output_file = snakemake.output.anno_indel
        # Log file from Snakemake (Python logging will go here via stderr redirection)
        # log_file = snakemake.log[0] # Not directly used by basicConfig, but good for context

        # Reconfigure logging to also go to Snakemake log file if desired,
        # but Snakemake's default stderr redirection to log file is usually sufficient.
        # logging.basicConfig(filename=log_file, level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

        annotate_variants(annotation_input_file, variant_matrix_input_file, output_file)

    except NameError:
        # (when 'snakemake' object is not defined)
        log.info("Running script in standalone mode (not via Snakemake). Using hardcoded file paths.")
        # 使用脚本顶部的常量定义的路径进行测试
        annotate_variants(ANNOTATION_FILE, VARIANT_MATRIX_FILE, OUTPUT_ANNOTATED_VARIANTS_FILE)