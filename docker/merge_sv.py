# scripts/merge_sv.py
import os
import shutil
import sys
from snakemake.script import snakemake
from SGVFinder2 import work_on_collection
import pandas as pd
import logging # Import the logging library
import warnings  # 1. 导入 warnings 模块


# --- Configure Logging ---
# Get log file path from Snakemake
log_file = snakemake.log[0]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,  # Or logging.DEBUG for more detail
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
# Optional: To also see logs on console during development/debugging, uncomment next line
# logging.getLogger().addHandler(logging.StreamHandler(sys.stderr))

# --- Script Start ---
logging.info("Starting SGVFinder2 merge script.")
logging.info(f"Input directory: {snakemake.params.input_dir}")
logging.info(f"Output directory: {snakemake.params.output_dir}")
logging.info(f"Frames path: {snakemake.params.frames_path}")
logging.info(f"Log file: {log_file}")

# --- Main Processing Logic ---
try: # Wrap everything in a try...except for robustness
    samp_to_map = snakemake.params.input_dir
    output_dir = snakemake.params.output_dir
    frames_path = snakemake.params.frames_path

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(frames_path, exist_ok=True)

    logging.info(f"Running work_on_collection on {samp_to_map}")
    
    with warnings.catch_warnings():
        # 3. 在这里设置过滤器，忽略特定类型的警告
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)
        # warnings.filterwarnings("ignore", message="Mean of empty slice") # 也可以按消息内容过滤
        vsgv, dsgv = work_on_collection(
            samp_to_map=samp_to_map,
            max_spacing=10,
            min_samp_cutoff=2,
            delsdetectthresh=0.25,
            real_del_thresh=0.95,
            dels_cooc_thresh=0.25,
            vsgv_dissim_thresh=0.125,
            vsgv_clip_quantile=0.02,
            vsgv_fit_interval=0.95,
            vsgv_fit_method='betaprime',
            x_coverage=0.01,
            rate_param=10,
            vsgv_dense_perc=85,
            browser_path=None,
            frames_path=frames_path
        )
    logging.info("work_on_collection completed.")

    if pd.api.types.is_object_dtype(vsgv.index.dtype):
        logging.info("Removing _paired suffix from vsgv index.")
        vsgv.index = vsgv.index.str.removesuffix('_paired')

    if pd.api.types.is_object_dtype(dsgv.index.dtype):
        logging.info("Removing _paired suffix from dsgv index.")
        dsgv.index = dsgv.index.str.removesuffix('_paired')

    logging.info(f"Writing vsgv.csv to {output_dir}")
    vsgv.to_csv(os.path.join(output_dir, 'vsgv.csv'), index=True) # Keep index
    logging.info(f"Writing dsgv.csv to {output_dir}")
    dsgv.to_csv(os.path.join(output_dir, 'dsgv.csv'), index=True) # Keep index

    if os.path.isdir(frames_path):
        logging.info(f"Processing .df files in: {frames_path}")
        processed_count = 0
        error_count = 0
        for filename in os.listdir(frames_path):
            if filename.endswith(".df") and os.path.isfile(os.path.join(frames_path, filename)):
                df_filepath = os.path.join(frames_path, filename)
                csv_filename = filename[:-3] + ".csv"
                csv_filepath = os.path.join(frames_path, csv_filename)
                logging.info(f"  Attempting to convert {filename} to {csv_filename}...")

                try:
                    df_content = pd.read_pickle(df_filepath)
                    logging.info(f"    Successfully read {filename}.")

                    df_content.to_csv(csv_filepath, index=True)  # FIXED: index=True
                    logging.info(f"    Successfully wrote {csv_filename}.")

                    os.remove(df_filepath)
                    logging.info(f"    Successfully removed {filename}.")
                    processed_count += 1

                except Exception as e:
                    error_count += 1
                    logging.warning(f"Warning: Error processing file {filename}: {e}")
                    # Retain .df file for debugging
        logging.info(f"Finished processing .df files. Converted: {processed_count}, Errors: {error_count}")

    # if os.path.isdir(samp_to_map):
    #     logging.info(f"Removing temporary directory: {samp_to_map}")
    #     try:
    #         shutil.rmtree(samp_to_map)
    #         logging.info(f"Successfully removed {samp_to_map}")
    #     except OSError as e:
    #         logging.warning(f"Warning: Could not delete {samp_to_map}: {e}")
except Exception as e:
    logging.exception("An unexpected error occurred during script execution.")
    sys.exit(1)