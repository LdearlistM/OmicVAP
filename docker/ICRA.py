import logging
import os
import sys
from snakemake.script import snakemake
from SGVFinder2 import single_file, get_sample_map
# Make sure pandas is available if to_pickle comes from it,
# or adjust if it's a different library's to_pickle
try:
    from pandas import to_pickle
except ImportError:
    import pickle
    # Define a simple to_pickle if pandas is not used/available
    def to_pickle(obj, filepath):
        with open(filepath, 'wb') as f:
            pickle.dump(obj, f)

# --- Configure Logging ---
log_file = snakemake.log[0]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,  # Or logging.DEBUG
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
# Optional: Log to console too
# logging.getLogger().addHandler(logging.StreamHandler(sys.stderr))

logging.info("Starting SGVFinder2 ICRA script.")

try:
    # --- Get Inputs/Parameters ---
    # read1 = snakemake.input.R1 # Use named input
    # read2 = snakemake.input.R2 # Use named input
    read1 = snakemake.params.fq1
    read2_param = snakemake.params.fq2
    read2 = None if read2_param.lower() == 'none' else read2_param
    db_index_main_file_in_container = snakemake.input.db_index_file
    DATABASE = db_index_main_file_in_container.replace(".1.bt2", "")
    
    # DATABASE = snakemake.params.database
    output_dir = snakemake.params.output_dir
    smp_output_file = snakemake.output.smp_file # Use named output
    
    

    

    logging.info(f"Input R1: {read1}")
    logging.info(f"Input R2: {read2}")
    logging.info(f"Database path prefix: {DATABASE}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Target SMP output file: {smp_output_file}")

    # Ensure output directory exists (Snakemake usually handles this, but belt-and-suspenders)
    os.makedirs(output_dir, exist_ok=True)

    # --- Core Logic ---
    logging.info("Running SGVFinder2 single_file...")
    # Note: Output printed by single_file itself will go to stdout/stderr
    # and be captured by Snakemake into the log file.
    jspi_file, jsdel_file = single_file(
        fq1=read1,
        fq2=read2,
        outfol=output_dir,
        dbpath=DATABASE
    )
    logging.info(f"single_file completed. Intermediate outputs: jspi={jspi_file}, jsdel={jsdel_file}")

    logging.info("Running SGVFinder2 get_sample_map...")
    # Construct the expected database length file path
    db_dlen_file = DATABASE + '.dlen'
    if not os.path.exists(db_dlen_file):
         logging.warning(f"Database length file may not exist: {db_dlen_file}")
    # Note: Output printed by get_sample_map itself will go to stdout/stderr
    sample_map = get_sample_map(jsdel_file, db_dlen_file)
    logging.info("get_sample_map completed.")

    logging.info(f"Saving sample map to pickle file: {smp_output_file}")
    to_pickle(sample_map, smp_output_file) # Save to the precise output path
    logging.info("Sample map successfully saved.")

    logging.info("ICRA script finished successfully.")

except Exception as e:
    logging.exception("An error occurred during ICRA script execution.") # Log exception with traceback
    sys.exit(1) # Exit with error status to indicate failure to Snakemake