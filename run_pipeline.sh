#!/bin/bash
set -e

declare -A TARGET_MAP
RESULTS_DIR="/pipeline/results"
VC_DIR="${RESULTS_DIR}/02-variant-calling"
TARGET_MAP["snp_gtpro"]="${VC_DIR}/SNP/GT-Pro/across-sample/snp.tsv"
TARGET_MAP["snp_midas"]="${RESULTS_DIR}/logs/MIDASv3/snp_done.txt ${VC_DIR}/SNP/MIDAS/across_sample/snps/snps_summary.tsv"
TARGET_MAP["indel"]="${VC_DIR}/INDEL/across-sample/indel_anno.tsv"
TARGET_MAP["sv_sgvfinder"]="${VC_DIR}/SV/dSV-vSV/across-sample/dsgv_anno.tsv ${VC_DIR}/SV/dSV-vSV/across-sample/vsgv_anno.tsv"
TARGET_MAP["sv_inversion"]="${VC_DIR}/SV/Inversion/across-sample/inversion_anno.tsv"
TARGET_MAP["sv_midas"]="${RESULTS_DIR}/logs/MIDASv3/cnv_done.txt ${VC_DIR}/SV/CNV/across_sample/merge.genes_copynum.tsv"

SNAKEMAKE_TARGETS=()
SNAKEMAKE_EXTRA_ARGS=()

for arg in "$@"; do
    if [[ -v TARGET_MAP["$arg"] ]]; then
        SNAKEMAKE_TARGETS+=(${TARGET_MAP[$arg]})
        echo "Recognized target key: '$arg'. Adding target(s): ${TARGET_MAP[$arg]}"
    elif [[ "$arg" == "all" ]]; then
        SNAKEMAKE_TARGETS=()
        echo "Recognized 'all' target. The entire workflow will be executed."
    else
        SNAKEMAKE_EXTRA_ARGS+=("$arg")
    fi
done

if [ ${#SNAKEMAKE_TARGETS[@]} -eq 0 ] && [ ${#SNAKEMAKE_EXTRA_ARGS[@]} -eq 0 ] && [ "$#" -gt 0 ] && [ "$1" != "all" ]; then
    SNAKEMAKE_TARGETS=()
    SNAKEMAKE_EXTRA_ARGS=("$@")
    echo "No predefined target keys matched. Passing all arguments directly to Snakemake."
elif [ "$#" -eq 0 ]; then
    echo "No arguments provided. Running the complete workflow (default)."
fi

echo "--- Starting Snakemake ---"
echo "Final targets: ${SNAKEMAKE_TARGETS[@]}"
echo "Extra arguments: ${SNAKEMAKE_EXTRA_ARGS[@]}"

snakemake \
    --cores 16 \
    --configfile /pipeline/config.yaml \
    --use-conda \
    --conda-prefix /home/data/CYM/pipeline/.snakemake/conda \
    "${SNAKEMAKE_TARGETS[@]}" \
    "${SNAKEMAKE_EXTRA_ARGS[@]}"

SNAKEMAKE_EXIT_CODE=$?

echo "--- Snakemake finished with exit code ${SNAKEMAKE_EXIT_CODE}. Setting permissions... ---"
chmod -R 777 /pipeline/results
echo "--- Permissions set. Workflow complete. ---"

exit ${SNAKEMAKE_EXIT_CODE}
