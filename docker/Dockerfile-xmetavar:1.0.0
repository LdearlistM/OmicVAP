FROM snakemake/snakemake:v8.25.3

WORKDIR /pipeline

RUN conda install -n base -c conda-forge -c bioconda pyranges -y && \
    conda clean --all --yes

COPY Snakefile .
COPY scripts/ ./scripts/
COPY envs/ ./envs/
COPY software/ ./software/
COPY conda/ /home/data/CYM/pipeline/.snakemake/conda

COPY run_pipeline.sh /run_pipeline.sh
RUN chmod +x /run_pipeline.sh
ENTRYPOINT ["/run_pipeline.sh"]
CMD ["--help"]
