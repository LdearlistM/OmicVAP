from snakemake.script import snakemake
from SGVFinder2 import work_on_collection
import pandas as pd

samp_to_map = snakemake.params.input_dir
output_dir = snakemake.params.output_dir
DATABASE = snakemake.params.database


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
    # browser_path="/home/data/CYM/pipeline/04-variant-calling/SV/SGVFinder/browser",
    #browser_path=None, 
    taxonomypath=DATABASE+'.taxonomy.df',
    genepospath=DATABASE+'genepos.df',
    frames_path="/home/data/CYM/pipeline/02-variant-calling/SV/SGVFinder2/merge"
)
vsgv.to_csv(f'{output_dir}/vsgv.csv')
dsgv.to_csv(f'{output_dir}/dsgv.csv')
# vsgv.to_pickle(f'{output_dir}/vsgv.df')
# dsgv.to_pickle(f'{output_dir}/dsgv.df')