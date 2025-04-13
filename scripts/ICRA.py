from snakemake.script import snakemake
from SGVFinder2 import single_file, get_sample_map
from pandas import to_pickle

read1 = snakemake.input[0]
read2 = snakemake.input[1]

# DATABASE = 'my_db_folder/my_custom_db'
DATABASE = snakemake.params.database

output_dir = snakemake.params.output_dir

jspi_file, jsdel_file = single_file(
    fq1=read1,
    fq2=read2,
    outfol=output_dir, 
    dbpath=DATABASE
)

sample_map = get_sample_map(jsdel_file, DATABASE + '.dlen') 

to_pickle(sample_map, jsdel_file.replace('.jsdel', '.smp'))