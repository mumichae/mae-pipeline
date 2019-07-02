import sys
# Add the folder path for the python parsing functions to the sys.path list
from config_parser import ConfigHelper

parser = ConfigHelper(config)

 
include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
#htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
htmlOutputPath = "Output/html"


rule test_allelic_counts:
    input: parser.getProcDataDir() + "/mae/MLL_11513--MLL_11513-M043.Rds"

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")
    
# overwriting wbuild rule output
rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

rule allelic_counts:
    input: expand(parser.getProcDataDir() + "/mae/{id}.Rds", id=parser.getMaeIDs())
