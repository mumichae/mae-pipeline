import sys
# Add the folder path for the python parsing functions to the sys.path list
sys.path.insert(0,'../genetic_diagnosis_modified/src/python') 
from config_helper import ConfigHelper

configfile: "wbuild.yaml" 
parser = ConfigHelper(config)

subworkflow variantsPipeline:
    workdir:
        "../variant-annotation-pipeline"
    snakefile:
        "../variant-annotation-pipeline/Snakefile"
    configfile:
        "../variant-annotation-pipeline/wbuild.yaml"
        
        
## Needed for MAE: set config variables for mae
vcfs, rnas = parser.getMaeIDs()
config["vcfs"] = vcfs
config["rnas"] = rnas
config["mae_ids"] = list(map('--'.join, zip(vcfs, rnas)))
print(config["mae_ids"])


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
    input: expand(parser.getProcDataDir() + "/mae/{id}.Rds", id=config["mae_ids"])
