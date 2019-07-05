### SNAKEFILE MONOALLELIC EXPRESSION

import os
from config_parser import ConfigHelper

parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
htmlOutputPath = config["htmlOutputPath"]

include: os.getcwd() + "/.wBuild/wBuild.snakefile"
# create temporary folder
if not os.path.exists('tmp'):
    os.makedirs('tmp')

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")
    
# overwriting wbuild rule output
rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

rule allelic_counts:
    input: expand(parser.getProcDataDir() + "/mae/{id}.Rds", id=parser.getMaeIDs())
