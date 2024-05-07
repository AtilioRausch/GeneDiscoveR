"Author: Atilio O. Rausch"
"Date: 2024-04-06"
"Email: atiliorau@gmail.com"

import numpy as np

"""This script is used to run Orthofinder in parallel for different inflation values.
    The inflation values are defined by the user in the config.yaml file.
    The parameters used for running Orthofinder are:
        - Minimum inflation value: {config["imin"]}
        - Maximum inflation value: {config["imax"]}
        - Inflation step: {config["istep"]}
        - Proteomes: {config["proteomes"]}
        - Threads per core: {config["threads"]}
"""
rule Orthofinder:
    input:
        proteomes = config["proteomes"]
    params:
        imin = config["imin"],
        imax = config["imax"],
        istep = config["istep"]
    threads: 
        threads = config["threads"]
    run:
        inflations = [round(x, 2) for x in np.arange(params.imin, params.imax, params.istep)]
        for inflation in inflations:
            shell(f'orthofinder -f {input.proteomes} -t {threads} -a {threads} -I {inflation} -n output_{inflation}')

#For parallel execution, execute the following command:
# snakemake --use-conda --cores 32 --snakefile run_orthofinder.py --configfile config.yaml