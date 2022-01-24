# pacqproc
Post-acquisition processing pipelines for light-sheet microscopy

This script is for performing routine post-acquisition processing of raw images obtained from light-sheet microscopes. 
The pipeline is built on <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">Snakemake</a> and
each step relies on other scripts developed in Chung Lab.

## Installation
```
conda create -n pacqproc --file environment.yaml
```
This command will create an isolated conda environment named `pacqproc` and install all the dependencies listed in `environment.yaml`.

## Running the pipeline
### Configuration
Please edit `config.yaml` to configure the pipeline. You can specify which workflow to use, data path, channels, labels, and so forth. 


### Running Snakefile
```
snakemake -j
``` 
This will simply run all the necessary steps with all available cores. If you want to limit the number of cores, use `-c [#cores]` flag instead of `-j`. 

Please refer to <a href="https://snakemake.readthedocs.io/en/stable/executing/cli.html" target="_blank">this</a> for more complicated options.





