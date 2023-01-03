
# Pro-differentiation Project: AI-guided identification and validation of a target to differentiate colorectal cancers?

## Details
* [scr](https://github.com/sahoo00/BoNE/blob/master/Prodiff/scr): A bash script to download all required Hegemon files from Hegemon.ucsd.edu and create all required files to build the Boolean Implication Network (in cells #3-8 in the Prodiff-Paper.ipynb) using [analyze.pl](https://github.com/sahoo00/BoNE/blob/master/analyze.pl) available in the parent directory.
* [explore.conf](https://github.com/sahoo00/BoNE/blob/master/Prodiff/explore.conf): A configuration file with datasets information (including file paths).
* [Prodiff-Paper.ipynb](https://github.com/sahoo00/BoNE/blob/master/Prodiff/Prodiff-Paper.ipynb): This Jupyter Notebook file contains python code for reproducing manuscript figures.
* [node-*.txt](https://github.com/sahoo00/BoNE/blob/master/Prodiff/): The txt outputs from the .ipynb file contain a list of genes in each cluster which can be found in Fig. S1-1 (E).

## Requirements
List of the required tools and dependencies for this project:
* python 3.6
* python packages: matplotlib, networkx, seaborn, pandas, ...
* java  1.8.0 25
* perl v5.26.1
* Hegemon
* StepMiner
* BoNE


## How to run
To setup the requirements run the following script:
```sh
$ ./scr
```

Edit the configuration file to add your files' paths. After that, you only need to run the Prodiff-Paper.ipynb to reproduce the important paper's figures. 

## Questions?
Please send your questions to [dsahoo@health.ucsd.edu](mailto:dsahoo@health.ucsd.edu) or  [sataheri@ucsd.edu](mailto:sataheri@ucsd.edu).
