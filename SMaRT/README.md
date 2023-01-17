
# Macrophage Polarization Project: Machine Learning Identifies Signatures of Macrophage Reactivity and Tolerance that Predict Disease Outcomes.

## Details
* [scr](https://github.com/sahoo00/BoNE/blob/master/SMaRT/scr): A bash script to download all required Hegemon files from hegemon.ucsd.edu and create all required files to build the Boolean Implication Network using [analyze.pl](https://github.com/sahoo00/BoNE/blob/master/analyze.pl) available in the parent directory.
* [macrophage.ipynb](https://github.com/sahoo00/BoNE/blob/master/SMaRT/macrophage.ipynb): This Jupyter Notebook file contains python code for reproducing manuscript figures.
* [node-\*.txt](https://github.com/sahoo00/BoNE/blob/master/SMaRT/): The txt file contain a list of genes in each node (cluster) of the macrophage network.

## Requirements
List of the required tools and dependencies for this project:
* python 3.6
* python packages: matplotlib, networkx, seaborn, pandas, ...
* java  1.8.0 25
* perl v5.26.1


## How to run
To build the macrophage network from scratch follow the bash script provided:
```sh
$ ./scr
```

To reproduce the main figure panels follow the code in macrophage.ipynb.
You need to run the jupyter using following commands:
```sh
$ export PYTHONHASHSEED=0
$ jupyter-notebook --no-browser
```

## Questions?
Please send your questions to [dsahoo@ucsd.edu](mailto:dsahoo@ucsd.edu)
