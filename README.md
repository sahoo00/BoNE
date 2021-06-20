# BoNE
Boolean Network Explorer

## Table of contents
* [General info](#general-info)
* [Systems requirements](#systems-requirements)
* [Installation guide](#installation-guide)
* [Demo](#demo)
* [Instruction for use](#instruction-for-use)

## General info
BoNE is a free, open source software
that can be used analyze biomedical datasets.
BoNE transforms Boolean implication relationships into simplified directed graph.
In the context of biological gene regulatory networks, the directed
graph posesses a new set of challenges. BoNE implements a set of
toolkit to explore and analyze the directed graph representation of biological
dataset and integrates with machine learning for computing predictive models.
BoNE simplify the Boolean implication netwok by clustering them first using
Boolean Equivalent relationships. The edges between clusters are
defined using overwhelming relationships observed between them. BoNE traverse
the graph to discover different directed paths and chose them using machine
learning framework to build predictive models.

## Systems requirements
Project is created with:
* java: 1.8.0 25
* perl: v5.26.1
* python: 3.6.9
  - lifelines: 0.22.8
  - matplotlib: 3.2.1
  - numpy: 1.19.4
* StepMiner: 1.1
* Hegemon: 1.1

## Installation guide

bash script to run the softwares:
scr

Prerequisite CPAN perl module:
Statistics::R

Software for Boolean analysis is available in github:
https://github.com/sahoo00/BooleanNet

A compiled version is provided in this directory: stepminer-1.1.jar

Software for Hegemon analysis is available in github:
https://github.com/sahoo00/Hegemon

Environment variables:
```
export JAVA_HOME="/booleanfs/sahoo/softwares/java/jdk1.8.0_45"
export PATH=$JAVA_HOME/bin:$PATH
export CLASSPATH="./stepminer-1.1.jar"
export HEGEMONPATH="/booleanfs2/sahoo/Hegemon"
stepminer="java -cp $CLASSPATH -Xms64m -Xmx10G tools.CustomAnalysis"
stepminer1="java -cp $CLASSPATH -Xms64m -Xmx10G tools.Analyze"
export PERL_HASH_SEED=0
```

Typical installation time on a normal computer varies from 5 mins to one hour.

## Demo

Building Boolean Implication Network:

```
FILE=ibd-network
rm -f $FILE.rl
${stepminer}
${stepminer} boolean bitMatrix $FILE.rl \
 peters-2017-ibd-bv.txt \
 $FILE.ph All 0.1 2.5 0.05
cp $FILE.rl $FILE-1.rl
${stepminer} boolean bitMatrixFill $FILE-1.rl
${stepminer} boolean bitMatrixFillStats $FILE-1.rl
${stepminer} boolean bitMatrixPrint $FILE-1.rl > $FILE-res.txt

FILE=ibd-network-2
rm -f $FILE.rl
${stepminer}
${stepminer} boolean bitMatrix $FILE.rl \
 arijs-2018-uc-bv.txt \
 $FILE.ph All 0.1 2.5 0.05
cp $FILE.rl $FILE-1.rl
${stepminer} boolean bitMatrixFill $FILE-1.rl
${stepminer} boolean bitMatrixFillStats $FILE-1.rl
${stepminer} boolean bitMatrixPrint $FILE-1.rl > $FILE-res.txt
```

Building Clustered Boolean Implication Network:

```
PERL_HASH_SEED=0 perl analyze.pl ibd eq > ibd-network-g-eq-4.txt
PERL_HASH_SEED=0 perl analyze.pl ibd eq-corr > ibd-network-g-eq-4-corr.txt
PERL_HASH_SEED=0 perl analyze.pl ibd cls > ibd-network-g-eq-cls-4.txt
PERL_HASH_SEED=0 perl analyze.pl ibd g-cls > ibd-network-g-eq-g-4.txt
```

Identifying Paths:

```
perl analyze.pl ibd genes > path-1.json
perl analyze.pl ibd genes C21orf33 NCF2 > path-2.json
```

Analysis of Paths:

python code with bone.py and comp.py

Expected runtimes for this job is 1 hour.

## Instruction for use

Instructions for use in Neuroblastoma,
Barrett's Esophagus, COVID-19 are provided in NB, BE and covid directories
respectively.

Jupyter Notebook is provided to reproduce the manuscript figures.


