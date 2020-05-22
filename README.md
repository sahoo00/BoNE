# BoNE
Boolean Network Explorer

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
stepminer="java -cp $CLASSPATH -Xms64m -Xmx10G tools.CustomAnalysis"
stepminer1="java -cp $CLASSPATH -Xms64m -Xmx10G tools.Analyze"
export PERL_HASH_SEED=0
```
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


