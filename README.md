# Gene-Caller

A set of tools to define gene boundaries *de novo* from only a nucleotide sequence. The basic caller works rather well for prokaryotic genomes shown below for the case of MG1655v3 (the most current genomic anotation for a common laboratory strain of *Escherichia coli*):

![alt text](https://github.com/adityaradhakrishnan/Gene-Caller/blob/master/MG1655-ROC-Curve.png "E. coli ROC Curve")

All genes are called with the exception of those with alternative start sites (GTG, TTG, CTG, ATT).
