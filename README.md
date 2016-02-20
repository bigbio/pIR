[![Build Status](https://travis-ci.org/ypriverol/pIR.svg?branch=master)](https://travis-ci.org/ypriverol/pIR)

[pIR](https://github.com/ypriverol/pIR)
======

An [R package](https://github.com/ypriverol/pIR) to analyze the isoelectric point of peptides and proteins based on experimental values and predicted using different functions. The package provides an statistical framework to analyze the correlation between predicted and expeted values, and it can be use in other contexts.

### Installation  

First, we need to install `devtools`:  

    install.packages("devtools")
    library(devtools)
   
Then we just call  

    install_github("ypriverol/pIR")
    library(pIR)

##Examples
=================

```R

# Other different pKSets that can be selected: rodwell, emboss, lehninger, grimsley, patrickios, DtaSelect

library(pIR)
seq <- "GLPRKILCAIAKKKGKCKGPLKLVCKC"
pI <- pIIterative(sequence = seq, pkSetMethod = "solomon")
print(pI)

#The result will be 10.526

```

### How to cite

Enrique Audain, Yassel Ramos, Henning Hermjakob, Darren R. Flower, Yasset Perez-Riverol. Accurate estimation of Isoelectric Point of Protein and Peptide based on Amino Acid Sequences. Bioinformatics (2015), doi: 10.1093/bioinformatics/btv674
