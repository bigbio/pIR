[![Build Status](https://travis-ci.org/ypriverol/pIR.svg?branch=master)](https://travis-ci.org/ypriverol/pIR)

[pIR](https://github.com/ypriverol/pIBenchmark)
======

An [R package](https://github.com/ypriverol/pIR) to analyze the isoelectric point of peptides and proteins based on experimental values and predicted using different functions. The package provides an statistical framework to analyze the correlation between predicted and expeted values, and it can be use in other contexts.

### Installation  

First, we need to install `devtools`:  

    install.packages("devtools")
    library(devtools)
   
Then we just call  

    install_github("ypriverol/pIBenchmark")
    library(prideR)

##Examples
=================

```R
 seq <- "GLPRKILCAIAKKKGKCKGPLKLVCKC"

pI <- pIIterative(sequence = seq, pkSetMethod = "solomon")
print(pI)

#The result will be 10.526

# Other different pKSets that can be selected: rodwell, emboss, lehninger, grimsley, patrickios, DtaSelect

```

### How to cite

* Enrique Audain Martinez, Yassel Ramos, Henning Hermjakob, Darren R. Flower and Yasset Perez-Riverol (2015). Comparative Benchmarking Methods Estimating Isoelectric Point (pI) of Peptide and Proteins from Amino Acid Sequence. Summitted to Bioinformatics.   
