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

Audain, E., Ramos, Y., Hermjakob, H., Flower, D. R., & Perez-Riverol, Y. (2015). Accurate estimation of isoelectric point of protein and peptide based on amino acid sequences. Bioinformatics, btv674. (http://bioinformatics.oxfordjournals.org/content/early/2015/12/09/bioinformatics.btv674.full)[article]

Perez-Riverol, Y., Audain, E., Millan, A., Ramos, Y., Sanchez, A., Vizcaíno, J. A., ... & González, L. J. (2012). Isoelectric point optimization using peptide descriptors and support vector machines. Journal of proteomics, 75(7), 2269-2274. (https://www.researchgate.net/profile/Vladimir_Besada/publication/221825414_Isoelectric_point_optimization_using_peptide_descriptors_and_support_vector_machines/links/09e41503561f3b0787000000.pdf)[article]




