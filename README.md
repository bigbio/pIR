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

Audain, E., Ramos, Y., Hermjakob, H., Flower, D. R., & Perez-Riverol, Y. (2015). Accurate estimation of isoelectric point of protein and peptide based on amino acid sequences. Bioinformatics, btv674. [article](http://bioinformatics.oxfordjournals.org/content/early/2015/12/09/bioinformatics.btv674.full)

Perez-Riverol, Y., Audain, E., Millan, A., Ramos, Y., Sanchez, A., Vizcaíno, J. A., ... & González, L. J. (2012). Isoelectric point optimization using peptide descriptors and support vector machines. Journal of proteomics, 75(7), 2269-2274. [article](https://www.researchgate.net/profile/Vladimir_Besada/publication/221825414_Isoelectric_point_optimization_using_peptide_descriptors_and_support_vector_machines/links/09e41503561f3b0787000000.pdf)

### This library has been used in:

- Ramos, Y., Gutierrez, E., Machado, Y., Sánchez, A., Castellanos-Serra, L., González, L.J., Fernández-de-Cossio, J., Pérez-Riverol, Y., Betancourt, L., Gil, J. and Padrón, G., 2008. Proteomics based on peptide fractionation by SDS-free PAGE. Journal of proteome research, 7(6), pp.2427-2434. [article](https://www.researchgate.net/profile/Vladimir_Besada/publication/5431019_Proteomics_based_on_peptide_fractionation_by_SDS-free_PAGE/links/0912f50355c4ac1a82000000.pdf)

- Ramos, Y., Garcia, Y., Pérez‐Riverol, Y., Leyva, A., Padrón, G., Sánchez, A., Castellanos‐Serra, L., González, L.J. and Besada, V., 2011. Peptide fractionation by acid pH SDS‐free electrophoresis. Electrophoresis, 32(11), pp.1323-1326. [article](https://www.researchgate.net/profile/Yasset_Perez-Riverol/publication/51094357_Peptide_fractionation_by_acid_pH_SDS-free_electrophoresis/links/53dd6e730cf2a76fb667c9f4.pdf)

