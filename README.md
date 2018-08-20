# sMD-PF-M1AAP
This repo containts the input files and scripts used in the paper: "STEERED MOLECULAR DYNAMIC SIMULATIONS REVEAL CRITICAL RESIDUES FOR (UN)BINDING OF SUBSTRATES, INHIBITORS AND A PRODUCT OF THE MALARIAL PFM1AAP" by Moore et.al.

### Additional Info

This package is currently in active developemnt. If additional help is needed in creating the Ligand Interaction Networks please contact the author. Additionaly, please contact myself or raise an issue if you detect any bugs with the library. 

### Install Notes

To install this package:

1) download or clone the repo
2) change to the directory
3) in an R console type:
       install.packages("./smdanalysis/", repos = NULL, type = "source")

### Creating a Ligand Occupancy Map

library(smdanalysis)

pmat <- matrix(rep(0, 16), ncol = 4)
diag(pmat) <- 1
lom <- lom_matrix(directory = "~/Documents/xyz_files/", pmat = view_mat, spacing = 75)
 
plot_lom(lom)
