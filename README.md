This repository contains the MATLAB codes used to generate astrocytes in the following manuscript:
Estimating the glutamate transporter surface density in mouse hippocampal astrocytes
Anca R. Rădulescu, Cassandra L. Williams,  Gabrielle C. Todd,  Alex A. Lemus, Haley E. Chesbro,  Annalisa Scimemi

doi: https://doi.org/10.1101/2021.05.08.443234

https://www.biorxiv.org/content/10.1101/2021.05.08.443234v1

The file named "cell_3D_octahedron" is the main MATLAB code. When run, it accesses the other files, and simulates a cell. It produces a figure with the 3D structure of an astrocyte with and branches color-coded according to their branching level. The code also calculates cell measures (Sholl profile, number of branches, average branch lengths, etc, as specified in the line by line comments). The code creates up to P=7 primary branches. To change the number of primary branches, change the value of P accordingly, and comment out the unused branches.
