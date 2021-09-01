# Moleculors

This repository contains a set of function for the calculation of molecular descriptor to be used for QSER analysis. The package requires a cartesian set of coordinates defining the target molecule structure. In this application, the cartesian file generated from AVOGADRO (https://avogadro.cc/) was employed. From such coordinates, graphical matrices are calculated and used for the computation of 0D, 1D, and 2D molecular descriptor to be used for later analysis. The package has been extensively tested for small molecule. Due to the use of cartesian coordinate as spatial input for the target molecule, it is recommended to check the distance matrix to avoid any error coming from twisted three-dimensional structures where bond distances may induce false bond information in graphical matrices due to the relative distance in three-dimensional space. 

The following material can be consulted for additional information regarding graphical matrices and molecular descriptors:
Janežič, D., 2015. Graph-Theoretical Matrices in Chemistry. Taylor & Francis Group.
Todeschini, R., Consonni, V., 2000. Handbook of Molecular Descriptors, 1st ed, Methods and Principles in Medicinal Chemistry. Wiley. https://doi.org/10.1002/9783527613106

