# Moleculors
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

This repository contains a set of function for the calculation of molecular descriptor to be used for QSER analysis. The package requires a cartesian set of coordinates defining the target molecule structure. In this application, the cartesian file generated from AVOGADRO (https://avogadro.cc/) was employed. From such coordinates, graphical matrices are calculated and used for the computation of 0D, 1D, and 2D molecular descriptor to be used for later analysis. The package has been extensively tested for small molecule. Due to the use of cartesian coordinate as spatial input for the target molecule, it is recommended to check the distance matrix to avoid any error coming from twisted three-dimensional structures where bond distances may induce false bond information in graphical matrices due to the relative distance in three-dimensional space. 

The following material can be consulted for additional information regarding graphical matrices and molecular descriptors:
- Janežič, D., 2015. Graph-Theoretical Matrices in Chemistry. CRC Press. https://doi.org/10.1201/b18389
- Todeschini, R., Consonni, V., 2000. Handbook of Molecular Descriptors, 1st ed, Methods and Principles in Medicinal Chemistry. Wiley. https://doi.org/10.1002/9783527613106
- Roy, K., and R. N. Das. "On some novel extended topochemical atom (ETA) parameters for effective encoding of chemical information and modelling of fundamental physicochemical properties." SAR and QSAR in Environmental Research 22.5-6 (2011): 451-472. https://doi.org/10.1080/1062936X.2011.569900
- Understanding the Basics of QSAR for Applications in Pharmaceutical Sciences and Risk Assessment. Roy, K., and R. N. Das. Academic Press <u>ISBN:978-0-12-801633-6<u>
