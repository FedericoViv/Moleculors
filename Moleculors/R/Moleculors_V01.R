#Federico Molecular Descriptor collection


# Moleculors envinronments. Every functions will be part and only part
# of this environment thus preventing global variables overlay. Graphic
# matrices will belong to the Mol_mat envinroment and molecular descriptor to
# Ouput_descp

Moleculors = new.env()
Mol_mat = new.env()
Output_descp = new.env()

source("R/Molecular_Input_v01.R")
source("R/Molecular_graph_V01.R")
source("R/GUI_V01.R")
source("R/zeroD_descriptor_V01.R")
source("R/twoD_descriptor_V01.R")
source("R/descriptor_launcher_V01.R")




Moleculors_GUI()
