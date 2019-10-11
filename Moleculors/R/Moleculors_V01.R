#Federico Molecular Descriptor collection

source("R/Molecular_Input_v01.R")
source("R/Molecular_graph_V01.R")
source("R/GUI_V01.R")


# Moleculors envinronment. Every function of output will be part and only part
# of this environment thus preventing global variables overlay.

Moleculors = new.env()

file_uploader()
