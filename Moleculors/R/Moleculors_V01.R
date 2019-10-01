#Federico Molecular Descriptor collection

source("R/Molecular_Input_V01")
source("R/Molecular_graph_V01")


# Moleculors envinronment. Every function of output will be part and only part
# of this environment thus preventing global variables overlay.

Moleculors = new.env()
