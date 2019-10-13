#Federico Molecular Descriptor collection


# Moleculors envinronments. Every functions will be part and only part
# of this environment thus preventing global variables overlay. Graphic
# matrices will belong to the Mol_mat envinroment and molecular descriptor to
# Ouput_descp

Moleculors_init = function(GUI = TRUE){

  assign("Moleculors", new.env(), envir = .GlobalEnv)
  assign("Mol_mat", new.env(), envir = .GlobalEnv)
  assign("Output_descp", new.env(), envir = .GlobalEnv)

  source("R/Molecular_Input_v01.R")
  source("R/Molecular_graph_V01.R")
  source("R/GUI_V01.R")
  source("R/zeroD_descriptor_V01.R")
  source("R/twoD_descriptor_V01.R")
  source("R/descriptor_launcher_V01.R")



  if (GUI == TRUE) {
    Moleculors_GUI()
  }
}
