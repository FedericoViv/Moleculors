#' include Molecular_Input_v01.R Molecular_graph_V01.R GUI_V01.R zeroD_descriptor_V01.R twoD_descriptor_V01.R descriptor_launcher_V01.R

#' Moleculors environment
#' @export
#'

Moleculors <- new.env()

#' Matrices environment
#' @export
#'
Mol_mat = new.env()

#' Descriptors environment
#' @export
#'
Output_descp = new.env()



#' Moleculors initialization and environments set.
#'
#' Every functions will be part and only part
#' of this environment thus preventing global variables overlay. Graphic
#' matrices will belong to the Mol_mat envinroment and molecular descriptor to
#' Ouput_descp
#'
#' @param GUI Graphical user interface option. If true launch the GUI.
#'
#' @return if GUI = TRUE Moleculors environment is initializated with a GUI, else just set the environment
#'
#' @examples
#' Moleculors_init()
#'
#' @export

Moleculors_init = function(GUI = TRUE){

  assign("Moleculors", new.env(), envir = .GlobalEnv)
  assign("Mol_mat", new.env(), envir = .GlobalEnv)
  assign("Output_descp", new.env(), envir = .GlobalEnv)

  if (GUI == TRUE) {
    Moleculors_GUI()
  }
}
