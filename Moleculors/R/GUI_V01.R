## This source contain the necessary function to build the GUI with gwidgets
# This code is not intended to be modified by the user.


# library loading

library(gWidgets2)
library(gWidgets2RGtk2)


# Gui initialization with loader and calculators.

file_uploader = function(){

  load_win = gwindow("Moleculors GUI")

  welcome_grp = ggroup(container = load_win, horizontal = FALSE)

  lbl_welcome_msg = glabel(
    "Welcome to the Moleculors loading service",
    container = welcome_grp
  )

  upload_grp = ggroup(container = welcome_grp)

  upload_btn = gbutton(
    "Upload Cartesian Coordinates",
    container = upload_grp,
    handler = function(h, ...){
      Moleculors$Molecular_input()
      cartesian_output = gtable(Moleculors$Input, container = welcome_grp)
      calc_matrix()

    }
  )

  decimal_place <- function()
  {
    unname(Sys.localeconv()["decimal_point"] == ",")
  }

  chk_eurostyle <- gcheckbox(
    text      = "Use comma for decimal place",
    checked   = decimal_place(),
    container = upload_grp
  )

  status_bar = gstatusbar("", container = load_win)

  #oD descriptors to be calculated with the cartesian loading
}


#Molecular graphical matricex calculation and visualization


calc_matrix = function(){
  graph_matrices_grp = ggroup(horizontal = FALSE,
                              container = upload_grp)
  matrices_calc_btn = gbutton(
    "Graphical matrices calculation",
    container = graph_matrices_grp,
    handler = function(h, ...){
      Moleculors$graphical_matrix()
      #matrix_picker = gmenu(c(),
                             #container = graph_matrices_grp)
      vis_descriptor()
    }
  )
}

# Molecular descriptos calculation and visualization

vis_descriptor = function() {
  descriptor_grp = ggroup(horizontal = FALSE,
                          container = upload_grp)
  descriptor_calc_btn = gbutton(
    "2D_Molecular descriptors calculation",
    container = descriptor_grp,
    handler = function(h, ...){
      # calculation function
      #descriptor_picker = gmenu(c(),
                            #container = descriptor_grp)
    }
  )
}
