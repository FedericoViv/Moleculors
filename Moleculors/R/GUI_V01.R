## This source contain the necessary function to build the GUI with gwidgets
# This code is not intended to be modified by the user.


# library loading

library(gWidgets2)
library(gWidgets2RGtk2)


# Gui initialization with loader and calculators.

Moleculors_GUI = function(){

  environment(calc_matrix) <- environment()
  environment(vis_descriptor) = environment()
  environment(file_dump) = environment()


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
      Moleculors$molecular_input()
      cartesian_output = gtable(Mol_mat$input, container = welcome_grp)
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

  calc_matrix()
  vis_descriptor()
  file_dump()


  status_bar = gstatusbar("", container = load_win)

}


# Molecular graphical matricex calculation and visualization


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

    }
  )
}

# Molecular descriptos calculation and visualization

vis_descriptor = function() {



  descriptor_grp = ggroup(horizontal = FALSE,
                          container = upload_grp)


  descriptor_calc_btn = gbutton(
    "Molecular descriptors calculation",
    container = descriptor_grp,
    handler = function(h, ...){
      Moleculors$descriptor_launcher()
      #descriptor_picker = gmenu(c(),
      #container = descriptor_grp)

    }
  )
}


# Output file saving options

file_dump = function(){

  dumper_grp = ggroup(horizontal = FALSE,
                          container = upload_grp)

  dumper_btn = gbutton(
    "Download output file",
    container = dumper_grp,
    handler = function(h, ...){

    }
  )
}


