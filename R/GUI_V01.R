## This source contain the necessary function to build the GUI with gwidgets
# This code is not intended to be modified by the user.


#' Moleculors GUI
#'
#' Basic GUI for moleculors functions User is not supposed to use this functions. GUI should be loaded using Moleculors_init
#'
#' @import gWidgets2
#' @import gWidgets2RGtk2
#'
#' @return Graphical user interface for Moleculors functions
#'
#' @examples
#' Moleculors_GUI()
#'
#' @export

Moleculors_GUI = function(){

  environment(calc_matrix) <- environment()
  environment(vis_descriptor) <- environment()
  environment(file_dump) <- environment()

  load_win = gwindow("Moleculors GUI")
  main_grp = ggroup(container = load_win, horizontal = FALSE)

  header = glayout(container = main_grp)
  body = glayout(container = main_grp)
  tail = glayout(container = main_grp)

  header[1,1] = glabel(
    "Welcome to the Moleculors loading service",
  )

  body[1,1] = gbutton(
    "Upload Cartesian Coordinates",
    handler = function(h, ...){
      text_show_load = capture.output(molecular_input(), type = "message")
      body[31:32,1:2] = glabel(text_show_load)
      body[2:30,1:2] = gtable(Mol_mat$input)
    }
  )

  decimal_place <- function()
  {
    unname(Sys.localeconv()["decimal_point"] == ",")
  }

  body[1,2] <- gcheckbox(
    text      = "Use comma for decimal place",
    checked   = decimal_place(),
  )

  calc_matrix()
  vis_descriptor()
  file_dump()


  status_bar = gstatusbar("", container = load_win)


}


# Molecular graphical matricex calculation and visualization


calc_matrix = function(){

  body[1,5] = gbutton(
    "Graphical matrices calculation",
    handler = function(h, ...){
      text_show_graph = capture.output(graphical_matrix(), type = "message")
      tail[2,30] = glabel(text_show_graph)
      picker_list = names(Mol_mat)
      body[2,5:10] = gcombobox(picker_list,
                            handler = function(h, ...){
                              body[3:20,5:20] = gtable(get(svalue(h$obj), envir = Mol_mat))
                            })
    }
  )
}

# Molecular descriptos calculation and visualization

vis_descriptor = function() {

  body[1,40] = gbutton(
    "Molecular descriptors calculation",
    handler = function(h, ...){
      text_show_desc = capture.output(descriptor_launcher(), type = "message")
      tail[2,60] = glabel(text_show_desc)
      picker_list = names(Output_descp)
      body[2,40] = gcombobox(picker_list,
                               handler = function(h, ...){
                                 body[3:20,40] = gtable(get(svalue(h$obj), envir = Output_descp))
                               })
    }
  )
}


# Output file saving options

file_dump = function(){

  tail[1,1] = gbutton(
    "Download output file",
    handler = function(h, ...){

    }
  )
}


