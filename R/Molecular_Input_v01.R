#' Moleculors input file loader.
#'
#' This function take a csv file containing 4 colums:
#' atom symbol, X,Y,Z. Any difference will be detected as
#' an error or warning depending on the situation. NOTE: the input
#' coordinates are supposed to be the optimized coordinates of the
#' target molecule! If this requirement is not fullfilled any descriptor
#' will not be related to the topological property at all.
#'
#'
#' @examples
#' molecular_input()
#'
#'
#' @export

molecular_input = function(){

  cartesian_csv = tryCatch({ cartesian_csv = read.csv(file.choose(),
                                                      header = FALSE)
                      },
                    warning = function(w){
                        warning(w)
                        message("csv file doesn't look properly formatted")
                    },

                    error = function(e){
                        message("Input file doesn't look like a csv file")
                        return(NA)
                    },

                    finally = { message("Always use cartesian coordinates as input!")

                    })

  if (ncol(cartesian_csv) != 4) {

    return(message("Input file has more/less column than expected"))

  }
  names(cartesian_csv) = c("Atom", "X", "Y", "Z")

  print(cartesian_csv)

  Mol_mat$input = cartesian_csv

  return(message("Loading successful"))

}



#' Moleculors multiple input file loader.
#'
#' This function take a csv file containing 4 colums:
#' atom symbol, X,Y,Z. Any difference will be detected as
#' an error or warning depending on the situation. NOTE: the input
#' coordinates are supposed to be the optimized coordinates of the
#' target molecule! If this requirement is not fullfilled any descriptor
#' will not be related to the topological property at all.
#'
#' @import tcltk
#'
#' @examples
#' molecular_input_multiple()
#'
#'
#' @export


molecular_input_multiple <- function(){

  data_names <- tk_choose.files()

  cartesian_coordinates <- list()

  for (i in 1:length(data_names)) {

    cartesian_coordinates[[i]] <- read.csv(data_names[i])
    names(cartesian_coordinates[[i]]) = c("Atom", "X", "Y", "Z")

  }

  Mol_mat$input_list = cartesian_coordinates

}


moleculors_multiple_descriptor <- function() {

  if (length(Mol_mat$input_list) > 0 ) {

    descriptor_matrix = matrix(ncol = length(Mol_mat$input_list), nrow = 10)

    for (i in 1:length(Mol_mat$input_list)) {

      Mol_mat$input <- Mol_mat$input_list[[i]]

      molecular_weight()

      N_atoms()

      graphical_matrix()

      descriptor_launcher()

      descriptor_matrix[,i] <- t(as.matrix(Data_summary))


    }

    assign("descriptor_matrix", descriptor_matrix, envir = globalenv())

  }
}
