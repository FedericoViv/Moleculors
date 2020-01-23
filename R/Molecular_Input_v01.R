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

    cartesian_coordinates[[i]] <- read.csv(data_names[i], header = FALSE)
    names(cartesian_coordinates[[i]]) = c("Atom", "X", "Y", "Z")

  }

  Mol_mat$input_list = cartesian_coordinates

}

#' prediction matrix loader.
#'
#' This function take as input the prediction vector to be used
#' for the training of the NN for the qsar study and store the information in
#' the prediction matrix
#'
#'
#'
#' @examples
#' prediction_parameter()
#'
#'
#' @export

prediction_parameter <- function(){

  prediction_matrix <- read.csv(file.choose(), header = FALSE)

  if (ncol(prediction_matrix) > 1) {
    return(message("Error in prediction matrix structure. COL > 1"))

  } else {

    assign("prediction_matrix", prediction_matrix, envir = globalenv())
  }
}


#' Moleculors multiple graph and descriptor calculator.
#'
#' This function loop through the multiple input selected and calculate
#' the graphical matrices and descriptor for each molecule. Results are stored
#' in the descriptor matrix where each column is a molecule and each row are the
#' orded descriptors
#'
#'
#'
#' @examples
#' moleculors_multiple_descriptor()
#'
#'
#' @export


moleculors_multiple_descriptor <- function() {

  if (length(Mol_mat$input_list) > 0 ) {

    descriptor_matrix = matrix(ncol = length(Mol_mat$input_list), nrow = 13)


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
