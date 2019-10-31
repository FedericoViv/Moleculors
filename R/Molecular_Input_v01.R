#' Moleculors input file loader.
#'
#' This function take a csv file containing 4 colums:
#' atom symbol, X,Y,Z. Any difference will be detected as
#' an error or warning depending on the situation. NOTE: the input
#' coordinates are supposed to be the optimized coordinates of the
#' target molecule! If this requirement is not fullfilled any descriptor
#' will not be related to the topological property at all.
#'
#' @name Moleculors$molecular_input
#'
#' @usage Moleculors$molecular_input()
#'
#' @return dataframe named Input inside Mol_mat environment to be used for matrices and descriptors calculation
#'
#' @examples
#' Moleculors$molecular_input()
#'
#' @export

Moleculors$molecular_input = function(){

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

