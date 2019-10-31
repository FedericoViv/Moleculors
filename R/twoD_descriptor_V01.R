## This file contain the two dimentional molecular descriptor
# Each function take as input a graphical matrix calculated in the Molecular_graph source file.
# Each molecular descriptor will be saved in the Output_descp environment for future call


#' Moleculors Weiner index calculator
#'
#' This function take as input the vertex distance graphical matrix and return the Wiener Index
#' this index is calculated as half of the summatory of the summatory of the distances in each
#' row of the graphical matrix
#'
#'
#' @return Weiner index. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Wiener_index_calc()
#'
#' @export
#'


Wiener_index_calc = function(){

  if (is.matrix(Mol_mat$graph_Vdistance_matrix)) {

    W = 0

    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {

      W = W + sum(Mol_mat$graph_Vdistance_matrix[i,])

    }

    Output_descp$W_index = 0.5 * W

    message("Wiener index ... OK")

  } else {
    message("Wiener index ... FAIL")
  }
}


#' Moleculors Platt number calculator
#'
#' This function take as input the edge adjacency matrix and return the Platt Number
#' this number is calculated as the summatory of the summatory of each elementij in each row
#' of the adjacency matrix
#'
#'
#' @return Platt number. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Platt_number_calc()
#'
#' @export
#'


Platt_number_calc = function(){

  if (is.matrix(Mol_mat$graph_Eadj_matrix)) {

    Platt = 0

    for (i in 1:nrow(Mol_mat$graph_Eadj_matrix)) {

      Platt = Platt + sum(Mol_mat$graph_Eadj_matrix[i,])

    }

    Output_descp$Platt_number = Platt

    message("Platt number ... OK")

  } else {
    message("Platt number ... FAIL")
  }
}

#' Moleculors Zagreb index calculator
#'
#' This function take as input the Laplacian Vertex matrix and return
#' the Zagreb index for the selected molecular graph as the summatory
#' of the squared degree for each atom of the matrix
#'
#'
#' @return Zagreb index. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Zagreb_index_calc()
#'
#' @export
#'


Zagreb_index_calc = function(){

  if (is.matrix(Mol_mat$graph_Vlaplacian_matrix)) {
    Zagreb = 0

    for (i in 1:nrow(Mol_mat$graph_Vlaplacian_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Vlaplacian_matrix)) {
        if (i == j) {
          Zagreb = Zagreb + (Mol_mat$graph_Vlaplacian_matrix[i,j] * Mol_mat$graph_Vlaplacian_matrix[i,j])
        }
      }
    }

    Output_descp$Zagreb_index = Zagreb

    message("Zagreb index ... OK")
  } else {
    message("Zagreb index ... FAIL")
  }
}

#' Moleculors Balaban index calculator
#'
#' This function take as input the Vertex adjacency matrix, the Vertex
#' distance matrix,the Edge adjacency matrix and compute the
#' Balaban index as Edges/ 'Edges - Vertex + 1' x summatory 'total_distancei * total_distancej'^-0.5
#' for any adjacence vertex. The number 'Edges - Vertex + 1' is denominated with u and is
#' called cyclomatic number
#'
#'
#' @return Baladan index. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Balaban_index_calc()
#'
#' @export
#'


Balaban_index_calc = function(){

  if (is.matrix(Mol_mat$graph_Vdistance_matrix) & is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Eadj_matrix)) {
    u = nrow(Mol_mat$graph_Eadj_matrix) - nrow(Mol_mat$graph_Vadj_matrix) + 1

    J = 0

    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
        if (j > i & Mol_mat$graph_Vadj_matrix[i,j] == 1) {
          J = J + (sum(Mol_mat$graph_Vdistance_matrix[i,]) * sum(Mol_mat$graph_Vdistance_matrix[j,]))^-0.5
        }
      }
    }

    J = (nrow(Mol_mat$graph_Eadj_matrix)/(u + 1)) * J

    Output_descp$Balaban_index = J

    message("Balaban index ... OK")
  } else {
    message("Balaban index ... FAIL")
  }
}

