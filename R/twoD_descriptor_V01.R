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

#' Moleculors Randic connectivity index calculator
#'
#' This function takes the Vertex adjacency matrix and the laplacian matrix
#' and return the randic connectivity index calculated as the summatory of the product
#' of the element ij of the adjacency matrix and the inverse of square root
#' of the product of the degree of the vertexes i j. NOTE this index can be considered
#' as the Kier and Hall index of the first order.
#'
#'
#' @return Randic connectivity index. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Randix_index_calc()
#'
#' @export
#'


Randic_index_calc <- function(){
  if (is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Vlaplacian_matrix)) {
    Randic = 0

    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      for (j in 1:ncol(Mol_mat$graph_Vadj_matrix)) {
        if (j > i) {
          Randic = Randic + Mol_mat$graph_Vadj_matrix[i,j]*(1/sqrt((Mol_mat$graph_Vlaplacian_matrix[i,i]*Mol_mat$graph_Vlaplacian_matrix[j,j])))
        }
      }
    }

    Output_descp$Randic_index = Randic

    message("Randic connectivity index ... OK")
  } else {
    message("Randic connectivity index ... FAIL")
  }
}


#' Moleculors Randic valence connectivity index calculator
#'
#' This function takes the Vertex adjacency matrix, the full adjacency matrix, the raw input and the input with H suppressed
#' and return the randic valence connectivity index calculated as the summatory of the product
#' of the valence corrected degree of each atom. This element calculated as Valence electrons - number of hydrogen over
#' total electrons of the atom - valence electrons - 1. In this way this element take into account the effect of heteroatoms
#' and multiple bonds.
#'
#'
#' @return Randic valence connectivity index. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Randix_valence_index_calc()
#'
#' @export
#'


Randic_valence_index_calc <- function(){
  if (is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Vadj_matrix_full) &
      is.data.frame(Mol_mat$input) & is.data.frame(input_H_suppressed)) {
    Randic_valence = 0
    d = c()
    Zval = c()
    Ztot = c()
    H = c()
    valence_electrons = read.csv("tables/valence_electrons_table.csv")
    for (h in 1:nrow(Mol_mat$input)) {
      if (Mol_mat$input$Atom[h] != "H") {
        H[h] <- sum(length(intersect(which(Mol_mat$graph_Vadj_matrix_full[,h] == 1), which(Mol_mat$input$Atom == "H"))))
      }
    }
    H = H[- which(is.na(H))]

    for (v in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      Ztot[v] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(input_H_suppressed$Atom[[v]]))] # total electrons
      Zval[v] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(input_H_suppressed$Atom[[v]]))] ##valence electrons
      }


    for (k in 1:length(Zval)) {
      d[k] <- (Zval[k] - H[k])/(Ztot[k] - Zval[k] - 1)
    }

    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      for (j in 1:ncol(Mol_mat$graph_Vadj_matrix)) {
        if (j > i) {
          Randic_valence = Randic_valence + Mol_mat$graph_Vadj_matrix[i,j]*(1/sqrt(d[i]*d[j]))
        }
      }
    }

    Output_descp$Randic_valence_index = Randic_valence

    message("Randic valence connectivity index ... OK")
  } else {
    message("Randic valence connectivity index ... FAIL")
  }
}

