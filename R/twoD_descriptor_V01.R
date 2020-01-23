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
      is.data.frame(Mol_mat$input) & is.data.frame(Mol_mat$input_H_suppressed)) {
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
    if (sum(is.na(H)) >= 1) {
      H = H[- which(is.na(H))]
    }


    for (v in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      Ztot[v] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
      Zval[v] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
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


#' Moleculors Kappa shape index first order
#'
#' This function takes the vertex adjacency matrix and the edge adjacency matrix and
#' return the first order Kappa shape index. This parameters gives information about the
#' shape of the molecule. As K is calculated as Pmin * Pmax /P where Pmin = V -1
#' Pmax = V * 'V-1' and P is the number of paths of the molecules aka edges.
#' This parameter gives information about the numbers of cycles inside the molecule. The lower
#' K compared to V the more cycle are present.
#'
#'
#' @return First order kappa shape index. Value is stored in Ouput_descp environment.
#'
#' @examples
#' First_order_kappa_calc()
#'
#' @export
#'

First_order_kappa_calc <- function() {
  if (is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Eadj_matrix)) {

    Pmin = nrow(Mol_mat$graph_Vadj_matrix) - 1
    Pmax = (nrow(Mol_mat$graph_Vadj_matrix) * (nrow(Mol_mat$graph_Vadj_matrix) - 1))/2

    kappa1 = 2*(Pmin * Pmax)/(nrow(Mol_mat$graph_Eadj_matrix)* nrow(Mol_mat$graph_Eadj_matrix))

    Output_descp$kappa_first_index = kappa1

    message("First order Kappa shape index ... OK")
  } else {
    message("First order Kappa shape index ... FAIL")
  }
}


#' Moleculors Electrotopological state index
#'
#' This function requires access to the vertex adjacency matrices, full and suppressed, to the
#' raw inputs, full and suppressed, to the vertex laplacian matrix and to the vertex distance matrix
#' in order to compute the E-state index for the whole molecule. This is computed by calculating first
#' the intrinsic factor for each atom as 2/principal quanto number ^2 x valence degree + 1 / degree of the atom
#' and correcting this factor by a Di equal to the summatory of intrinsic state of atom i - intrinsic state of atom j
#' over the topological distanze r + 1^2. The obtained E-state index I + Di can be used for the single atoms of can be summed
#' to get the E-state of substituent or, as was done in this function, to compute the whole molecule E-state
#'
#' @return E-state index for the selected molecule. Value is stored in Ouput_descp environment.
#'
#' @examples
#' E_state()
#'
#' @export
#'

E_state <- function() {
  if (is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Vadj_matrix_full) &
      is.data.frame(Mol_mat$input) & is.data.frame(Mol_mat$input_H_suppressed) &
      is.matrix(Mol_mat$graph_Vdistance_matrix) & is.matrix(Mol_mat$graph_Vlaplacian_matrix)) {

    d = c()
    Zval = c()
    Ztot = c()
    H = c()
    intrinsic_state = c()
    di = c(rep(0, nrow(Mol_mat$graph_Vadj_matrix)))
    total_Estate = 0

    valence_electrons = read.csv("tables/valence_electrons_table.csv")
    for (h in 1:nrow(Mol_mat$input)) {
      if (Mol_mat$input$Atom[h] != "H") {
        H[h] <- sum(length(intersect(which(Mol_mat$graph_Vadj_matrix_full[,h] == 1), which(Mol_mat$input$Atom == "H"))))
      }
    }
    if (sum(is.na(H)) >= 1) {
      H = H[- which(is.na(H))]
    }

    for (v in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      Ztot[v] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
      Zval[v] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
    }


    for (k in 1:length(Zval)) {
      d[k] <- (Zval[k] - H[k])/(Ztot[k] - Zval[k] - 1)
    }


    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      intrinsic_state[i] <- (((2/valence_electrons$Principal_quantum_number[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[i]]))])^2) *
                                      d[i] + 1)/Mol_mat$graph_Vlaplacian_matrix[i,i]
    }

    for (y in 1:length(intrinsic_state)) {
      for (j in 1:length(intrinsic_state)) {
        di[y] <- di[y] + (intrinsic_state[y] - intrinsic_state[j])/(Mol_mat$graph_Vdistance_matrix[y,j] + 1)^2
      }
    }

    for (z in 1:length(intrinsic_state)) {
      total_Estate = total_Estate + intrinsic_state[z] + di[z]
    }

    Output_descp$total_Estate_index = total_Estate

    message("E-state index ... OK")
  } else {
    message("E-state index ... FAIL")
  }
}



#' Moleculors Molecular bulk index
#'
#' This function define the molecular bulk of the selected molecule as
#' the summatory of Z - Zv / Z * 1/ PN -1
#' where Z is the atomic number of the selected atom, Zv is the number of
#' valence electrons and PN is the periodic number of the atom.
#' The only required matrix is the input_H_suppressed as Hydrogen volume is taken
#' as reference and thus its value is set to zero.
#'
#' @return Molecular bulk for the selected molecule. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Bulk_index_calc()
#'
#' @export
#'

Bulk_index_calc <- function() {
  if ( is.data.frame(Mol_mat$input_H_suppressed)) {

    d <- c()
    Zval <- c()
    Ztot <- c()
    total_bulk <- 0
    PN <- c()

    valence_electrons <- read.csv("tables/valence_electrons_table.csv")

    for (v in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      Ztot[v] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
      Zval[v] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
      PN[v] <- valence_electrons$Principal_quantum_number[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
    }


    for (k in 1:length(Zval)) {
      d[k] <- ((Ztot[k] - Zval[k])/Zval[k])*(1/(PN[k] - 1))
    }

    total_bulk <- sum(d)

    Output_descp$Bulk_index = total_bulk

    message("Bulk index ... OK")
  } else {
    message("Bulk index ... FAIL")
  }
}

#' Moleculors hydrophilicity factor
#'
#' This function define the hydrophilicity factor of the input molecule as:
#' Hy = ""1+Nhy" *log2 "1+Nhy" + Nc "1/A log2 1/A" + sqrt"Nhy/A^2""/log"1+A"
#' where Nhy is the number of hydrogen attached to oxygen sulfur or nitrogen groups,
#' Nc is the number of carbons and A is the number of non hydrogen atoms
#'
#' @return Hydrophilicity factor for the input molecule. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Hydro_factor_calc()
#'
#' @export
#'

Hydro_factor_calc <- function(){
  if (is.data.frame(Mol_mat$input) & is.matrix(Mol_mat$graph_Vadj_matrix_full)) {

    Nc <- sum(length(which(Mol_mat$input$Atom == "C")))
    A <- sum(length(which(Mol_mat$input$Atom != "H")))
    Nhy <- 0

    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix_full)) {
      if (Mol_mat$input$Atom[i] == "O" | Mol_mat$input$Atom[i] == "S" | Mol_mat$input$Atom[i] == "N") {
        for (j in 1:ncol(Mol_mat$graph_Vadj_matrix_full)) {
          if (Mol_mat$graph_Vadj_matrix_full[i,j] == 1 & Mol_mat$input$Atom[j] == "H") {
            Nhy <- Nhy + 1
          }
        }
      }
    }

    Hy <- ((1 + Nhy)*log2(1 + Nhy) + Nc*((1/A)*log2(1/A)) + sqrt(Nhy/(A^2)))/log2(1 + A)

    Output_descp$Hydro_factor <- Hy

    message("Hydrophilicity factor ... OK")
  } else {
    message("Hydrophilicity factor ... FAIL")
  }

}


#' Moleculors Unsaturation Index
#'
#' This function define the unsaturation index as follow:
#' UI = log2 "1 + SBO + nB" where nD is the number of bonds
#' and SBO is the sum of bonds order in the molecule. The H suppressed matrix is used for
#' the computation.
#'
#' @return Unsaturation index of the input molecule. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Unsaturation_index_calc()
#'
#'
#'

##### function under development

Unsaturation_index_calc <- function() {

  if (is.data.frame(Mol_mat$input) & is.matrix(Mol_mat$graph_Vadj_matrix_full) &
      is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Eadj_matrix_full)) {

    nB <- nrow(Mol_mat$graph_Eadj_matrix_full)
    SBO <- 0


    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix_full)) {
      counter <- 0
      nDB <- 0
      nTB <- 0
      Bonds <- 0
      if (Mol_mat$input$Atom[i] != "H") {
        counter <- sum(Mol_mat$graph_Vadj_matrix_full[i,])
        if (Mol_mat$input$Atom[i] == "C" & 4-counter == 1) {
          nDB <- 1
        } else if (Mol_mat$input$Atom[i] == "C" & 4-counter == 2){
          nTB <- 2
        } else if (Mol_mat$input$Atom[i] %in% c("O","N","S") & 2-counter == 1){
          nDB <- 1
        }
        Bonds <- counter + nDB + nTB

        SBO <- SBO + Bonds/counter
      }

    }

    UI <- log2(1 + SBO + nB)


  }

}

#' Moleculors mean electronic distribution index
#'
#' This function define the mean electronic distribution as follow:
#' mED = 'number of coniugated double bonds/ number of atoms in the H suppresed matrix' +
#' 'number of electron donators groups'/number of the double bonds * max distance from the double bonds  -
#' 'number of electron actractor groups'/number of the double bonds * max distance from the double bonds
#' a donator group is: C-C/OH/OR/NH2/NHR/NR2/SH/SR while an actractor is:
#' NO2/CN/F/Cl/Br/I/C---C
#'
#' @return mean electronic distribution index. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Unsaturation_index_calc()
#'
#'
#'

MED_index_calc <- function(){
  if (is.data.frame(Mol_mat$input) & is.matrix(Mol_mat$graph_Vlaplacian_full_matrix) &
      is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Vadj_matrix_full) &
      is.matrix(Mol_mat$input_H_suppressed) & is.matrix(Mol_mat$graph_Vdistance_matrix)){

    bonds_library <- data.frame(atoms = c("H", "C", "N", "O", "S", "Cl", "Br", "F", "I"),
                                bonds = c(1, 4, 3, 2, 2, 1, 1, 1, 1))

    Nconj <- 0
    Nat <- nrow(Mol_mat$input_H_suppressed)
    dbonds <-0
    skiprow <- c()

    for (i in 1:nrow(Mol_mat$graph_Vlaplacian_full_matrix)) {
      if (bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[i])),2] - Mol_mat$graph_Vlaplacian_full_matrix[i,i] == 1 & !i %in% skiprow) {
        if (is.vector(skiprow)) {
          skiprow_memory <- append(skiprow_memory,skiprow)
        } else {
          skiprow_memory <- c()
        }
        holder <- c()
        skiprow <- c()
        dbonds <- dbonds + 1
        holder <- append(holder, which(Mol_mat$graph_Vlaplacian_full_matrix[i,] == -1))
        for (j in 1:length(holder)) {
          if ((bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[holder[j]])),2] - Mol_mat$graph_Vlaplacian_full_matrix[holder[j],holder[j]]) == 1) {
            skiprow <- append(skiprow, holder[j])
          }
        }
          Nconj <- length(unique(skiprow_memory))
      }
    }






    message("mED index ... OK")
  } else {
    message("mED index ... FAIL")
  }
}
