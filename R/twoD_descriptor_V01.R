## This file contain the two dimentional molecular descriptor
# Each function take as input a graphical matrix calculated in the Molecular_graph source file.
# Each molecular descriptor will be saved in the Output_descp environment for future call


#' Moleculors Weiner indexes calculator
#'
#' This function take as input the vertex distance graphical matrix and the V distance distance matrix
#' and return the Wiener Indexes calculated from both
#' Wiener index is calculated as half of the summatory of the summatory of the distances in each
#' row of the graphical matrix.
#' average distance distance degree is calculated as the summatory of the summatory of the distances in each
#' row of the distance distance graphical matrix divided by the number of atoms.
#' D/D index is the defined as half of the summatory of the summatory of the distances in each
#' row of the graphical distance distance matrix. This last two give informations about
#' molecule folding.
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

  if (is.matrix(Mol_mat$graph_Vdistance_matrix) & is.matrix(Mol_mat$graph_Vdistance_distance_matrix)) {

    W = 0

    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {

      W = W + sum(Mol_mat$graph_Vdistance_matrix[i,])

    }

    ADDD = 0

    for (i in 1:nrow(Mol_mat$graph_Vdistance_distance_matrix)) {

      ADDD = ADDD + sum(Mol_mat$graph_Vdistance_distance_matrix[i,])

    }


    Output_descp$W_index = 0.5 * W
    Output_descp$ADDD_index = ADDD/nrow(Mol_mat$graph_Vdistance_distance_matrix)
    Output_descp$D_D_index = 0.5 * ADDD

    message("Wiener indexes ... OK")

  } else {
    message("Wiener indexes ... FAIL")
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
    Output_descp$mean_Estate <- sum(intrinsic_state)/nrow(Mol_mat$graph_Vadj_matrix)

    message("E-state index ... OK")
  } else {
    message("E-state index ... FAIL")
  }
}



#' Moleculors Molecular bulk indexes and Electronegativity indexes
#'
#' This function define the molecular bulk of the selected molecule as
#' the summatory of alpha, with alpha Z - Zv / Z * 1/ PN -1
#' where Z is the atomic number of the selected atom, Zv is the number of
#' valence electrons and PN is the periodic number of the atom.
#' The comparison of alpha summatory with a reference molecules
#' without heteroatoms allows the calculation of two other indexes
#' dalphaA = 'summatory alpha - summatory alpha reference' / Natoms without hydrogens
#' This index gives informations about the presence of heteroatoms
#' dalphaB = 'sumamtory alpha reference - summatory alpha' /Natoms without hydrogens
#' This index gives a rough count of hydrogen bond acceptor atoms and polar surface areas.
#'
#' This function defines also the electronegativity index and other correlated indexes.
#' epsilon = - alpha + 0.3* Zv Electronegativity of each atoms
#' epsilon1 = summatory epsilon / Nof atoms including hydrogens. A measure of electronegativity atoms count
#' epsilon2 = summatory of epsilon excluding hydrogens / Nof atoms excluding hydrogens
#' epsilon3 = summatory of epsilon reference / Nof reference atoms
#' epsilon4 = summatory of epsilon of the saturated carbon skeleton /Nof atoms saturated carbon skeleton
#' epsilon5 = summatory of epsilon excluding hydrogens + summatory of epsilon of hydrogens attached to heteroatoms / Nof atoms excluding hydrogens + Nof hydrogens attached to heteroatoms
#'
#' DepsilonA =epsilon1 - epsilon3 measure of contribution of unsaturation and electronegative atom count
#' DepsilonB = epsilon1 - epsilon4 measure of contribuion of unsaturation
#' DepsilonC = epsilon3 - epsilon4 a measure of contribuion of electronegativity
#' DepsilonD = epsilon2 - epsilon5 a measure of contribution of hydrogen-bond donor atoms
#'
#'
#' @return 13 different indexes regarding: bulk and electronic distribution.  Values are stored in Ouput_descp environment.
#'
#' @examples
#' Bulk_Electronegativity_indexes_calc()
#'
#' @export
#'

Bulk_Electronegativity_indexes_calc <- function() {
  if ( is.data.frame(Mol_mat$input_H_suppressed) & is.data.frame(Mol_mat$input) &
       is.matrix(Mol_mat$graph_Vlaplacian_full_matrix) & is.matrix(Mol_mat$graph_Vadj_matrix_full)) {

    alpha <- c()
    epsilon <- c()
    alphaR <- c()
    epsilonR <- c()
    Zval <- c()
    Ztot <- c()
    total_bulk <- 0
    PN <- c()
    ZvalR <- c()
    ZtotR <- c()
    total_bulkR <- 0
    PNR <- c()

    reference_alcane <- Mol_mat$input_H_suppressed

    for (i in 1:nrow(reference_alcane)){
      if (reference_alcane$Atom[i] != "C" ) {
        reference_alcane$Atom[i] = "C"
      }
    }




    valence_electrons <- read.csv("tables/valence_electrons_table.csv")

    for (v in 1:nrow(Mol_mat$input_H_suppressed)) {
      Ztot[v] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
      Zval[v] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
      PN[v] <- valence_electrons$Principal_quantum_number[which(valence_electrons$Symbol == as.character(Mol_mat$input_H_suppressed$Atom[[v]]))]
      ZtotR[v] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(reference_alcane$Atom[[v]]))]
      ZvalR[v] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(reference_alcane$Atom[[v]]))]
      PNR[v] <- valence_electrons$Principal_quantum_number[which(valence_electrons$Symbol == as.character(reference_alcane$Atom[[v]]))]

      alpha[v] <- ((Ztot[v] - Zval[v])/Zval[v])*(1/(PN[v] - 1))
      alphaR[v] <- ((ZtotR[v] - ZvalR[v])/ZvalR[v])*(1/(PNR[v] - 1))
      epsilon[v] <- -alpha[v] +0.3*Zval[v]
      epsilonR[v] <- -alphaR[v] +0.3*ZvalR[v]
    }

    alpha_full <- c()
    epsilon_full <- c()
    Zval_full <- c()
    Ztot_full <- c()
    PN_full <- c()
    alpha_sat <- c()
    epsilon_sat <- c()
    Zval_sat <- c()
    Ztot_sat <- c()
    PN_sat <- c()


    saturated_skeleton <- Mol_mat$input
    levels(saturated_skeleton$Atom) = c(levels(saturated_skeleton$Atom),"H")
    adding_hydrogen <- 0

    for (h in 1:nrow(saturated_skeleton)) {
      if (Mol_mat$graph_Vlaplacian_full_matrix[h,h] != 4 & as.character(saturated_skeleton$Atom[h]) != "H") {
        adding_hydrogen <- adding_hydrogen + 4 - Mol_mat$graph_Vlaplacian_full_matrix[h,h]
      }
    }

    for (l in 1:adding_hydrogen) {
      saturated_skeleton[nrow(saturated_skeleton) +1,] = c("H", 0, 0, 0)
    }

    for (j in 1:nrow(Mol_mat$input)) {
      Ztot_full[j] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input$Atom[[j]]))]
      Zval_full[j] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(Mol_mat$input$Atom[[j]]))]
      PN_full[j] <- valence_electrons$Principal_quantum_number[which(valence_electrons$Symbol == as.character(Mol_mat$input$Atom[[j]]))]
      if (as.character(Mol_mat$input$Atom[j]) == "H") {
        alpha_full[j] = 0
      } else {
        alpha_full[j] <- ((Ztot_full[j] - Zval_full[j])/Zval_full[j])*(1/(PN_full[j] - 1))
      }

      epsilon_full[j] <- -alpha_full[j] +0.3*Zval_full[j]
    }

    for (j in 1:nrow(saturated_skeleton)) {
      Ztot_sat[j] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == as.character(saturated_skeleton$Atom[[j]]))]
      Zval_sat[j] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == as.character(saturated_skeleton$Atom[[j]]))]
      PN_sat[j] <- valence_electrons$Principal_quantum_number[which(valence_electrons$Symbol == as.character(saturated_skeleton$Atom[[j]]))]

      if (as.character(saturated_skeleton$Atom[j]) == "H") {
        alpha_sat[j] = 0
      } else {
        alpha_sat[j] <- ((Ztot_sat[j] - Zval_sat[j])/Zval_sat[j])*(1/(PN_sat[j] - 1))
      }
      epsilon_sat[j] <- -alpha_sat[j] +0.3*Zval_sat[j]
    }


    hydrogen_to_heteroatom <- 0

    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix_full)) {
      if (as.character(Mol_mat$input$Atom[[i]]) == "H") {
        for (j in 1:ncol(Mol_mat$graph_Vadj_matrix_full)) {
          if (Mol_mat$graph_Vadj_matrix_full[i,j] == 1 & as.character(Mol_mat$input$Atom[[j]]) != "C" ) {
            hydrogen_to_heteroatom <- hydrogen_to_heteroatom + 1
          }
        }
      }
    }


    total_bulk <- sum(alpha)
    total_bulkR <- sum(alphaR)

    dalphaA <- (total_bulk - total_bulkR)/nrow(Mol_mat$input_H_suppressed)
    dalphaB <- (total_bulkR - total_bulk)/nrow(Mol_mat$input_H_suppressed)

    epsilon2 <- sum(epsilon)/nrow(Mol_mat$input_H_suppressed)
    epsilon3 <- sum(epsilonR)/nrow(Mol_mat$input_H_suppressed)
    epsilon1 <- sum(epsilon_full)/nrow(Mol_mat$input)
    epsilon4 <- sum(epsilon_sat)/nrow(saturated_skeleton)
    epsilon5 <- (sum(epsilon) + hydrogen_to_heteroatom*(0.3))/(nrow(Mol_mat$input_H_suppressed) + hydrogen_to_heteroatom)

    depsilonA <- epsilon1 - epsilon3
    depsilonB <- epsilon1 - epsilon4
    depsilonC <- epsilon3 - epsilon4
    depsilonD <- epsilon2 - epsilon5


    Output_descp$Bulk_index <- total_bulk
    Output_descp$dalphaA <- dalphaA
    Output_descp$dalphaB <- dalphaB
    Output_descp$epsilon1 <- epsilon1
    Output_descp$epsilon2 <- epsilon2
    Output_descp$epsilon3 <- epsilon3
    Output_descp$epsilon4 <- epsilon4
    Output_descp$epsilon5 <- epsilon5
    Output_descp$depsilonA <- depsilonA
    Output_descp$depsilonB <- depsilonB
    Output_descp$depsilonC <- depsilonC
    Output_descp$depsilonD <- depsilonD

    message("ETA indexes ... OK")
  } else {
    message("ETA indexes ... FAIL")
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


#' Moleculors eccentric connectivity index
#'
#' This function define the eccentric connectivity index of the input molecule as:
#' Ec_index = sum"csii*sigmai" where csii is the max topological distance from the atom
#' i and sigmai is is vertex degree
#'
#' @return eccentric connectivity index for the input molecule. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Ec_index_calc()
#'
#' @export
#'

Ec_index_calc <- function(){
  if (is.matrix(Mol_mat$graph_Vdistance_matrix) & is.matrix(Mol_mat$graph_Vlaplacian_matrix)) {
    csi <- 0
    sigma <- 0
    Ec <- 0

    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
      csi <- max(Mol_mat$graph_Vdistance_matrix[i,])
      sigma <- Mol_mat$graph_Vlaplacian_matrix[i,i]
      Ec <- Ec + csi*sigma
    }
    Output_descp$Ec_index <- Ec
    message("Eccentric connectivity index ... OK")
  } else {
    message("Eccentric connectivity index... FAIL")
  }
}

#' Moleculors adj matrix eigenvalues descriptors
#'
#' This function define several descriptors from the eigenvalue of the adjmatrix:
#' Epi represents the energy of the energy of the pi electrons of the molecule. It is define as the
#' summatory of the eigenvalue i times the occupation number of the MO i, this can be 0,1 or 2. Since usually
#' the orbital contains 2 electrons when the eigenvalue is positive and 0 when negative. We can redefine the Epi
#' value as the summatory of the absolute value of the eigenvalue. Since the number of eigenvalues is equal to the number of
#' atoms in the H depleted matrxi, if the number is ODD we will have 1 occupation number, this is automatically taken into account
#' when adding the absolute values. e.g. Epi = 2 x +y1 + 2 x +y2 + 1 x +y3 + 0 x -y4 0 x -y5
#' LEDM "leading value of distance matrix" represent the max value of the eigenvalue of the Vdistance matrix.
#' This descriptor can be a discriminant of molecules size in molecules series.
#' Sparse CRI "characteristic root index" represent the sum of the positive eigenvalues of the sparse csi
#' matrix. Such value is sensitive to the presence of heteroatoms in the molecule.
#' Extadj_Vdegree_index and Extadj_Vdegree_max_index represent respectively the sum of the eigenvalues and the
#' maximum value of the eigenvalues of the extended Vadj matrix for vertex degree. This descriptors are
#' usually well correlated with physico-chemical properties and biological activities of organic compounds.
#' Extadj_Veln_index and Extadj_Veln_max_index represent respectively the sum of the eigenvalues and the
#' maximum value of the eigenvalues of the extended Vadj matrix for vertex electronegativity. This descriptors are
#' usually well correlated with physico-chemical properties and biological activities of organic compounds and affected by
#' the presence of heteroatoms.
#' @return eigenvalues descriptors for the input molecule. Value is stored in Ouput_descp environment.
#'
#' @examples
#' eigenvalues_descp_cacl()
#'
#' @export
#'

eigenvalues_descp_cacl <- function(){
  if (is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Vdistance_matrix) &
      is.matrix(Mol_mat$graph_Vsparsecsi_matrix) & is.matrix(Mol_mat$graph_Extended_Vadj_degree_matrix) &
      is.matrix(Mol_mat$graph_Extended_Vadj_eln_matrix) & is.matrix(Mol_mat$graph_Vdistance_distance_matrix)) {

    eigenvalues_Vadj <- eigen(Mol_mat$graph_Vadj_matrix)$values
    if (round(sum(eigenvalues_Vadj)) == 0) {
      Epi <- sum(abs(eigenvalues_Vadj))
    } else {
      Epi <- NA
      message("Anormalities in pi energies calculation. Probable non equal distribution in molecular orbitals electrons")
    }
    eigenvalues_Vdist <- eigen(Mol_mat$graph_Vdistance_matrix)$values
    LEDM <- max(eigenvalues_Vdist)

    eigenvalues_Vsparcecsi <- eigen(Mol_mat$graph_Vsparsecsi_matrix)$value
    eigenvalues_Vsparcecsi <- eigenvalues_Vsparcecsi[which(eigenvalues_Vsparcecsi > 0)]
    sparse_CRI <- sum(eigenvalues_Vsparcecsi)


    eigenvalues_extendedVadj_degree <- eigen(Mol_mat$graph_Extended_Vadj_degree_matrix)
    Extadj_Vdegree <- sum(abs(eigenvalues_extendedVadj_degree$values))
    Extadj_Vdegree_max <- max(abs(eigenvalues_extendedVadj_degree$values))

    eigenvalues_extendedVadj_eln <- eigen(Mol_mat$graph_Extended_Vadj_eln_matrix)
    Extadj_Veln <- sum(abs(eigenvalues_extendedVadj_eln$values))
    Extadj_Veln_max <- max(abs(eigenvalues_extendedVadj_eln$values))

    eigenvalues_distance_distance_matrix <- eigen(Mol_mat$graph_Vdistance_distance_matrix)
    folding_degree <- max(eigenvalues_distance_distance_matrix$values)/nrow(Mol_mat$graph_Vdistance_distance_matrix)


    Output_descp$Epi_index <- Epi
    Output_descp$LEDM_index <- LEDM
    Output_descp$sparse_CRI_index <- sparse_CRI
    Output_descp$Extadj_Vdegree_index <- Extadj_Vdegree
    Output_descp$Extadj_Vdegree_max_index <- Extadj_Vdegree_max
    Output_descp$Extadj_Veln_index <- Extadj_Veln
    Output_descp$Extadj_Veln_max_index <- Extadj_Veln_max
    Output_descp$folding_degree_index <- folding_degree
    message("Eigenvalues indexes... OK")
  } else {
    message("Eigenvalues indexes... FAIL")
  }
}


#' Moleculors Q polarity index
#'
#' This function define the Q polarity index of the input molecule as:
#' Q = ' A^2 sumIrefskeleton/sumI^2 '
#' where A is the numbero of atoms in the H suppresed input, Irefskeleton is the
#' intrinsic state of each atom in the reference alcane chain and I is the intrisic state
#' in the actual molecule. Q lies between to boundary value, one of total polarity approximated
#' by the square of the number of atoms A  and one of zero polarity given by the sp3 skeleton
#'
#'
#' @return Q polarity index for the input molecule. Value is stored in Ouput_descp environment.
#'
#' @examples
#' Q_polarity_calc()
#'
#' @export
#'

Q_polarity_calc <- function(){
  if (is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Vadj_matrix_full) &
      is.data.frame(Mol_mat$input) & is.data.frame(Mol_mat$input_H_suppressed) &
      is.matrix(Mol_mat$graph_Vdistance_matrix) & is.matrix(Mol_mat$graph_Vlaplacian_matrix)) {

    d = c()
    Zval = c()
    Ztot = c()
    H = c()
    intrinsic_state = c()
    ref_intrinsic_state = c()
    di = c(rep(0, nrow(Mol_mat$graph_Vadj_matrix)))

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


    saturated_skeleton <- Mol_mat$input
    levels(saturated_skeleton$Atom) <- c(levels(saturated_skeleton$Atom),"H")


    for (i in 1:nrow(saturated_skeleton)) {
      if (!(Mol_mat$input$Atom[i] %in% c("C", "H"))) {
        saturated_skeleton$Atom[i] = "C"
      }
      if (Mol_mat$graph_Vlaplacian_full_matrix[i,i] != 4 & as.character(saturated_skeleton$Atom[i]) != "H") {
        H[i] <- H[i] + 4 - Mol_mat$graph_Vlaplacian_full_matrix[i,i]
      }
    }

    for (h in 1:nrow(Mol_mat$input)) {
      if (Mol_mat$input$Atom[h] != "H") {
        H[h] <- H[h] + sum(length(intersect(which(Mol_mat$graph_Vadj_matrix_full[,h] == 1), which(Mol_mat$input$Atom == "H"))))
      }
    }

    if (sum(is.na(H)) >= 1) {
      H = H[- which(is.na(H))]
    }

    for (v in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      Ztot[v] <- valence_electrons$Total_electrons[which(valence_electrons$Symbol == "C")]
      Zval[v] <- valence_electrons$valence_electrons[which(valence_electrons$Symbol == "C")]
    }


    for (k in 1:length(Zval)) {
      d[k] <- (Zval[k] - H[k])/(Ztot[k] - Zval[k] - 1)
    }


    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      ref_intrinsic_state[i] <- (((2/valence_electrons$Principal_quantum_number[which(valence_electrons$Symbol == "C")])^2) *
                               d[i] + 1)/4
    }



    Output_descp$Q_polarity_index <- (nrow(Mol_mat$input_H_suppressed)^2)*sum(ref_intrinsic_state)/(sum(intrinsic_state)^2)

    message("Q polarity index ... OK")
  } else {
    message("Q polarity index ... FAIL")
  }
}





##########################TO BE IMPLEMENTED#####################

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
#'

# function under development

MED_index_calc <- function(){
  if (is.data.frame(Mol_mat$input) & is.matrix(Mol_mat$graph_Vlaplacian_full_matrix) &
      is.matrix(Mol_mat$graph_Vadj_matrix) & is.matrix(Mol_mat$graph_Vadj_matrix_full) &
      is.data.frame(Mol_mat$input_H_suppressed) & is.matrix(Mol_mat$graph_Vdistance_matrix)){

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

    Nactr <- 0

    Ndon <- 0

    for (i in 1:length(Mol_mat$input)) {

      if (as.character(Mol_mat$input$Atom[i]) %in% c("Br", "Cl", "F", "I")) {
        Nactr <- Nactr + 1
      } else if (as.character(Mol_mat$input$Atom[i]) == "C"){
        if (bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[i])),2] - Mol_mat$graph_Vlaplacian_full_matrix[i,i] == 2) {
          Nactr <- Nactr + 1
        } else if (bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[i])),2] - Mol_mat$graph_Vlaplacian_full_matrix[i,i] == 0) {
          Ndon <- Ndon + 1
        }
      } else if (as.character(Mol_mat$input$Atom[i]) == "N"){
        if (bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[i])),2] - Mol_mat$graph_Vlaplacian_full_matrix[i,i] == 0 &
            ! "O" %in% Mol_mat$input$Atom[Mol_mat$graph_Vadj_matrix_full[i,]] |
            bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[i])),2] - Mol_mat$graph_Vlaplacian_full_matrix[i,i] == 1) {
          Ndon <- Ndon + 1
        } else {
          Nactr <- Nactr + 1
        }
      } else if (as.character(Mol_mat$input$Atom[i]) == "S"){
        if (bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[i])),2] - Mol_mat$graph_Vlaplacian_full_matrix[i,i] == 0){
          Ndon <- Ndon + 1
        } else {
          Nactr <- Nactr + 1
        }
      } else if (as.character(Mol_mat$input$Atom[i]) == "O")
        if (bonds_library[which(bonds_library$atoms == as.character(Mol_mat$input$Atom[i])),2] - Mol_mat$graph_Vlaplacian_full_matrix[i,i] == 0){
          Ndon <- Ndon + 1
        }
    }

    mED <- (Nconj + Ndon/2 - Nactr)/Nat

    Output_descp$mED <- mED


    message("mED index ... OK")
  } else {
    message("mED index ... FAIL")
  }
}

