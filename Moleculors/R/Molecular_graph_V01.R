# This function take the cartesian coordinates of a molecule and return the graphical matrices
# of its structure. To do this first of all remove all the hydrogens leaving only the carbons
# and heteroatoms. Then it compute the relative distance for each atom from each other and
# round it up. In the end we will have a matrix n*n where n is the number of atoms
# and each cell contain an integer pointing to the number of bonds intercurring
# between the diagonal element and the others. The adjacency matrix is then computed by checking if
# the distance is 1 or not.
# For the edge matrix, the edge coordinate is computed from the input coordinates as the mid point
# between two adjacent atoms and then the distance and adjacency is computed using the previuos
# algorithm


Moleculors$graphical_matrix = function(){

  if(is.data.frame(Moleculors$Input)){

    hydrogen_vector = c()

    for (i in 1:nrow(Moleculors$Input)) {

      if (as.character(Moleculors$Input$Atom[i] == "H")) {

        hydrogen_vector = append(hydrogen_vector, i)

      }
      if (i == nrow(Moleculors$Input)) {

        Input_H_suppressed = Moleculors$Input[-hydrogen_vector,]

      }
    }


    Moleculors$Vmatrices(Input_H_suppressed)

    Moleculors$Ematrices(Input_H_suppressed, Moleculors$graph_Vadj_matrix)

    Moleculors$V_Harary_matrix(Moleculors$graph_Vdistance_matrix)

    Moleculors$E_Harary_matrix(Moleculors$graph_Edistance_matrix)

    Moleculors$Rec_VCdistance(Moleculors$graph_VCdistance_matrix)

    Moleculors$Rec_ECdistance(Moleculors$graph_ECdistance_matrix)

    Moleculors$Complementary_Vdistance(Moleculors$graph_Vdistance_matrix)

    Moleculors$Complementary_Edistance(Moleculors$graph_Edistance_matrix)

    Moleculors$Rec_Complementary_Vdistance(Moleculors$graph_VCOMPdistance_matrix)

    Moleculors$Rec_Complementary_Edistance(Moleculors$graph_ECOMPdistance_matrix)

    Moleculors$Distance_path_Vdistance(Moleculors$graph_Vdistance_matrix)

    Moleculors$Distance_path_Edistance(Moleculors$graph_Edistance_matrix)


  } else {

    message("No Input file detected")

    return("No graphical matrices were computed")

  }

  return(message("Computing graphical Matrices... OK"))

}



# This function calculate the V matrices using as Input the Hydrogen suppressed
# cartesian matrix. In order to compute the distances it starts by calculating the
# magnitude of the vector from each atom to one another. Then after rounding (too many digits
# may lead to bad calculation in the following steps) it choose the smallest value (different
# from 0) has the 1 distance value. A loop is used to assure that every value of lenght 1 has
# the same value for the normalization step following later.
# a normalization is the computed and the values are rounded to the upper interger (has the module
# will be for sure lesser than the real distance). This algoright should be usable for small
# molecule without many folding. It has been tested for small branched molecule and some cyclic one.
# What needs to be tested: heteroatoms may change lenght distance, branched molecule with atoms > 20 - 30

Moleculors$Vmatrices = function(Cart_Input_Hsupp){

  graph_Vdistance_matrix = matrix(nrow = nrow(Cart_Input_Hsupp), ncol = nrow(Cart_Input_Hsupp))
  graph_Vadj_matrix = matrix(nrow = nrow(Cart_Input_Hsupp), ncol = nrow(Cart_Input_Hsupp))
  graph_VCdistance_matrix = matrix(nrow = nrow(Cart_Input_Hsupp), ncol = nrow(Cart_Input_Hsupp))

  for (i in 1:nrow(graph_Vdistance_matrix)) {

    for (j in 1:ncol(graph_Vdistance_matrix)) {

      graph_Vdistance_matrix[i,j] = sqrt((Cart_Input_Hsupp$X[j] - Cart_Input_Hsupp$X[i])^2 +
                                           (Cart_Input_Hsupp$Y[j] - Cart_Input_Hsupp$Y[i])^2 +
                                           (Cart_Input_Hsupp$Z[j] - Cart_Input_Hsupp$Z[i])^2)





    }
  }

  graph_Vdistance_matrix = apply(graph_Vdistance_matrix, 2, round, 2)

  for (i in 1:nrow(graph_Vdistance_matrix)) {
    for (j in 1:nrow(graph_Vdistance_matrix)) {
      if ((graph_Vdistance_matrix[i,j] - min(graph_Vdistance_matrix[1,-1])) <= 0.05 & i != j) {
        graph_Vdistance_matrix[i,j] = min(graph_Vdistance_matrix[1,-1])
      }
    }
  }

  graph_Vdistance_matrix = ceiling(apply(graph_Vdistance_matrix, 2, `/`, min(graph_Vdistance_matrix[1,-1])))
  graph_VCdistance_matrix = abs(apply(graph_Vdistance_matrix, 2, `-`, nrow(graph_Vdistance_matrix)))

  for (i in 1:nrow(graph_Vdistance_matrix)) {
    for (j in 1:ncol(graph_Vdistance_matrix)) {
      if (i == j) {
        graph_Vadj_matrix[i,j] = 0
        graph_VCdistance_matrix[i,j] = 0
      } else if (graph_Vdistance_matrix[i,j] == 1){
        graph_Vadj_matrix[i,j] = 1
      } else {
        graph_Vadj_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Vdistance_matrix = graph_Vdistance_matrix
  Moleculors$graph_Vadj_matrix = graph_Vadj_matrix
  Moleculors$graph_VCdistance_matrix = graph_VCdistance_matrix

}



# This function require the H suppressed cartesian matrix
# and the adjagency matrix for the vertexes in order to compute
# the edge cartesian matrix and the relation edges matrices.
# each matrix will be save in the moleculors scope to be used for
# future descriptor calculation.


Moleculors$Ematrices = function(Cart_Input_Hsupp, Vadjmatrix){
  edge_matrix = matrix(nrow = nrow(Cart_Input_Hsupp),
                       ncol = (ncol(Cart_Input_Hsupp)-1))
  h = 1
  counter = 1

  for (i in 1:nrow(Vadjmatrix)) {
    for (j in h:nrow(Vadjmatrix)) {
      if (Vadjmatrix[i,j] == 1) {
        edge_matrix[counter, 1] = (Cart_Input_Hsupp$X[i] + Cart_Input_Hsupp$X[j])/2
        edge_matrix[counter, 2] = (Cart_Input_Hsupp$Y[i] + Cart_Input_Hsupp$Y[j])/2
        edge_matrix[counter, 3] = (Cart_Input_Hsupp$Z[i] + Cart_Input_Hsupp$Z[j])/2
        counter = counter + 1
      }
    }
    h = h + 1
  }

  graph_Edistance_matrix = matrix(nrow = nrow(edge_matrix), ncol = nrow(edge_matrix))


   for (i in 1:nrow(graph_Edistance_matrix)) {
     for (j in 1:ncol(graph_Edistance_matrix)) {
       graph_Edistance_matrix[i,j] = sqrt((edge_matrix[j,1] - edge_matrix[i,1])^2 +
                                                  (edge_matrix[j,2] - edge_matrix[i,2])^2 +
                                                  (edge_matrix[j,3] - edge_matrix[i,3])^2)
     }
   }

  NA_vector = c()

    for (j in 1:ncol(graph_Edistance_matrix)) {
      if (is.na(graph_Edistance_matrix[1,j]) & !(j %in% NA_vector)) {
        NA_vector = append(NA_vector, j)
      }
    }

  graph_Edistance_matrix = graph_Edistance_matrix[-NA_vector, -NA_vector]
  graph_Eadj_matrix = matrix(nrow = nrow(graph_Edistance_matrix), ncol = nrow(graph_Edistance_matrix))
  graph_ECdistance_matrix = matrix(nrow = nrow(graph_Edistance_matrix), ncol = nrow(graph_Edistance_matrix))

  for (i in 1:nrow(graph_Edistance_matrix)) {
    for (j in 1:nrow(graph_Edistance_matrix)) {
      if ((graph_Edistance_matrix[i,j] - min(graph_Edistance_matrix[1,-1])) <= 0.2 & i != j) {
        graph_Edistance_matrix[i,j] = min(graph_Edistance_matrix[1,-1])
      }
    }
  }

  graph_Edistance_matrix = ceiling(apply(graph_Edistance_matrix, 2, `/`, min(graph_Edistance_matrix[1,-1])))
  graph_ECdistance_matrix = abs(apply(graph_Edistance_matrix, 2, `-`, nrow(graph_Edistance_matrix)))

  for (i in 1:nrow(graph_Edistance_matrix)) {
    for (j in 1:ncol(graph_Edistance_matrix)) {
      if (i == j) {
        graph_Eadj_matrix[i,j] = 0
        graph_ECdistance_matrix[i,j] = 0
      } else if (graph_Edistance_matrix[i,j] == 1){
        graph_Eadj_matrix[i,j] = 1
      } else {
        graph_Eadj_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Edistance_matrix = graph_Edistance_matrix
  Moleculors$graph_Eadj_matrix = graph_Eadj_matrix
  Moleculors$graph_ECdistance_matrix = graph_ECdistance_matrix
}

# This function compute the vertex version of the so called
# Harary matrix. It take has input the Vdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$V_Harary_matrix = function(Vdistance_matrix){

  graph_VHarary_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))
  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (i != j) {
        graph_VHarary_matrix[i,j] = 1 /Vdistance_matrix[i,j]
      } else {
        graph_VHarary_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_VHarary_matrix = graph_VHarary_matrix

}

# This function compute the edge version of the so called
# Harary matrix. It take has input the Edistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$E_Harary_matrix = function(Edistance_matrix){

  graph_EHarary_matrix = matrix(nrow = nrow(Edistance_matrix), ncol = nrow(Edistance_matrix))
  for (i in 1:nrow(Edistance_matrix)) {
    for (j in 1:nrow(Edistance_matrix)) {
      if (i != j) {
        graph_EHarary_matrix[i,j] = 1 /Edistance_matrix[i,j]
      } else {
        graph_EHarary_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_EHarary_matrix = graph_EHarary_matrix

}

# This function compute the vertex version of the reciprocal
# complement vertex matrix . It take has input the VCdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$Rec_VCdistance = function(VCdistance_matrix){

  graph_RecVCdistance_matrix = matrix(nrow = nrow(VCdistance_matrix), ncol = nrow(VCdistance_matrix))
  for (i in 1:nrow(VCdistance_matrix)) {
    for (j in 1:nrow(VCdistance_matrix)) {
      if (i != j) {
        graph_RecVCdistance_matrix[i,j] = 1 /VCdistance_matrix[i,j]
      } else {
        graph_RecVCdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_RecVCdistance_matrix = graph_RecVCdistance_matrix
}




# This function compute the edge version of the reciprocal
# complement edge matrix . It take has input the ECdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the edge distance matrix
#

Moleculors$Rec_ECdistance = function(ECdistance_matrix){

  graph_RecECdistance_matrix = matrix(nrow = nrow(ECdistance_matrix), ncol = nrow(ECdistance_matrix))
  for (i in 1:nrow(ECdistance_matrix)) {
    for (j in 1:nrow(ECdistance_matrix)) {
      if (i != j) {
        graph_RecECdistance_matrix[i,j] = 1 /ECdistance_matrix[i,j]
      } else {
        graph_RecECdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_RecECdistance_matrix = graph_RecECdistance_matrix
}


# This function compute the vertex version of the complementary
# vertex matrix . It take has input the Vdistance_matrix and for each
# element of the matrix if i != j it returns the  value of dmin + dmax - elemij
# were dmin is the minimal lenght in the molecular graph (1 for simple molecules)
# dmax is the maximum lenght in the system and elemij is the value of the element i j
# of the vertex distance matrix
#


Moleculors$Complementary_Vdistance = function(Vdistance_matrix){

  graph_VCOMPdistance_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))
  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (i != j) {
        graph_VCOMPdistance_matrix[i,j] = 1 + max(Vdistance_matrix) - Vdistance_matrix[i,j]
      } else {
        graph_VCOMPdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_VCOMPdistance_matrix = graph_VCOMPdistance_matrix
}


# This function compute the edge version of the complementary
# edge matrix . It take has input the Edistance_matrix and for each
# element of the matrix if i != j it returns the  value of dmin + dmax - elemij
# were dmin is the minimal lenght in the molecular graph (1 for simple molecules)
# dmax is the maximum lenght in the system and elemij is the value of the element i j
# of the edge distance matrix
#


Moleculors$Complementary_Edistance = function(Edistance_matrix){

  graph_ECOMPdistance_matrix = matrix(nrow = nrow(Edistance_matrix), ncol = nrow(Edistance_matrix))
  for (i in 1:nrow(Edistance_matrix)) {
    for (j in 1:nrow(Edistance_matrix)) {
      if (i != j) {
        graph_ECOMPdistance_matrix[i,j] = 1 + max(Edistance_matrix) - Edistance_matrix[i,j]
      } else {
        graph_ECOMPdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_ECOMPdistance_matrix = graph_ECOMPdistance_matrix
}

# This function compute the vertex version of the reciprocal
# complementary vertex matrix . It take has input the Complementary_Vdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$Rec_Complementary_Vdistance = function(Complentary_Vdistance_matrix){

  graph_VRecCOMPdistance_matrix = matrix(nrow = nrow(Complentary_Vdistance_matrix), ncol = nrow(Complentary_Vdistance_matrix))
  for (i in 1:nrow(Complentary_Vdistance_matrix)) {
    for (j in 1:nrow(Complentary_Vdistance_matrix)) {
      if (i != j) {
        graph_VRecCOMPdistance_matrix[i,j] = 1 /Complentary_Vdistance_matrix[i,j]
      } else {
        graph_VRecCOMPdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_VRecCOMPdistance_matrix = graph_VRecCOMPdistance_matrix
}

# This function compute the edge version of the reciprocal
# complementary edge matrix . It take has input the Complementary_Edistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$Rec_Complementary_Edistance = function(Complentary_Edistance_matrix){

  graph_ERecCOMPdistance_matrix = matrix(nrow = nrow(Complentary_Edistance_matrix), ncol = nrow(Complentary_Edistance_matrix))
  for (i in 1:nrow(Complentary_Edistance_matrix)) {
    for (j in 1:nrow(Complentary_Edistance_matrix)) {
      if (i != j) {
        graph_ERecCOMPdistance_matrix[i,j] = 1 /Complentary_Edistance_matrix[i,j]
      } else {
        graph_ERecCOMPdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_ERecCOMPdistance_matrix = graph_ERecCOMPdistance_matrix
}

# This function compute the vertex version of the distance-path-matrix
# It take has input the Vdistance_matrix and for each
# element of the matrix if i != j it returns elemij*(elemij + 1)/2 value of
# the distance matrix
#

Moleculors$Distance_path_Vdistance = function(Vdistance_matrix){

  graph_Vpathdistance_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))
  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (i != j) {
        graph_Vpathdistance_matrix[i,j] = Vdistance_matrix[i,j]*((Vdistance_matrix[i,j] + 1)/2)
      } else {
        graph_Vpathdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Vpathdistance_matrix = graph_Vpathdistance_matrix
}

# This function compute the edge version of the distance-path-matrix
# It take has input the Edistance_matrix and for each
# element of the matrix if i != j it returns elemij*(elemij + 1)/2 value of
# the distance matrix
#

Moleculors$Distance_path_Edistance = function(Edistance_matrix){

  graph_Epathdistance_matrix = matrix(nrow = nrow(Edistance_matrix), ncol = nrow(Edistance_matrix))
  for (i in 1:nrow(Edistance_matrix)) {
    for (j in 1:nrow(Edistance_matrix)) {
      if (i != j) {
        graph_Epathdistance_matrix[i,j] = Edistance_matrix[i,j]*((Edistance_matrix[i,j] + 1)/2)
      } else {
        graph_Epathdistance_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Epathdistance_matrix = graph_Epathdistance_matrix
}


# This function return the laplacian matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return
# -1 if elemij = 1, sum(elemi, == 1) if i = j (which is the angle of the vertex)
# and 0 for every other element



Moleculors$Laplacian_Vdistance = function(Vdistance_matrix){

  graph_Vlaplacian_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))
  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (Vdistance_matrix[i,j] == 1) {
        graph_Vlaplacian_matrix[i,j] = - 1
      } else if (i == j){
        graph_Vlaplacian_matrix[i,j] = sum(Vdistance_matrix[i,] == 1)
      } else {
        graph_Vlaplacian_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Vlaplacian_matrix = graph_Vlaplacian_matrix
}


# This function return the laplacian matrix for the edge distance matrix
# it take as input the Edistance_matrix and for each element return
# -1 if elemij = 1, sum(elemi, == 1) if i = j (which is the angle of the edge)
# and 0 for every other element



Moleculors$Laplacian_Edistance = function(Edistance_matrix){

  graph_Elaplacian_matrix = matrix(nrow = nrow(Edistance_matrix), ncol = nrow(Edistance_matrix))
  for (i in 1:nrow(Edistance_matrix)) {
    for (j in 1:nrow(Edistance_matrix)) {
      if (Edistance_matrix[i,j] == 1) {
        graph_Elaplacian_matrix[i,j] = - 1
      } else if (i == j){
        graph_Elaplacian_matrix[i,j] = sum(Edistance_matrix[i,] == 1)
      } else {
        graph_Elaplacian_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Elaplacian_matrix = graph_Elaplacian_matrix
}



########################################TO BE IMPLEMENTED ################################


# This function return the detour vertex matrix, which contain the longest
# distance between two vertexes.

Moleculors$Detour_Vdistance = function(Vdistance_matrix){

}

# This function return the detour edges matrix, which contain the longest
# distance between two edges.

Moleculors$Detour_Edistance = function(Edistance_matrix){

}

# This function return the complement V matrix of the detour matrix as V - dij where V
# is the number of vertexes

Moleculors$Detour_VCdistance = function(Vdistance_matrix){

}


# This function return the complement E matrix of the detour matrix as E - dij where E
# is the number of edges

Moleculors$Detour_ECdistance = function(Edistance_matrix){

}

# This function return the reverse V matrix of the detour matrix as dmax - Vdetourij where dmax
# is the longest detour distance (max(Vdetour)) for vertexes

Moleculors$RevDetour_Vdistance = function(detour_Vdistance_matrix){

}

# This function return the reverse E matrix of the detour matrix as dmax - Edetourij where dmax
# is the longest detour distance (max(Edetour)) for edges


Moleculors$RevDetour_Edistance = function(detour_Edistance_matrix){

}

# This function contain information on both V detour and V distance matrices
# It returns Vdetourij if i < j and Vdistance if if i > j (i = j = 0)


Moleculors$Detour_Distance_Vdistance = function(Vdistance_matrix, detour_Vdistance_matrix){

}


# This function contain information on both E detour and E distance matrices
# It returns Edetourij if i < j and Edistance if if i > j (i = j = 0)


Moleculors$Detour_Distance_Edistance = function(Edistance_matrix, detour_Edistance_matrix){

}
