## This source file contain all the necessary functions to compute several
# graphical matrix necessary for the further computation of twoD descriptors (TOPOLOGICAL)
# The main function graphical_matrix() Calls every grpahical function and warn the user
# if any of the computation failed. Each matrix is then stored in the Moleculors environment
# for future use.





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

  if(is.data.frame(Mol_mat$input)){

    hydrogen_vector = c()

    for (i in 1:nrow(Mol_mat$input)) {

      if (as.character(Mol_mat$input$Atom[i] == "H")) {

        hydrogen_vector = append(hydrogen_vector, i)

      }
      if (i == nrow(Mol_mat$input)) {

        input_H_suppressed = Mol_mat$input[-hydrogen_vector,]

      }
    }


    Moleculors$Vadj_matrix(input_H_suppressed)

    Moleculors$Vdistance_matrix()

    Moleculors$VCdistance_matrix()

    Moleculors$Eadj_matrix(input_H_suppressed)

    Moleculors$Edistance_matrix()

    Moleculors$ECdistance_matrix()

    Moleculors$V_Harary_matrix()

    Moleculors$E_Harary_matrix()

    Moleculors$Rec_VCdistance()

    Moleculors$Rec_ECdistance()

    Moleculors$Complementary_Vdistance()

    Moleculors$Complementary_Edistance()

    Moleculors$Rec_Complementary_Vdistance()

    Moleculors$Rec_Complementary_Edistance()

    Moleculors$Distance_path_Vdistance()

    Moleculors$Distance_path_Edistance()

    Moleculors$Laplacian_Vdistance()

    Moleculors$Laplacian_Edistance()

    Moleculors$Vadj_kneigh_matrix(2)

    Moleculors$Vadj_highpower_matrix(2)



  } else {

    message("No Input file detected")

    return("No graphical matrices were computed")

  }

  return(message("Computing graphical Matrices... OK"))

}



# This function calculate the V adjacency matrix using as Input the Hydrogen suppressed
# cartesian matrix. In order to compute the distances it starts by calculating the
# magnitude of the vector from each atom to one another. Then after rounding (too many digits
# may lead to bad calculation in the following steps) it choose the smallest value (different
# from 0) has the 1 distance value. A loop is used to assure that every value of lenght 1 has
# the same value for the normalization step following later.
# a normalization is the computed and  every other value different from 1 is set to 0
#
#
#

Moleculors$Vadj_matrix = function(Cart_Input_Hsupp){

  graph_Vadj_matrix = matrix(nrow = nrow(Cart_Input_Hsupp), ncol = nrow(Cart_Input_Hsupp))

  for (i in 1:nrow(graph_Vadj_matrix)) {

    for (j in 1:ncol(graph_Vadj_matrix)) {

      graph_Vadj_matrix[i,j] = sqrt((Cart_Input_Hsupp$X[j] - Cart_Input_Hsupp$X[i])^2 +
                                      (Cart_Input_Hsupp$Y[j] - Cart_Input_Hsupp$Y[i])^2 +
                                      (Cart_Input_Hsupp$Z[j] - Cart_Input_Hsupp$Z[i])^2)

    }
  }

  graph_Vadj_matrix = apply(graph_Vadj_matrix, 2, round, 2)


  for (i in 1:nrow(graph_Vadj_matrix)) {
    for (j in 1:nrow(graph_Vadj_matrix)) {
      if ((graph_Vadj_matrix[i,j] - min(graph_Vadj_matrix[1,-1])) <= 0.1 & i != j) {
        graph_Vadj_matrix[i,j] = min(graph_Vadj_matrix[1,-1])
      }
    }
  }

  graph_Vadj_matrix = apply(graph_Vadj_matrix, 2, `/`, min(graph_Vadj_matrix[1,-1]))

  for (i in 1:nrow(graph_Vadj_matrix)) {
    for (j in 1:nrow(graph_Vadj_matrix)) {
      if (graph_Vadj_matrix[i,j] != 1) {
        graph_Vadj_matrix[i,j] = 0
      }
    }
  }

  Mol_mat$graph_Vadj_matrix = graph_Vadj_matrix

  message("Vertex adjacency matrix ... OK")
}



# This function take the adjacency matrix as input and generate a connection list where the index
# is related to the atom and the value in the index list are the atoms to which that atom is connected
# For each element with d != 1 is then calculated the distance by taking the value in the connection list
# copying those into the connection vector and if the interested atom (j) is not in the connection vector
# a new vector is created with all the connection of the atom to which the first atoms was connected
# this is looped increasing distance at each failed loop to detect the interested element.

Moleculors$Vdistance_matrix = function(){

  if (is.matrix(Mol_mat$graph_Vadj_matrix)) {
    graph_Vdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Vadj_matrix), ncol = nrow(Mol_mat$graph_Vadj_matrix))

    connections = list()
    index_counter = c()

    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      connections = append(connections, list(which(Mol_mat$graph_Vadj_matrix[i,] == 1)))
    }


    for (i in 1:nrow(graph_Vdistance_matrix)) {
      for (j in 1:nrow(graph_Vdistance_matrix)) {
        d = 1
        if (Mol_mat$graph_Vadj_matrix[i,j] == 1) {
          graph_Vdistance_matrix[i,j] = 1
        } else if (i == j) {
          graph_Vdistance_matrix[i,j] = 0
        } else if (i != j & Mol_mat$graph_Vadj_matrix[i,j] == 0) {
          connection_vector = connections[[i]]
          while (!(j %in% connection_vector)) {
            connections_holder = c()
            for (k in 1:length(connection_vector)) {
              connections_holder = append(connections_holder, connections[[connection_vector[k]]])
            }
            d = d + 1
            connection_vector = connections_holder
          }
          graph_Vdistance_matrix[i,j] = d
        }
      }
    }

    Mol_mat$graph_Vdistance_matrix = graph_Vdistance_matrix

    message("Vertex distance matrix ... OK")
  } else {
    message("Vertex distance matrix ... FAIL")
  }
}

# This function requires the Vdistance_matrix in order to compute the complement
# of the distance matrix. Each element is set has nrow(Vdistance) - Vdistanceij
#
#

Moleculors$VCdistance_matrix = function(){

  if (is.matrix(Mol_mat$graph_Vdistance_matrix)) {
    graph_VCdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Vdistance_matrix), ncol = nrow(Mol_mat$graph_Vdistance_matrix))

    for (i in 1:nrow(graph_VCdistance_matrix)) {
      for (j in 1:nrow(graph_VCdistance_matrix)) {
        if (i != j) {
          graph_VCdistance_matrix[i,j] = nrow(Mol_mat$graph_Vdistance_matrix) - Mol_mat$graph_Vdistance_matrix[i,j]
        } else {
          graph_VCdistance_matrix[i,j] = 0
        }
      }
    }
    Mol_mat$graph_VCdistance_matrix = graph_VCdistance_matrix

    message("Vertex complement distance matrix ... OK")
  } else {
    message("Vertex complement distance matrix ... FAIL")
  }
}


# This function require the H suppressed cartesian matrix
# and the adjagency matrix for the vertexes in order to compute
# the edge cartesian matrix and the adjacency edges matrices.
#
#


Moleculors$Eadj_matrix = function(Cart_Input_Hsupp){

  if (is.matrix(Mol_mat$graph_Vadj_matrix)) {
    edge_matrix = matrix(nrow = nrow(Cart_Input_Hsupp),
                         ncol = (ncol(Cart_Input_Hsupp)-1))
    h = 1
    counter = 1

    for (i in 1:nrow(Mol_mat$graph_Vadj_matrix)) {
      for (j in h:nrow(Mol_mat$graph_Vadj_matrix)) {
        if (Mol_mat$graph_Vadj_matrix[i,j] == 1) {
          edge_matrix[counter, 1] = (Cart_Input_Hsupp$X[i] + Cart_Input_Hsupp$X[j])/2
          edge_matrix[counter, 2] = (Cart_Input_Hsupp$Y[i] + Cart_Input_Hsupp$Y[j])/2
          edge_matrix[counter, 3] = (Cart_Input_Hsupp$Z[i] + Cart_Input_Hsupp$Z[j])/2
          counter = counter + 1
        }
      }
      h = h + 1
    }

    graph_Eadj_matrix = matrix(nrow = nrow(edge_matrix), ncol = nrow(edge_matrix))


    for (i in 1:nrow(graph_Eadj_matrix)) {
      for (j in 1:ncol(graph_Eadj_matrix)) {
        graph_Eadj_matrix[i,j] = sqrt((edge_matrix[j,1] - edge_matrix[i,1])^2 +
                                        (edge_matrix[j,2] - edge_matrix[i,2])^2 +
                                        (edge_matrix[j,3] - edge_matrix[i,3])^2)
      }
    }

    NA_vector = c()

    for (j in 1:ncol(graph_Eadj_matrix)) {
      if (is.na(graph_Eadj_matrix[1,j]) & !(j %in% NA_vector)) {
        NA_vector = append(NA_vector, j)
      }
    }

    if (!(is.null(NA_vector))) {
      graph_Eadj_matrix = graph_Eadj_matrix[-NA_vector, -NA_vector]
    }


    graph_Eadj_matrix = apply(graph_Eadj_matrix, 2, round, 2)

    for (i in 1:nrow(graph_Eadj_matrix)) {
      for (j in 1:nrow(graph_Eadj_matrix)) {
        if ((graph_Eadj_matrix[i,j] - min(graph_Eadj_matrix[1,-1])) <= 0.1 & i != j) {
          graph_Eadj_matrix[i,j] = min(graph_Eadj_matrix[1,-1])
        }
      }
    }

    graph_Eadj_matrix = apply(graph_Eadj_matrix, 2, `/`, min(graph_Eadj_matrix[1,-1]))

    for (i in 1:nrow(graph_Eadj_matrix)) {
      for (j in 1:nrow(graph_Eadj_matrix)) {
        if (graph_Eadj_matrix[i,j] != 1) {
          graph_Eadj_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_Eadj_matrix = graph_Eadj_matrix

    message("Edge adjacency matrix ... OK")
  } else {
    message("Edge adjacency matrix ... FAIL")
  }

}

# This function takes as input the Eadjancency matrix and compute, using
# the same algorithm used for the vertexes, the distance matrix for the edges.
# This is done by computing a connection vector that states with which edge each edge is connected
# than it loop through each connection to find the wanted connection (e.g. elemij It loop through
# the connection of 1 until it find the element j, increasing the distance at each loop)

Moleculors$Edistance_matrix = function(){

  if (is.matrix(Mol_mat$graph_Eadj_matrix)) {
    graph_Edistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Eadj_matrix), ncol = nrow(Mol_mat$graph_Eadj_matrix))

    connections = list()
    index_counter = c()

    for (i in 1:nrow(Mol_mat$graph_Eadj_matrix)) {
      connections = append(connections, list(which(Mol_mat$graph_Eadj_matrix[i,] == 1)))
    }


    for (i in 1:nrow(graph_Edistance_matrix)) {
      for (j in 1:nrow(graph_Edistance_matrix)) {
        d = 1
        if (Mol_mat$graph_Eadj_matrix[i,j] == 1) {
          graph_Edistance_matrix[i,j] = 1
        } else if (i == j) {
          graph_Edistance_matrix[i,j] = 0
        } else if (i != j & Mol_mat$graph_Eadj_matrix[i,j] == 0) {
          connection_vector = connections[[i]]
          while (!(j %in% connection_vector)) {
            connections_holder = c()
            for (k in 1:length(connection_vector)) {
              connections_holder = append(connections_holder, connections[[connection_vector[k]]])
            }
            d = d + 1
            connection_vector = connections_holder
          }
          graph_Edistance_matrix[i,j] = d
        }
      }
    }

    Mol_mat$graph_Edistance_matrix = graph_Edistance_matrix

    message("Edge distance matrix ... OK")
  } else {
    message("Edge distance matrix ... FAIL")
  }
}



# This function takes as input the Edistance matrix and compute for each element
# the complement as nrow(Edistance_matrix) - Edistance_matrix[i,j]

Moleculors$ECdistance_matrix = function(){

  if (is.matrix(Mol_mat$graph_Edistance_matrix)) {
    graph_ECdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Edistance_matrix), ncol = nrow(Mol_mat$graph_Edistance_matrix))

    for (i in 1:nrow(graph_ECdistance_matrix)) {
      for (j in 1:nrow(graph_ECdistance_matrix)) {
        if (i != j) {
          graph_ECdistance_matrix[i,j] = nrow(Mol_mat$graph_Edistance_matrix) - Mol_mat$graph_Edistance_matrix[i,j]
        } else {
          graph_ECdistance_matrix[i,j] = 0
        }
      }
    }
    Mol_mat$graph_ECdistance_matrix = graph_ECdistance_matrix

    message("Edge distance complement matrix ... OK")
  } else {
    message("Edge distance complement matrix ... FAIL")
  }
}


# This function compute the vertex version of the so called
# Harary matrix. It take has input the Vdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$V_Harary_matrix = function(){

  if (is.matrix(Mol_mat$graph_Vdistance_matrix)) {
    graph_VHarary_matrix = matrix(nrow = nrow(Mol_mat$graph_Vdistance_matrix), ncol = nrow(Mol_mat$graph_Vdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
        if (i != j) {
          graph_VHarary_matrix[i,j] = 1 /Mol_mat$graph_Vdistance_matrix[i,j]
        } else {
          graph_VHarary_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_VHarary_matrix = graph_VHarary_matrix

    message("Vertex Harary matrix ... OK")
  } else {
    message("Vertex Harary matrix ... FAIL")
  }
}

# This function compute the edge version of the so called
# Harary matrix. It take has input the Edistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$E_Harary_matrix = function(){

  if (is.matrix(Mol_mat$graph_Edistance_matrix)) {
    graph_EHarary_matrix = matrix(nrow = nrow(Mol_mat$graph_Edistance_matrix), ncol = nrow(Mol_mat$graph_Edistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
        if (i != j) {
          graph_EHarary_matrix[i,j] = 1 /Mol_mat$graph_Edistance_matrix[i,j]
        } else {
          graph_EHarary_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_EHarary_matrix = graph_EHarary_matrix

    message("Edge Harary matrix ... OK")
  } else {
    message("Edge Harary matrix ... FAIL")
  }
}

# This function compute the vertex version of the reciprocal
# complement vertex matrix . It take has input the VCdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$Rec_VCdistance = function(){

  if (is.matrix(Mol_mat$graph_VCdistance_matrix)) {
    graph_RecVCdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_VCdistance_matrix), ncol = nrow(Mol_mat$graph_VCdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_VCdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_VCdistance_matrix)) {
        if (i != j) {
          graph_RecVCdistance_matrix[i,j] = 1 /Mol_mat$graph_VCdistance_matrix[i,j]
        } else {
          graph_RecVCdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_RecVCdistance_matrix = graph_RecVCdistance_matrix

    message("Reciprocal vertex complement matrix ... OK")
  } else {
    message("Reciprocal vertex complement matrix ... FAIL")
  }
}




# This function compute the edge version of the reciprocal
# complement edge matrix . It take has input the ECdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the edge distance matrix
#

Moleculors$Rec_ECdistance = function(){

  if (is.matrix(Mol_mat$graph_ECdistance_matrix)) {
    graph_RecECdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_ECdistance_matrix), ncol = nrow(Mol_mat$graph_ECdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_ECdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_ECdistance_matrix)) {
        if (i != j) {
          graph_RecECdistance_matrix[i,j] = 1 /Mol_mat$graph_ECdistance_matrix[i,j]
        } else {
          graph_RecECdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_RecECdistance_matrix = graph_RecECdistance_matrix

    message("Reciprocal edge complement matrix ... OK")
  } else {
    message("Reciprocal edge complement matrix ... FAIL")
  }
}


# This function compute the vertex version of the complementary
# vertex matrix . It take has input the Vdistance_matrix and for each
# element of the matrix if i != j it returns the  value of dmin + dmax - elemij
# were dmin is the minimal lenght in the molecular graph (1 for simple molecules)
# dmax is the maximum lenght in the system and elemij is the value of the element i j
# of the vertex distance matrix
#


Moleculors$Complementary_Vdistance = function(){

  if (is.matrix(Mol_mat$graph_Vdistance_matrix)) {
    graph_VCOMPdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Vdistance_matrix), ncol = nrow(Mol_mat$graph_Vdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
        if (i != j) {
          graph_VCOMPdistance_matrix[i,j] = 1 + max(Mol_mat$graph_Vdistance_matrix) - Mol_mat$graph_Vdistance_matrix[i,j]
        } else {
          graph_VCOMPdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_VCOMPdistance_matrix = graph_VCOMPdistance_matrix

    message("Complementary vertex matrix ... OK")
  } else {
    message("Complementary vertex matrix ... FAIL")
  }
}


# This function compute the edge version of the complementary
# edge matrix . It take has input the Edistance_matrix and for each
# element of the matrix if i != j it returns the  value of dmin + dmax - elemij
# were dmin is the minimal lenght in the molecular graph (1 for simple molecules)
# dmax is the maximum lenght in the system and elemij is the value of the element i j
# of the edge distance matrix
#


Moleculors$Complementary_Edistance = function(){

  if (is.matrix(Mol_mat$graph_Edistance_matrix)) {
    graph_ECOMPdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Edistance_matrix), ncol = nrow(Mol_mat$graph_Edistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
        if (i != j) {
          graph_ECOMPdistance_matrix[i,j] = 1 + max(Mol_mat$graph_Edistance_matrix) - Mol_mat$graph_Edistance_matrix[i,j]
        } else {
          graph_ECOMPdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_ECOMPdistance_matrix = graph_ECOMPdistance_matrix
    message("Complementary edge matrix ... OK")
  } else {
    message("Complementary edge matrix ... FAIL")
  }
}

# This function compute the vertex version of the reciprocal
# complementary vertex matrix . It take has input the Complementary_Vdistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$Rec_Complementary_Vdistance = function(){

  if (is.matrix(Mol_mat$graph_VCOMPdistance_matrix)) {
    graph_VRecCOMPdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_VCOMPdistance_matrix), ncol = nrow(Mol_mat$graph_VCOMPdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_VCOMPdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_VCOMPdistance_matrix)) {
        if (i != j) {
          graph_VRecCOMPdistance_matrix[i,j] = 1 /Mol_mat$graph_VCOMPdistance_matrix[i,j]
        } else {
          graph_VRecCOMPdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_VRecCOMPdistance_matrix = graph_VRecCOMPdistance_matrix
    message("Reciprocal complementary vertex matrix ... OK")
  } else {
    message("Reciprocal complementary vertex matrix ... FAIL")
  }
}

# This function compute the edge version of the reciprocal
# complementary edge matrix . It take has input the Complementary_Edistance_matrix and for each
# element of the matrix if i != j it returns the reciprocal value of
# the distance matrix
#

Moleculors$Rec_Complementary_Edistance = function(){

  if (is.matrix(Mol_mat$graph_ECOMPdistance_matrix)) {
    graph_ERecCOMPdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_ECOMPdistance_matrix), ncol = nrow(Mol_mat$graph_ECOMPdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_ECOMPdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_ECOMPdistance_matrix)) {
        if (i != j) {
          graph_ERecCOMPdistance_matrix[i,j] = 1 /Mol_mat$graph_ECOMPdistance_matrix[i,j]
        } else {
          graph_ERecCOMPdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_ERecCOMPdistance_matrix = graph_ERecCOMPdistance_matrix
    message("Reciprocal complementary edge matrix ... OK")
  } else {
    message("Reciprocal complementary edge matrix ... FAIL")
  }
}

# This function compute the vertex version of the distance-path-matrix
# It take has input the Vdistance_matrix and for each
# element of the matrix if i != j it returns elemij*(elemij + 1)/2 value of
# the distance matrix
#

Moleculors$Distance_path_Vdistance = function(){

  if (is.matrix(Mol_mat$graph_Vdistance_matrix)) {
    graph_Vpathdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Vdistance_matrix), ncol = nrow(Mol_mat$graph_Vdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
        if (i != j) {
          graph_Vpathdistance_matrix[i,j] = Mol_mat$graph_Vdistance_matrix[i,j]*((Mol_mat$graph_Vdistance_matrix[i,j] + 1)/2)
        } else {
          graph_Vpathdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_Vpathdistance_matrix = graph_Vpathdistance_matrix
    message("Distance path vertex matrix ... OK")
  } else {
    message("Distance path vertex matrix ... FAIL")
  }
}

# This function compute the edge version of the distance-path-matrix
# It take has input the Edistance_matrix and for each
# element of the matrix if i != j it returns elemij*(elemij + 1)/2 value of
# the distance matrix
#

Moleculors$Distance_path_Edistance = function(){

  if (is.matrix(Mol_mat$graph_Edistance_matrix)) {
    graph_Epathdistance_matrix = matrix(nrow = nrow(Mol_mat$graph_Edistance_matrix), ncol = nrow(Mol_mat$graph_Edistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
        if (i != j) {
          graph_Epathdistance_matrix[i,j] = Mol_mat$graph_Edistance_matrix[i,j]*((Mol_mat$graph_Edistance_matrix[i,j] + 1)/2)
        } else {
          graph_Epathdistance_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_Epathdistance_matrix = graph_Epathdistance_matrix
    message("Distance path vertex matrix ... OK")
  } else {
    message("Distance path vertex matrix ... FAIL")
  }
}


# This function return the laplacian matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return
# -1 if elemij = 1, sum(elemi, == 1) if i = j (which is the angle of the vertex)
# and 0 for every other element



Moleculors$Laplacian_Vdistance = function(){

  if (is.matrix(Mol_mat$graph_Vdistance_matrix)) {
    graph_Vlaplacian_matrix = matrix(nrow = nrow(Mol_mat$graph_Vdistance_matrix), ncol = nrow(Mol_mat$graph_Vdistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
        if (Mol_mat$graph_Vdistance_matrix[i,j] == 1) {
          graph_Vlaplacian_matrix[i,j] = - 1
        } else if (i == j){
          graph_Vlaplacian_matrix[i,j] = sum(Mol_mat$graph_Vdistance_matrix[i,] == 1)
        } else {
          graph_Vlaplacian_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_Vlaplacian_matrix = graph_Vlaplacian_matrix
    message("Laplacian vertex matrix ... OK")
  } else {
    message("Laplacian vertex matrix ... FAIL")
  }
}


# This function return the laplacian matrix for the edge distance matrix
# it take as input the Edistance_matrix and for each element return
# -1 if elemij = 1, sum(elemi, == 1) if i = j (which is the angle of the edge)
# and 0 for every other element



Moleculors$Laplacian_Edistance = function(){

  if (is.matrix(Mol_mat$graph_Edistance_matrix)) {
    graph_Elaplacian_matrix = matrix(nrow = nrow(Mol_mat$graph_Edistance_matrix), ncol = nrow(Mol_mat$graph_Edistance_matrix))
    for (i in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Edistance_matrix)) {
        if (Mol_mat$graph_Edistance_matrix[i,j] == 1) {
          graph_Elaplacian_matrix[i,j] = - 1
        } else if (i == j){
          graph_Elaplacian_matrix[i,j] = sum(Mol_mat$graph_Edistance_matrix[i,] == 1)
        } else {
          graph_Elaplacian_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_Elaplacian_matrix = graph_Elaplacian_matrix
    message("Laplacian edge matrix ... OK")
  } else {
    message("Laplacian edge matrix ... FAIL")
  }
}

# This function take as input the Vdistance matrix and an interger < then
# nrow(vdistance) and return the adjmatrix for the k(n) neighbour

Moleculors$Vadj_kneigh_matrix = function(n) {

  if (is.matrix(Mol_mat$graph_Vdistance_matrix) & n < nrow(Mol_mat$graph_Vdistance_matrix) & n%%1 == 0){
    graph_Vadj_kneigh_matrix = matrix(nrow = nrow(Mol_mat$graph_Vdistance_matrix), ncol = nrow(Mol_mat$graph_Vdistance_matrix))


    for (i in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
      for (j in 1:nrow(Mol_mat$graph_Vdistance_matrix)) {
        if (i != j & Mol_mat$graph_Vdistance_matrix[i,j] == n){
          graph_Vadj_kneigh_matrix[i,j] = 1
        } else {
          graph_Vadj_kneigh_matrix[i,j] = 0
        }
      }
    }

    Mol_mat$graph_Vadj_kneigh_matrix = graph_Vadj_kneigh_matrix

    message(paste('Vadj_kneigh_matrix computed for n = ', n, sep = ''))

  } else {
    message("no adjacency K neighbour matrix was computed")
  }
}


# This function takes as input the adjacency matrix and an integer n and
# return the n power of the adjacency matrix.

Moleculors$Vadj_highpower_matrix = function(n){

  if (is.matrix(Mol_mat$graph_Vadj_matrix) & n%%1 == 0){

    graph_Vadj_highpower_matrix = is.matrix(Mol_mat$graph_Vdistance_matrix)

    for (i in 2:n) {

      graph_Vadj_highpower_matrix = graph_Vadj_highpower_matrix %*% is.matrix(Mol_mat$graph_Vdistance_matrix)
    }

    Mol_mat$graph_vadj_highpower_matrix = graph_Vadj_highpower_matrix

    message(paste('Vadj_highpower_matrix computed for n = ', n, sep =''))

  } else {
    message("no adj highpower matrix was computed")
  }
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
