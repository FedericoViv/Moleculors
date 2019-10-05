# This source contain prototypes of functions to be implemented.
# Once implemented will be removed from this file


# Vadj prototypes

####BUONA####
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
      if ((graph_Vadj_matrix[i,j] - min(graph_Vadj_matrix[1,-1])) <= 0.05 & i != j) {
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

  Moleculors$graph_Vadj_matrix = graph_Vadj_matrix
}

#prototype of VCdistances
####BUONA####
Moleculors$VCdistance_matrix = function(Vdistance_matrix){

  graph_VCdistance_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))

  for (i in 1:nrow(graph_VCdistance_matrix)) {
    for (j in 1:nrow(graph_VCdistance_matrix)) {
      if (i != j) {
        graph_VCdistance_matrix[i,j] = nrow(Vdistance_matrix) - Vdistance_matrix[i,j]
      } else {
        graph_VCdistance_matrix[i,j] = 0
      }
    }
  }
  Moleculors$VCdistance_matrix = graph_VCdistance_matrix
}


## prototype of vdistances
######BUONA, need testing#####

# This function take the adjacency matrix as input and generate a connection list where the index
# is related to the atom and the value in the index list are the atoms to which that atom is connected
# For each element with d != 1 is then calculated the distance by taking the value in the connection list
# copying those into the connection vector and if the interested atom (j) is not in the connection vector
# a new vector is created with all the connection of the atom to which the first atoms was connected
# this is looped increasing distance at each failed loop to detect the interested element.

Moleculors$Vdistance_prototype = function(Vadj_Input_matrix){

  graph_Vdistance_matrix = matrix(nrow = nrow(Vadj_Input_matrix), ncol = nrow(Vadj_Input_matrix))

  connections = list()
  index_counter = c()

  for (i in 1:nrow(Vadj_Input_matrix)) {
    connections = append(connections, list(which(Vadj_Input_matrix[i,] == 1)))
  }


  for (i in 1:nrow(graph_Vdistance_matrix)) {
    for (j in 1:nrow(graph_Vdistance_matrix)) {
      d = 1
      if (graph_Vadj_matrix[i,j] == 1) {
        graph_Vdistance_matrix[i,j] = 1
      } else if (i == j) {
        graph_Vdistance_matrix[i,j] = 0
      } else if (i != j & graph_Vadj_matrix[i,j] == 0) {
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

  Moleculors$Vdistance_matrix = graph_Vdistance_matrix
}






# This function return the Szeged matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return
# -1 if elemij = 1 and 0 if i = j
#



Moleculors$Szeged_Vdistance = function(Vdistance_matrix){

  graph_Vszeged_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))
  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (Vdistance_matrix[i,j] == 1) {
        graph_Vszeged_matrix[i,j] = - 1
      } else if (i == j){
        graph_Vszeged_matrix[i,j] = sum(Vdistance_matrix[i,] == 1)
      } else {
        graph_Vszeged_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Vszeged_matrix = graph_Vszeged_matrix
}


# This function return the Cluj matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return the
# number of vertex close to i after removing any vertex in the path i-j

Moleculors$cluj_Vdistance = function(Vdistance_matrix){

  graph_Vcluj_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))
  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (Vdistance_matrix[i,j] == 1) {
        graph_Vcluj_matrix[i,j] = - 1
      } else if (i == j){
        graph_Vcluj_matrix[i,j] = sum(Vdistance_matrix[i,] == 1)
      } else {
        graph_Vcluj_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Vcluj_matrix = graph_Vcluj_matrix
}


# This function return the allpath matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return the
# number all the possible path to get from i to j

Moleculors$allpath_Vdistance = function(Vdistance_matrix){

  graph_Vallpath_matrix = matrix(nrow = nrow(Vdistance_matrix), ncol = nrow(Vdistance_matrix))
  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (Vdistance_matrix[i,j] == 1) {
        graph_Vallpath_matrix[i,j] = - 1
      } else if (i == j){
        graph_Vallpath_matrix[i,j] = sum(Vdistance_matrix[i,] == 1)
      } else {
        graph_Vallpath_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Vallpath_matrix = graph_Vallpath_matrix
}
