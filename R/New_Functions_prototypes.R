# This source contain prototypes of functions to be implemented.
# Once implemented will be removed from this file


Hosoya_index_calc = function(Vadj_Input_matrix, Eadj_Input_matrix){

  Z = 1

  for (i in ceiling(1:(nrow(Vadj_Input_matrix)/2))) {

    if (i == 1) {

      Z = Z + nrow(Eadj_Input_matrix)

    } else {

      Z = Z

    }
  }

}





# This function return the Szeged matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return
# -1 if elemij = 1 and 0 if i = j
#



Szeged_Vdistance = function(Vdistance_matrix){

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

  graph_Vszeged_matrix = graph_Vszeged_matrix
}


# This function return the Cluj matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return the
# number of vertex close to i after removing any vertex in the path i-j

cluj_Vdistance = function(Vdistance_matrix){

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

  graph_Vcluj_matrix = graph_Vcluj_matrix
}


# This function return the allpath matrix for the vertex distance matrix
# it take as input the Vdistance_matrix and for each element return the
# number all the possible path to get from i to j

allpath_Vdistance = function(Vdistance_matrix){

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

  graph_Vallpath_matrix = graph_Vallpath_matrix
}
