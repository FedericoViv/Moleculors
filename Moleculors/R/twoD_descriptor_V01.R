## This file contain the two dimentional molecular descriptor
# Each function take as input a graphical matrix calculated in the Molecular_graph source file.
# Each molecular descriptor will be saved the Moleculors environment for future call




# This function take as input the vertex distance graphical matrix and return the Wiener Index
# this index is calculated as half of the summatory of the summatory of the distances in each
# row of the graphical matrix

Moleculors$Wiener_index_calc = function(Vdistance_matrix){

  W = 0

  for (i in 1:nrow(Vdistance_matrix)) {

      W = W + sum(Vdistance_matrix[i,])

  }

  Moleculors$W_index = 0.5 * W
}


# This function take as input the edge adjacency matrix and return the Platt Number
# this number is calculated as the summatory of the summatory of each elementij in each row
# of the adjacency matrix

Moleculors$Platt_number_calc = function(Eadj_Input_matrix){

  Platt = 0

  for (i in 1:nrow(Eadj_Input_matrix)) {

    Platt = Platt + sum(Eadj_Input_matrix[i,])

  }

  Moleculors$Platt_number = Platt
}



# This function take as input the Laplacian Vertex matrix and return
# the Zagreb index for the selected molecular graph as the summatory
# of the squared degree for each atom of the matrix

Moleculors$Zagreb_index_calc = function(VLaplacian_matrix){

  Zagreb = 0

  for (i in 1:nrow(VLaplacian_matrix)) {
    for (j in 1:nrow(VLaplacian_matrix)) {
      if (i == j) {
        Zagreb = Zagreb + (VLaplacian_matrix[i,j] * VLaplacian_matrix[i,j])
      }
    }
  }

  Moleculors$Zagreb_index = Zagreb
}



# This function take as input the Vertex adjacency matrix, the Vertex
# distance matrix,the Edge adjacency matrix and compute the
# Balaban index as Edges/ (Edges - Vertex + 1) * summatory (total_distancei * total_distancej)^-0.5
# for any adjacence vertex. The number (Edges -Vertex + 1) is denominated with u and is
# called cyclomatic number


Moleculors$Balaban_index_calc = function(Vdistance_matrix, Vadj_Input_matrix, Eadj_Input_matrix){

  u = nrow(Eadj_Input_matrix) - nrow(Vadj_Input_matrix) + 1

  J = 0

  for (i in 1:nrow(Vdistance_matrix)) {
    for (j in 1:nrow(Vdistance_matrix)) {
      if (j > i & Vadj_Input_matrix[i,j] == 1) {
        J = J + (sum(Vdistance_matrix[i,]) * sum(Vdistance_matrix[j,]))^-0.5
      }
    }
  }

  J = (nrow(Eadj_Input_matrix)/(u + 1)) * J

  Moleculors$Balaban_index = J

}
