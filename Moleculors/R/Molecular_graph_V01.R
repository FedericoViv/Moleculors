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


    Moleculors$Vmatrices_copy(Input_H_suppressed)

    Moleculors$Ematrices(Input_H_suppressed, Moleculors$graph_Vadj_matrix)


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
  graph_Eadj_matrix = matrix(nrow = nrow(edge_matrix), ncol = nrow(edge_matrix))

  for (i in 1:nrow(graph_Edistance_matrix)) {
    for (j in 1:ncol(graph_Edistance_matrix)) {
      graph_Edistance_matrix[i,j] = floor(sqrt((edge_matrix[j,1] - edge_matrix[i,1])^2 +
                                                 (edge_matrix[j,2] - edge_matrix[i,2])^2 +
                                                 (edge_matrix[j,3] - edge_matrix[i,3])^2))

      if (i == j) {
        graph_Eadj_matrix[i,j] = 0
      } else if (Vadjmatrix[i,j] == 1){
        graph_Eadj_matrix[i,j] = 1
      } else {
        graph_Eadj_matrix[i,j] = 0
      }
    }
  }

  Moleculors$graph_Edistance_matrix = graph_Edistance_matrix
  Moleculors$graph_Eadj_matrix = graph_Eadj_matrix
}


