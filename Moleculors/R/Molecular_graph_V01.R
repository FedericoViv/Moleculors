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

    graph_Vdistance_matrix = matrix(nrow = nrow(Input_H_suppressed), ncol = nrow(Input_H_suppressed))
    graph_Vadj_matrix = matrix(nrow = nrow(Input_H_suppressed), ncol = nrow(Input_H_suppressed))
    graph_VCdistance_matrix = matrix(nrow = nrow(Input_H_suppressed), ncol = nrow(Input_H_suppressed))


    for (i in 1:nrow(graph_Vdistance_matrix)) {

      for (j in 1:ncol(graph_Vdistance_matrix)) {

        graph_Vdistance_matrix[i,j] = round(sqrt((Input_H_suppressed$X[j] - Input_H_suppressed$X[i])^2 +
                                                   (Input_H_suppressed$Y[j] - Input_H_suppressed$Y[i])^2 +
                                                   (Input_H_suppressed$Z[j] - Input_H_suppressed$Z[i])^2),0)

        graph_VCdistance_matrix[i,j] = nrow(Input_H_suppressed) - round(sqrt((Input_H_suppressed$X[j] - Input_H_suppressed$X[i])^2 +
                                                                               (Input_H_suppressed$Y[j] - Input_H_suppressed$Y[i])^2 +
                                                                               (Input_H_suppressed$Z[j] - Input_H_suppressed$Z[i])^2),0)

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

    edge_matrix = matrix(nrow = nrow(Input_H_suppressed),
                         ncol = (ncol(Input_H_suppressed)-1))
    h = 1
    counter = 1

    for (i in 1:nrow(graph_Vdistance_matrix)) {
      for (j in h:nrow(graph_Vdistance_matrix)) {
        if (graph_Vdistance_matrix[i,j] == 1) {
          edge_matrix[counter, 1] = (Input_H_suppressed$X[i] + Input_H_suppressed$X[j])/2
          edge_matrix[counter, 2] = (Input_H_suppressed$Y[i] + Input_H_suppressed$Y[j])/2
          edge_matrix[counter, 3] = (Input_H_suppressed$Z[i] + Input_H_suppressed$Z[j])/2
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
        } else if (graph_Vdistance_matrix[i,j] == 1){
          graph_Eadj_matrix[i,j] = 1
        } else {
          graph_Eadj_matrix[i,j] = 0
        }
      }
    }

    Moleculors$graph_Vdistance_matrix = graph_Vdistance_matrix
    Moleculors$graph_Vadj_matrix = graph_Vadj_matrix
    Moleculors$graph_VCdistance_matrix = graph_VCdistance_matrix
    Moleculors$graph_Edistance_matrix = graph_Edistance_matrix
    Moleculors$graph_Eadj_matrix = graph_Eadj_matrix


  } else {

    message("No Input file detected")

    return("No graphical matrix was computed")

  }

  return(message("Computing graphical Matrix... OK"))

}



# This function calculate the Vdistance_matrix using as Input the Hydrogen suppressed
# cartesian matrix

Vdistance_matrix = function(Cart_Input_Hsupp){

}
