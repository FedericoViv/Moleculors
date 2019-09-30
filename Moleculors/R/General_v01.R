##Moleculors general

# Moleculors envinronment. Every function of output will be part and only part
# of this environment thus preventing global variables overlay.

Moleculors = new.env()


# This function take a csv file containing 4 colums:
# atom symbol, X,Y,Z. Any difference will be detected as
# an error or warning depending on the situation.
#
#
#

Moleculors$Molecular_input = function(){

  cartesian_csv = tryCatch({ cartesian_csv = read.csv(file.choose(),
                                                      header = FALSE)
                      },
                    warning = function(w){
                        warning(w)
                        message("csv file doesn't look properly formatted")
                    },

                    error = function(e){
                        message("Input file doesn't look like a csv file")
                        return(NA)
                    },

                    finally = { message("Always use cartesian coordinates as input!")

                    })

  if (ncol(cartesian_csv) != 4) {

    return(message("Input file has more/less column than expected"))

  }
  names(cartesian_csv) = c("Atom", "X", "Y", "Z")

  print(cartesian_csv)

  Moleculors$Input = cartesian_csv

  return(message("Loading successful"))

}

# This function compute the molecular weight of the input molecule
# by checking the atom symbol in the input matrix/df and crossing it
# with an internal library containing atoms weight.
#
#

Moleculors$molecular_weight = function(){

  if (is.data.frame(Moleculors$Input)) {

    Weight_library = read.csv("tables/weight_table.csv")

  } else {

    message("No Input file detected")

    return(message("Weight... FAILED"))
  }

  weight_vector = vector()

  for (i in 1:nrow(Moleculors$Input)) {

    weight_vector[i] = Weight_library$Weight[Weight_library$Symbol == as.character(Moleculors$Input$Atom[i])]

  }

  Moleculors$Weight = sum(weight_vector)

  return(message("Weight... Ok"))

}

# This function calculate the number of atoms of the
# molecular input by simply calculating the number of rows of the input matrix
#
#
#

Moleculors$N_atoms = function(){

  if(is.data.frame(Moleculors$Input)){

    Moleculors$Natoms = nrow(Moleculors$Input)

  } else {

    message("No Input file detected")

    return("N° of atoms... FAILED")

  }

  return(message("N° of atoms... OK"))

}



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


    for (i in 1:nrow(graph_Vdistance_matrix)) {

      for (j in 1:ncol(graph_Vdistance_matrix)) {

        graph_Vdistance_matrix[i,j] = floor(sqrt((Input_H_suppressed$X[j] - Input_H_suppressed$X[i])^2 +
                                                   (Input_H_suppressed$Y[j] - Input_H_suppressed$Y[i])^2 +
                                                   (Input_H_suppressed$Z[j] - Input_H_suppressed$Z[i])^2))

        if (i == j) {
          graph_Vadj_matrix[i,j] = 0
        } else if (graph_Vdistance_matrix[i,j] == 1){
          graph_Vadj_matrix[i,j] = 1
        } else {
          graph_Vadj_matrix[i,j] = 0
        }
      }
    }
    edge_matrix = matrix(nrow = (nrow(Input_H_suppressed) - 1), ncol = (ncol(Input_H_suppressed)-1))   ##MUST BE FIXED try using the adjacency matrix to check if two atoms are adjacent
    for (i in 1:nrow(edge_matrix)) {
      for (j in 1:ncol(edge_matrix)) {
        simplified_input_H_suppressed = Input_H_suppressed[-i,(j+1)]
        edge_matrix[i,j] = (Input_H_suppressed[i,(j+1)] +
                              simplified_input_H_suppressed[which.min(abs(simplified_input_H_suppressed -
                                                                        Input_H_suppressed[i,(j+1)]))])/2
      }

    }

    graph_Edistance_matrix = matrix(nrow = nrow(edge_matrix), ncol = nrow(edge_matrix))
    graph_Eadj_matrix = matrix(nrow = nrow(edge_matrix), ncol = nrow(edge_matrix))

    for (i in 1:nrow(graph_Edistance_matrix)) {
      for (j in 1:ncol(graph_Edistance_matrix)) {
        graph_Edistance_matrix[i,j] = floor(sqrt((edge_matrix[j,1] - edge_matrix[i,1])^2 +
                                                   (edge_matrix[j,2] - edge_matrix[j,2])^2 +
                                                   (edge_matrix[j,3] - edge_matrix[j,3])^2))

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
    Moleculors$graph_Edistance_matrix = graph_Edistance_matrix
    Moleculors$graph_Eadj_matrix = graph_Eadj_matrix


  } else {

    message("No Input file detected")

    return("No graphical matrix was computed")

  }

  return(message("Computing graphical Matrix... OK"))

}

