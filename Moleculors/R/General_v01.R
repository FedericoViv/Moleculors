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



#this function take the cartesian coordinates of a molecule and return the graphical matrix
#of its structure. To do this first of all remove all the hydrogens leaving only the carbons
# and heteroatoms. Then it compute the relative distance for each atom from each other and
# round it up. In the end we will have a matrix n*n where n is the number of atoms
# and each cell contain an integer pointing to the number of bonds intercurring
# between the diagonal element and the others.
#
#


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

    graph_matrix = matrix(nrow = nrow(Input_H_suppressed), ncol = nrow(Input_H_suppressed))

    Input_H_suppressed_norm = Input_H_suppressed

    for (i in 1:nrow(graph_matrix)) {
      for (h in (nrow(Input_H_suppressed):1)) {

        Input_H_suppressed_norm$X[h] = Input_H_suppressed$X[h] - Input_H_suppressed$X[i]
        Input_H_suppressed_norm$Y[h] = Input_H_suppressed$Y[h] - Input_H_suppressed$Y[i]
        Input_H_suppressed_norm$Z[h] = Input_H_suppressed$Z[h] - Input_H_suppressed$Z[i]
      }

      distance_vector = sqrt(apply(apply(Input_H_suppressed_norm[,-1], 1, `^`,2), 2, sum))

      normalizer = min(distance_vector[-i])

      for (t in 1:length(distance_vector)) {

        distance_vector[t] = ceiling(distance_vector[t]/normalizer)

      }

      for (j in 1:ncol(graph_matrix)) {


        graph_matrix[i,j] = abs(distance_vector[j] - distance_vector[i])

        if (i != j & distance_vector[j] - distance_vector[i] == 0) {

          graph_matrix[i,j] = distance_vector[i]
        }

      }
    }

    Moleculors$graph_matrix = graph_matrix


  } else {

    message("No Input file detected")

    return("No graphical matrix was computated")

  }

  return(message("Computing graphical Matrix... OK"))

}
