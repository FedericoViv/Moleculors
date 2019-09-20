##Moleculors general

Moleculors = new.env()

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

  return("Loading successful")

}

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

  return("Weight... Ok")

}

Moleculors$N_atoms = function(){

  if(is.data.frame(Moleculors$Input)){

    Moleculors$Natoms = nrow(Moleculors$Input)

  } else {

    message("No Input file detected")

    return("N° of atoms... FAILED")

  }

  return("N° of atoms... OK")

}

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

    for (i in (nrow(Input_H_suppressed):1)) {

      Input_H_suppressed$X[i] = Input_H_suppressed$X[i] - Input_H_suppressed$X[1]
      Input_H_suppressed$Y[i] = Input_H_suppressed$Y[i] - Input_H_suppressed$Y[1]
      Input_H_suppressed$Z[i] = Input_H_suppressed$Z[i] - Input_H_suppressed$Z[1]
      print(Input_H_suppressed)
    }

    distance_vector = sqrt(apply(apply(Input_H_suppressed[,-1], 1, `^`,2), 2, sum))

  } else {

    message("No Input file detected")

    return("No graphical matrix was computated")

  }

}
