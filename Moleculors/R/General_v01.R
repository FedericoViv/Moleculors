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

    return(message("No Input file detected"))
  }

  weight_vector = vector()

  for (i in 1:nrow(Moleculors$Input)) {

    weight_vector[i] = Weight_library$Weight[Weight_library$Symbol == as.character(Moleculors$Input$Atom[i])]

  }

  Moleculors$Weight = sum(weight_vector)

  return("Weight .... Ok")

}
