#his function compute the molecular weight of the input molecule
# by checking the atom symbol in the input matrix/df and crossing it
# with an internal library containing atoms weight.
#
#

Moleculors$molecular_weight = function(){

  if (is.data.frame(Mol_mat$input)) {

    Weight_library = read.csv("tables/weight_table.csv")

  } else {

    message("No Input file detected")

    return(message("Weight... FAILED"))
  }

  weight_vector = vector()

  for (i in 1:nrow(Mol_mat$input)) {

    weight_vector[i] = Weight_library$Weight[Weight_library$Symbol == as.character(Mol_mat$input$Atom[i])]

  }

  Output_descp$Weight = sum(weight_vector)

  return(message("Weight... Ok"))

}

# This function calculate the number of atoms of the
# molecular input by simply calculating the number of rows of the input matrix
#
#
#

Moleculors$N_atoms = function(){

  if(is.data.frame(Mol_mat$input)){

    Output_descp$Natoms = nrow(Mol_mat$input)

  } else {

    message("No Input file detected")

    return("N° of atoms... FAILED")

  }

  return(message("N° of atoms... OK"))

}


