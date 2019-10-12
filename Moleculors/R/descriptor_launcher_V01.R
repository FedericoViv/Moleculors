# This source file contain the Molecular descriptor functions caller and handler
# If this package is intended to be used in its basic application the user is suggested
# to start directly the function in this file to compute every available descriptor simultaneously


Moleculors$descriptor_launcher = function(){

  if (length(Mol_mat) > 0) {

    Moleculors$molecular_weight()

    Moleculors$N_atoms()

    Moleculors$Wiener_index_calc()

    Moleculors$Platt_number_calc()

    Moleculors$Zagreb_index_calc()

    Moleculors$Balaban_index_calc()

  } else {
    message("No molecular descriptors were computed")
  }

}

