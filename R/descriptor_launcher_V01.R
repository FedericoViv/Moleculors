#' Moleculors descriptor calculator launcher
#'
#' This function runs the Molecular descriptors functions and handle each of
#' them automatically. The user is suggested to start this function rather than each one
#' of the manually if all the descriptors are required.
#'
#' @name Moleculors$descriptor_launcher
#'
#' @usage Moleculors$descriptor_launcher()
#'
#' @return Returns the molecular descriptors '0D, 1D, 2D' for the selected molecule. Results are stored in Output_descp environment.
#'
#' @examples
#' Moleculors$descriptor_launcher()
#'
#' @export

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

