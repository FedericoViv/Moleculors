#' Moleculors descriptor calculator launcher
#'
#' This function runs the Molecular descriptors functions and handle each of
#' them automatically. The user is suggested to start this function rather than each one
#' of the manually if all the descriptors are required.
#'
#'
#' @return Returns the molecular descriptors '0D, 1D, 2D' for the selected molecule. Results are stored in Output_descp environment.
#'
#' @examples
#' descriptor_launcher()
#'
#' @export

descriptor_launcher = function(){

  if (length(Mol_mat) > 0) {

    molecular_weight()

    N_atoms()

    Wiener_index_calc()

    Platt_number_calc()

    Zagreb_index_calc()

    Balaban_index_calc()

  } else {
    message("No molecular descriptors were computed")
  }

}

