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

    Randic_index_calc()

    Randic_valence_index_calc()

    First_order_kappa_calc()

    E_state()

    Bulk_index_calc()

    Hydro_factor_calc()

    MED_index_calc()

    Data_summary <- data.frame(1)

    for (i in 1:length(names(Output_descp))) {

      Data_summary[i] <- get(names(Output_descp)[i], envir = Output_descp)

      names(Data_summary)[i] <- names(Output_descp)[i]
    }

    assign("Data_summary", Data_summary, envir = globalenv())

  } else {
    message("No molecular descriptors were computed")
  }

}

