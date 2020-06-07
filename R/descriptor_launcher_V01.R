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

    Bulk_Electronegativity_indexes_calc()

    Hydro_factor_calc()

    eigenvalues_descp_cacl()

    Q_polarity_calc()


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


#' Moleculors multiple graph and descriptor calculator.
#'
#' This function loop through the multiple input selected and calculate
#' the graphical matrices and descriptor for each molecule. Results are stored
#' in the descriptor matrix where each column is a molecule and each row are the
#' orded descriptors
#'
#'
#'
#' @examples
#' moleculors_multiple_descriptor()
#'
#'
#' @export


moleculors_multiple_descriptor <- function() {

  if (length(Mol_mat$input_list) > 0 ) {

    descriptor_matrix = matrix(ncol = length(Mol_mat$input_list), nrow = 35)


    for (i in 1:length(Mol_mat$input_list)) {

      Mol_mat$input <- Mol_mat$input_list[[i]]

      molecular_weight()

      N_atoms()

      graphical_matrix()

      descriptor_launcher()

      descriptor_matrix[,i] <- t(as.matrix(Data_summary))

    }

    assign("descriptor_matrix", descriptor_matrix, envir = globalenv())

  }
}

