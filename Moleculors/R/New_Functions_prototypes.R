# This source contain prototypes of functions to be implemented.
# Once implemented will be removed from this file


Moleculors$Vmatrices_prototype = function(Cart_Input_Hsupp){

  graph_Vdistance_matrix = matrix(nrow = nrow(Cart_Input_Hsupp), ncol = nrow(Cart_Input_Hsupp))
  graph_Vadj_matrix = matrix(nrow = nrow(Cart_Input_Hsupp), ncol = nrow(Cart_Input_Hsupp))
  graph_VCdistance_matrix = matrix(nrow = nrow(Cart_Input_Hsupp), ncol = nrow(Cart_Input_Hsupp))

  for (i in 1:nrow(graph_Vdistance_matrix)) {

    for (j in 1:ncol(graph_Vdistance_matrix)) {

      graph_Vdistance_matrix[i,j] = sqrt((Cart_Input_Hsupp$X[j] - Cart_Input_Hsupp$X[i])^2 +
                                           (Cart_Input_Hsupp$Y[j] - Cart_Input_Hsupp$Y[i])^2 +
                                           (Cart_Input_Hsupp$Z[j] - Cart_Input_Hsupp$Z[i])^2)





    }
  }




}
