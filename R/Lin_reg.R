## This source file contain all the necessary functions to perform a linear regression
# of the dataset.
#

#' Moleculors linear regression
#'
#' Compute a linear model of the descriptors and selected output. Data are automatically splitted
#' in train/validation with default 0.7 ratio.
#'
#' @return Linear regression model of the descriptors and output.
#'
#' @examples
#' Moleculors_lin_reg()
#'
#' @export
#'

Moleculors_lin_reg <- function(descp, output, row_wise = T, split_ratio = 0.7) {
  if ((ncol(descp) != ncol(output) & row_wise == T) | (nrow(descp) != nrow(output) & row_wise == F)) {
    message("Descriptor matrix and output show different size")
    return()
  } else if(split_ratio > 1 | split_ratio < 0) {
    message("Wrong split_ratio definition")
    return()
  } else {
    if (row_wise == F) {
      descp = t(as.matrix(descp))
      output = t(as.matrix(output))
      tdescp = raw_data_prep(descp,output)

    } else {
      tdescp = raw_data_prep(descp,output)
    }

    tdescp_refined_cormat <- dimension_reduction(tdescp)
    tdescp_refined <- tdescp_refined_cormat[[1]]
    model_data <- train_validation_split(tdescp_refined, split_ratio)
    model_coefficient <- train(model_data[[1]])
    validation_datas <- validation(model_data, model_coefficient)




    sink_data()


  }
}

# Prepare data for linear regresion model

raw_data_prep <- function(descriptors, target){
  tdescp = t(as.matrix(descriptors))
  tdescp = as.data.frame(tdescp)
  tdescp$output  = t(as.matrix(target))
  names(tdescp) = c(names(Output_descp), "output")
  return(tdescp)
}

#reduce the number of descriptors by removing those with correlation higher than the set threshold

dimension_reduction <- function(full_matrix, threshold = 0.9){
  cormat = cor(full_matrix[,-ncol(full_matrix)])
  covmat = cov(full_matrix[,-ncol(full_matrix)])
  highcor_index = c()
  for (i in 1:nrow(cormat)) {
    for (j in 1:nrow(cormat)) {
      if (is.na(cormat[i,j])) {
        cormat[i,j] = 0
      }
      if (is.na(covmat[i,j])) {
        covmat[i,j] = 0
      }
    }
  }

  for (i in 1:nrow(cormat)) {
    for (j in 1:nrow(cormat)) {
      if (i != j & !(j %in% highcor_index) & !(i %in% highcor_index) & (cormat[i,j] > threshold | cormat[i,j] < -threshold)) {
        highcor_index = append(highcor_index, j)
      }
    }
  }
  for (i in 1:nrow(covmat)) {
    if (covmat[i,i] > -0.005 & covmat[i,i] < 0.005 & !(i %in% highcor_index)) {
      highcor_index = append(highcor_index, i)
    }
  }
  tdescp_refined = full_matrix[,-highcor_index]
  return(list(tdescp_refined, cormat, covmat))
}

# data_split for model training

train_validation_split <- function(dataset, split_ratio.){
  split_size = round(nrow(dataset)*split_ratio.)
  train_data_index = sample(nrow(dataset), size = split_size)
  train_data = dataset[train_data_index,]
  test_data = dataset[-train_data_index,]
  return(list(train_data, test_data))
}

# model training

train <- function(train_validation_dataset){
  y = train_validation_dataset[,ncol(train_validation_dataset)]
  X = train_validation_dataset[,-ncol(train_validation_dataset)]
  X$b0 = c(rep(1, nrow(train_validation_dataset)))
  b_full = (invMat(t(as.matrix(X)) %*% as.matrix(X))%*%t(as.matrix(X)))%*%y
  rownames(b_full) = names(X)
  return(b_full)
}

# model validation

validation <- function(train_test_list_dataset, model_coeff){
  train_prediction = c(rep(0,nrow(train_test_list_dataset[[1]])))
  soq_temp = c(rep(0,nrow(train_test_list_dataset[[1]])))
  train_test_list_dataset[[1]]$b0 = c(rep(1, nrow(train_test_list_dataset[[1]])))
  train_test_list_dataset[[1]] <- train_test_list_dataset[[1]][, c(1:(ncol(train_test_list_dataset[[1]])-2), ncol(train_test_list_dataset[[1]]), (ncol(train_test_list_dataset[[1]]) -1))]
  for (i in 1:nrow(train_test_list_dataset[[1]])) {
    train_prediction[i] = as.matrix(train_test_list_dataset[[1]][i,-ncol(train_test_list_dataset[[1]])]) %*% as.matrix(model_coeff)
    soq_temp[i] = (train_prediction[i] - as.matrix(train_test_list_dataset[[1]][i,ncol(train_test_list_dataset[[1]])]))^2
    if (i == nrow(train_test_list_dataset[[1]])) {
      train_LOF = sum(soq_temp)/(nrow(train_test_list_dataset[[1]]) - length(model_coeff))
    }
  }
  train_comparison = as.data.frame(as.matrix(train_test_list_dataset[[1]][,ncol(train_test_list_dataset[[1]])]), names(c("reference")))
  train_comparison$train_pred = train_prediction

  test_prediction = c(rep(0,nrow(train_test_list_dataset[[2]])))
  soq_temp = c(rep(0,nrow(train_test_list_dataset[[2]])))
  train_test_list_dataset[[2]]$b0 = c(rep(1, nrow(train_test_list_dataset[[2]])))
  train_test_list_dataset[[2]] <- train_test_list_dataset[[2]][, c(1:(ncol(train_test_list_dataset[[2]])-2), ncol(train_test_list_dataset[[2]]), (ncol(train_test_list_dataset[[2]]) -1))]
  for (i in 1:nrow(train_test_list_dataset[[2]])) {
    test_prediction[i] = as.matrix(train_test_list_dataset[[2]][i,-ncol(train_test_list_dataset[[2]])]) %*% as.matrix(model_coeff)
    soq_temp[i] = (train_prediction[i] - as.matrix(train_test_list_dataset[[2]][i,ncol(train_test_list_dataset[[2]])]))^2
    if (i == nrow(train_test_list_dataset[[2]])) {
      test_LOF = sum(soq_temp)/(nrow(train_test_list_dataset[[2]]) - length(model_coeff))
    }
  }
  test_comparison = as.data.frame(as.matrix(train_test_list_dataset[[2]][,ncol(train_test_list_dataset[[2]])]), names(c("reference")))
  test_comparison$test_pred = test_prediction

  return(list(train_comparison, test_comparison))


}

#sink results to a folder

sink_data <- function(saving_path = getwd()) {
  if (getwd() != paste0(saving_path, "/Moleculors_results")) {
    dir.create(paste0(saving_path, "/Moleculors_results"))
  }
  setwd(paste0(saving_path, "/Moleculors_results"))
  tdescp <- get("tdescp", parent.frame())
  write.csv(tdescp, file = "Complete_descp_and_output_table.csv")
  tdescp_refined_plus_cormat <- get("tdescp_refined_cormat", parent.frame())
  write.csv(tdescp_refined_plus_cormat[[3]], file = "Covariance_matrix.csv")
  write.csv(tdescp_refined_plus_cormat[[2]], file = "Correlation_matrix.csv")
  write.csv(tdescp_refined_plus_cormat[[1]], file = "Refined_descp_output_table.csv")
  train_test_dataset <- get("model_data", parent.frame())
  write.csv(train_test_dataset[[1]], file = "Training_dataset.csv")
  write.csv(train_test_dataset[[2]], file = "Test_dataset.csv")
  model_coeff <- get("model_coefficient", parent.frame())
  write.csv(model_coeff, file = "Model_coefficients.csv")
  validation_results <- get("validation_datas", parent.frame())
  write.csv(validation_results[[1]], file = "Training_comparison.csv")
  write.csv(validation_results[[2]], file = "Test_comparison.csv")


  setwd(saving_path)
}

