#' Moleculors Neural network.
#'
#' This function takes as input the descriptor matrix, the predicted values
#' and return the performance in prediction of a neural network with one hidden layer
#' the train_split argument will determine how data are divided for train/test/validation.
#'
#' @param train_split train test split ratio for the input data
#' @param hl number of nodes in the hidden layer
#'
#' @examples
#' QSAR_NN_model()
#'
#'
#' @export

QSAR_NN_model <- function(descriptor, predictor, train_split = 0.75, hl = 4, lr = 0.1){

  descriptor <- t(descriptor)

  mean_vector <- apply(descriptor, 2, mean)

  for (i in 1:nrow(descriptor)) {
    for (j in 1:ncol(descriptor)) {
      descriptor[i,j] <- (descriptor[i,j] - mean_vector[j])/mean_vector[j]
    }
  }

  mean_predictor <- mean(predictor)

  for (i in 1:length(predictor)) {
    predictor[i] <- (predictor[i] - mean_predictor)/mean_predictor
  }


  rand_vector <- runif(ncol(descriptor) * hl)

  rand_matrix <- matrix(
    rand_vector,
    nrow = ncol(descriptor),
    ncol = nrow(descriptor),
    byrow = TRUE
  )

  my_nn <- list(
    input = descriptor,

    weights1 = rand_matrix,

    weights2 = matrix(runif(hl), ncol = 1),

    y = predictor,

    output = matrix(
      rep(0, times = ncol(descriptor)),
      ncol = 1
    )
  )

  tanh <- function(x) {
    (1.0 - exp(-2*x))  / (1.0 + exp(-2*x))
  }
  tanh_derivative <- function(x) {
    (1.0 - x^2)
  }

  loss_function <- function(nn) {
    sum((nn$y - nn$output) ^ 2)
  }

  feedforward <- function(nn) {

    nn$layer1 <- tanh(nn$input %*% nn$weights1)
    nn$output <- tanh(nn$layer1 %*% nn$weights2)

    nn
  }

  backprop <- function(nn) {

    d_weights2 <- (
      t(nn$layer1) %*%

        (2 * (nn$y - nn$output) *
           tanh_derivative(nn$output))
    )

    d_weights1 <- ( 2 * (nn$y - nn$output) * tanh_derivative(nn$output)) %*%
      t(nn$weights2)
    d_weights1 <- d_weights1 * tanh_derivative(nn$layer1)
    d_weights1 <- t(nn$input) %*% d_weights1

    nn$weights1 <- nn$weights1 + lr*d_weights1
    nn$weights2 <- nn$weights2 + lr*d_weights2

    nn
  }

  n <- 1500

  loss_df <- data.frame(
    iteration = 1:n,
    loss = vector("numeric", length = n)
  )

  for (i in seq_len(1500)) {
    my_nn <- feedforward(my_nn)
    my_nn <- backprop(my_nn)

    loss_df$loss[i] <- loss_function(my_nn)
  }

  nn_stats <-data.frame(
    "Predicted" = round(my_nn$output, 3),
    "Actual" = predictor,
    "Denorm_pred" = round(my_nn$output, 3) * mean_predictor + mean_predictor,
    "Denorm_actual" = predictor * mean_predictor + mean_predictor
  )
  assign("nn_stats", nn_stats, envir = globalenv())
}


## prototype version ##



library(ggplot2)

ggplot(data = loss_df, aes(x = iteration, y = loss)) +
  geom_line()
