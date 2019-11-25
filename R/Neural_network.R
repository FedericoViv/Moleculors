#' Moleculors input file loader.
#'
#' This function take a csv file containing 4 colums:
#' atom symbol, X,Y,Z. Any difference will be detected as
#' an error or warning depending on the situation. NOTE: the input
#' coordinates are supposed to be the optimized coordinates of the
#' target molecule! If this requirement is not fullfilled any descriptor
#' will not be related to the topological property at all.
#'
#'
#' @examples
#' molecular_input()
#'
#'
#' @export

QSAR_NN_model <- function(descriptor, predictor, train_split = 0.75, hl = 4){

  descriptor <- t(descriptor)

  mean_vector <- apply(descriptor, 2, mean)
  sd_vector <- apply(descriptor, 2, sd)

  for (i in 1:nrow(descriptor)) {
    for (j in 1:ncol(descriptor)) {
      descriptor[i,j] <- (descriptor[i,j] - mean_vector[j])/mean_vector[j]
    }
  }

  mean_predictor <- mean(predictor)
  sd_predictor <- sd(predictor)

  for (i in 1:length(predictor)) {
    predictor[i] <- (predictor[i] - mean_predictor)/mean_predictor
  }

  rand_vector <- runif(ncol(descriptor) * nrow(descriptor))

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
    # stores the predicted outcome
    output = matrix(
      rep(0, times = ncol(descriptor)),
      ncol = 1
    )
  )

  tanh <- function(x) {
    2.0 / (1.0 + exp(-2*x)) -1
  }
  tanh_derivative <- function(x) {
       (1.0 - tanh(x)^2)
    }
  # sigmoid <- function(x) {
  #   1.0 / (1.0 + exp(-x))
  # }
  #
  # sigmoid_derivative <- function(x) {
  #   sigmoid(x) * (1.0 - sigmoid(x))
  # }

  loss_function <- function(nn) {
    sum((nn$y - nn$output) ^ 2)
  }

  feedforward <- function(nn) {

    nn$layer1 <- tanh(nn$input %*% nn$weights1)
    nn$output <- tanh(nn$layer1 %*% nn$weights2)

    nn
  }

  backprop <- function(nn) {

    # application of the chain rule to find derivative of the loss function with
    # respect to weights2 and weights1
    d_weights2 <- (
      t(nn$layer1) %*%
        # `2 * (nn$y - nn$output)` is the derivative of the sigmoid loss function
        (2 * (nn$y - nn$output) *
           tanh_derivative(nn$output))
    )

    d_weights1 <- ( 2 * (nn$y - nn$output) * tanh_derivative(nn$output)) %*%
      t(nn$weights2)
    d_weights1 <- d_weights1 * tanh_derivative(nn$layer1)
    d_weights1 <- t(nn$input) %*% d_weights1

    # update the weights using the derivative (slope) of the loss function
    nn$weights1 <- nn$weights1 + d_weights1
    nn$weights2 <- nn$weights2 + d_weights2

    nn
  }

  # number of times to perform feedforward and backpropagation
  n <- 1500

  # data frame to store the results of the loss function.
  # this data frame is used to produce the plot in the
  # next code chunk
  loss_df <- data.frame(
    iteration = 1:n,
    loss = vector("numeric", length = n)
  )

  for (i in seq_len(1500)) {
    my_nn <- feedforward(my_nn)
    my_nn <- backprop(my_nn)

    # store the result of the loss function.  We will plot this later
    loss_df$loss[i] <- loss_function(my_nn)
  }

  # print the predicted outcome next to the actual outcome
  data.frame(
    "Predicted" = round(my_nn$output, 3),
    "Actual" = predictor,
    "Denorm_pred" = round(my_nn$output, 3) * mean_predictor + mean_predictor,
    "Denorm_actual" = predictor * mean_predictor + mean_predictor
  )

}
