#This source file contains all the necessary algebraic function in order to let the moleculors package works



#' Moleculors Ones Matrix
#'
#' Generates a ixj matrix filled with ones
#'
#' @return Generates a ixj matrix filled with ones.
#'
#' @examples
#' OneMat(2)
#'
#' @export
#'

OneMat = function(i,j=i) {
  Y = matrix(data = c(rep(1, i*j)), nrow = i, ncol = j)
  return(Y)
}

#' Moleculors Zeros Matrix
#'
#' Generates a ixj matrix filled with zeros
#'
#' @return Generates a ixj matrix filled with zero.
#'
#' @examples
#' ZeroMat(2)
#'
#' @export
#'

ZeroMat = function(i,j=i){
  Y = matrix(data = c(rep(0, i*j)), nrow = i, ncol = j)
  return(Y)
}


#' Moleculors inverse matrix
#'
#' Compute the inverse matrix of the input A using the determinant and the adj matrix
#' this function is not computing efficient due to the computation of the determinant
#' is supposed to be used for small matrices
#'
#' @return Inverse matrix of the input matrix A
#'
#' @examples
#' A = matrix(1,2,3,5)
#' invMat_old_not_supported(A)
#'
#' @export
#'

invMat_old_not_supported = function(A) {
  if (is.matrix(A) == FALSE) {
    return("Pointed object is not a matrix")
  } else if (det(A) == 0) {
    return("Inversion can't be performed on a singular matrix")
  }
  INV = ZeroMat(nrow(A), ncol(A))

  for (i1 in 1:nrow(A)) {
    for (j1 in 1:ncol(A)) {
      B = A[-i1,-j1]
      if (is.matrix(B) == FALSE) {
        INV[j1,i1] = (-1)^(i1+j1)*B
      } else {
        INV[j1,i1] = (-1)^(i1+j1)*det(B)
      }
    }
  }
  INV = (1/det(A))*INV
  return(INV)
}

#' Moleculors inverse matrix with gauss jordan algorithm
#'
#' Compute the inverse matrix of the input A using the gauss jordan algorithm
#' Given A*X-1 = I we solve with gauss algorithm the augmented matrix 'AI'
#' by forward elimination and backward we can manage to obtain 'IX-1'
#'
#' @return Inverse matrix of the input matrix A
#'
#' @examples
#' A = matrix(1,2,3,5)
#' invMat(A)
#'
#' @export
#'
#'

invMat <- function(A){
  if (is.matrix(A) == FALSE) {
    return("Pointed object is not a matrix")
  } else if (nrow(A) != ncol(A)) {
    return("Pointed object is not a square matrix")
  }
  I = diag(nrow(A))

  for (i in 1:nrow(A)) {
    if (A[i,i] != max(A[i:nrow(A),i])) {
      if (i == nrow(A)) {
        return("Pointed object is not a non-singular matrix")
      }
      for (k in (i+1):nrow(A)) {
        if (A[k,i] == max(A[k:nrow(A),i])) {
          temp_row = A[k,]
          A[k,] = A[i,]
          A[i,] = temp_row
          temp_row = I[k,]
          I[k,] = I[i,]
          I[i,] = temp_row
          break
        }
      }
    }
    I[i,] = I[i,]/A[i,i]
    A[i,] = A[i,]/A[i,i]


    for (j in 1:nrow(A)) {
      if (i != j & A[j,i] != 0) {
        ratio = A[j,i]
        for (h in 1:ncol(A)) {
          A[j,h] = A[j,h] -ratio*A[i,h]
          I[j,h] = I[j,h] -ratio*I[i,h]
        }
      }
    }
  }
  return(I)
}



#' Moleculors Hilbert matrix
#'
#' Compute a nxn hilbert matrix. It can be used to check the quality of the
#' solution algorithm in a quasi singular matrix. How to use: generate a nxn hilbert matrix
#' compute the b coefficient directly by multiplying the hylbert matrix by a vector x with n random
#' values. Solve the equation Hx=b with the preferred solution algorithm. The better xcalc is closed to
#' the x generated, the better the algorightm.
#'
#' @return Hilbert matrix, randomized x and computed b
#'
#' @examples
#' H = Hilb(10)
#'
#'
#' @export
#'
#'

Hilb <- function(n){
  H <- ZeroMat(n)
  for (i in 1:n) {
    for (j in 1:n) {
      H[i,j] = 1/(i+j-1)
    }
  }
  x <- runif(n)
  b <- H%*%x
  output <- list(H,x,b)
  names(output) <- c("H","x","b")
  return(output)
}

#' Moleculors Orthonormal Matrix
#'
#' Compute the orthonormal matrix of the input matrix A supposed non singular
#' and composed of linearly independent column. Gram-Schimdt algorithm
#' is used to compute the QR factorization of the input matrix A. Where
#' Q is the orthonormal matrix and R is a triangular matrix computed
#' following Gram-Schimdt algorigthm. The orthonormal bassis are computed
#' from the matrix column in the following way.
#' A1 = r11u1
#' A2 = r12u1 + r22u2
#' A3 = r13u1 + r23u2 + r33u3
#'
#' where r11 is the dot product of the column 1 of the matrix A,
#' r12 is computed as the dot product of u1 and A2 and rjj is computed
#' as sqrt"Aj^2 - r1,j^2 ... r1,j-1,j^2"
#' uj is then easily computed as uj = Aj - r1ju1 .... rj-1,juj-1/rjj
#'
#'
#' @return The orthonormal basis set of the input matrix A
#'
#' @examples
#' Obasis = Orthonormal(A)
#'
#'
#' @export
#'
#'

Orthonormal <- function(A){
  R = A
  for (i in 1:nrow(A)) {
    Rii = A[i,]%*%A[,i]
    A[,i] = A[,i]/as.numeric(sqrt(Rii))
    for (j in i:ncol(A)) {
      A[,j] = A[,j] - A[,i]
    }
  }

  return(A)
}
