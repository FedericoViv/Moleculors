% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Molecular_graph_V01.R
\name{Vadj_matrix}
\alias{Vadj_matrix}
\title{Moleculors vertex adjacency matrix}
\usage{
Vadj_matrix(Cart_Input_Hsupp)
}
\arguments{
\item{Cart_Input_Hsupp}{Cartesian coordinates of the molecule without the hydrogen atoms}
}
\value{
Vertex adjacency matrix for the loaded molecule. Matrix is stored in Mol_mat environment.
}
\description{
This function calculate the V adjacency matrix using as Input the Hydrogen suppressed
cartesian matrix. In order to compute the distances it starts by calculating the
magnitude of the vector from each atom to one another. Then after rounding 'too many digits
may lead to bad calculation in the following steps' it choose the smallest value 'different
from 0' has the 1 distance value. A loop is used to assure that every value of lenght 1 has
the same value for the normalization step following later.
a normalization is the computed and  every other value different from 1 is set to 0
}
\examples{
Vadj_matrix(input_H_suppressed)

}
