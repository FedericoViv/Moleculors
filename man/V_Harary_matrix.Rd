% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Molecular_graph_V01.R
\name{V_Harary_matrix}
\alias{V_Harary_matrix}
\title{Moleculors vertex harary matrix}
\usage{
V_Harary_matrix()
}
\value{
Vertex Harary matrix for the loaded molecule. Matrix is stored in Mol_mat environment.
}
\description{
This function compute the vertex version of the so called
Harary matrix. It take has input the Vdistance_matrix and for each
element of the matrix if i != j it returns the reciprocal value of
the distance matrix
}
\examples{
V_Harary_matrix()

}
