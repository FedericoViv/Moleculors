% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Molecular_graph_V01.R
\name{Rec_VCdistance}
\alias{Rec_VCdistance}
\title{Moleculors vertex reciprocal complement matrix}
\usage{
Rec_VCdistance()
}
\value{
Vertex reciprocal complement matrix for the loaded molecule. Matrix is stored in Mol_mat environment.
}
\description{
This function compute the vertex version of the reciprocal
complement vertex matrix . It take has input the VCdistance_matrix and for each
element of the matrix if i != j it returns the reciprocal value of
the distance matrix
}
\examples{
Rec_VCdistance()

}
