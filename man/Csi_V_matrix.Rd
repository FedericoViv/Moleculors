% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Molecular_graph_V01.R
\name{Csi_V_matrix}
\alias{Csi_V_matrix}
\title{Moleculors csi vertex matrix}
\usage{
Csi_V_matrix()
}
\value{
Csi vertex matrix for the loaded molecule. Matrix is stored in Mol_mat environment.
}
\description{
This function return the csi vertex matrix for the selected molecule
it take as input the Vadj_matrix, the laplacian_Vertex matrix and the Vdistance matrix
 and for each element return wij = "di x dk... x dj"^-1/2 where di is the vertex degree
 of the element i in the shorter path di-dj. By definition wij with i=j will be 0
}
\examples{
Csi_V_matrix()

}
