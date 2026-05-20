\name{geodesic_pca_sphere}
\alias{geodesic_pca_sphere}
\title{
BIC and ICL to choose the number of components in a model based clustering with elliptically symmetric distributions
}

\description{
BIC and ICL to choose the number of components in a model based clustering with elliptically symmetric distributions
}

\usage{
geodesic_pca_sphere(x, k = 2)
}

\arguments{
\item{x}{
A matrix containing directional data.
}
\item{k}{
The maximum number of clusters to be tested. Default value is 5.
}
}

\details{
The function computes the BIC and ICL to decide on the optimal number of clusters when
using mixtures of SESPC or mixtures of ESAG distributions.
}

\value{
A list including:
\item{x_reduced}{
The data projected onto the sphere (or a lower dimension hyper-sphere).
}
\item{eigenvalues}{
The BIC values for all the models tested.
}
\item{means}{
The BIC values for all the models tested.
}
\item{tangent_coords}{
The ICL values for all the models tested.
}
}

\references{
Biernacki C., Celeux G. and Govaert, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. IEEE Transactions on Pattern Analysis and Machine Intelligence, 22(7): 719--725.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{frechet_mean_sphere}, \link{log_map_sphere}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
x <- x / sqrt( rowSums(x^2) )
a <- geodesic_pca_sphere(x)
}

