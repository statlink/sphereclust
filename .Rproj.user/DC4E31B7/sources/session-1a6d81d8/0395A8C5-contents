\name{Interactive 3D plot of spherical data}
\alias{sphereplot}
\title{
Interactive 3D plot of spherical data
}

\description{
Interactive 3D plot of spherical data.
}

\usage{
sphereplot(dat, col = NULL, bgcol = "snow")
}

\arguments{
\item{dat}{
A matrix with three columns, unit-vectors, spherical data.
}
\item{col}{
If you want the points to appear with different colours put numbers here, otherwise leave it NULL.
}
\item{bgcol}{
The color of the surface of the sphere.
}
}

%\details{
%}

\value{
An interactive 3D plot of the spherical data will appear.
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{mixsespc.contour}
}
}

\examples{
x <- as.matrix( iris[, 1:3] )
x <- x / sqrt( rowSums(x^2) )
sphereplot(x, col = as.numeric(iris[, 5]))
}

