\name{Simulation of random values from a mixture of elliptically symmetric distributions}
\alias{rmixsespc}
\alias{rmixesag}
\title{
Simulation of random values from a mixture of elliptically symmetric distributions
}

\description{
The function simulates random values simulation from a given mixture of elliptically symmetric distributions.
}

\usage{
rmixsespc(n, probs, mu, theta)
rmixesag(n, probs, mu, gam)
}

\arguments{
\item{n}{
The sample size.
}
\item{probs}{
This is avector with the mixing probability of each group.
}
\item{mu}{
A matrix with the mean direction of each group.
}
\item{theta}{
A matrix with the \eqn{\theta} parameters of the SESPC distribution of each group.
}
\item{gam}{
A matrix with the \eqn{\gamma} parameters of the ESAG distribution of each group.
}
}

\details{
The function simulates random values simulation from a given mixture of SESPC
or a mixture of ESAG distributions.
}

\value{
A list including:
\item{id}{
An indicator of the group of each simulated vector.
}
\item{x}{
A matrix with the simulated data.
}
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Theodoros Perdikis
\email{theoperdikis92@gmail.com}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mixsespc.mle}, \link{bic.mixsespc}
}
}

\examples{
probs <- c(0.5, 0.5)
mu <- matrix(rnorm(6, 10, 5), ncol = 3)
theta <- matrix( c(-2, 1, -2, 1), ncol = 2 )
x <- rmixsespc(200, probs, mu, theta)$x
}
