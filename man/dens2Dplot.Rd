\name{dens2Dplot}
\alias{dens2Dplot}
\title{Two dimensional density plot}
\description{
  This function generates a scatter plot enhancing density with different colors.
}
\usage{
dens2Dplot(x, y, nbins=100, pal=NULL, xlab='', ylab='', xlim=c(NA,NA),
ylim=c(NA,NA), ...)
}
\arguments{
  \item{x, y}{Variables for the scatter plot.}
  \item{nbins}{Number of bins in X and Y coordinates of the scatter plot.}
  \item{pal}{Color palette}
  \item{xlab, ylab}{Labels for X and Y axis}
  \item{xlim, ylim}{Limits for X and Y axis}
  \item{\dots}{Other graphic parameters}
}
\details{
  This function has been inspired by Elzizi's answer at http://stackoverflow.com/questions/18089752/r-generate-2d-histogram-from-raw-data
  The user can add a grid, title and other details to the scatter plot.
}
\examples{
n=1000000; x=rnorm(n,15,4); y=x+rnorm(n,5,6)
dens2Dplot(x,y,xlab='Variable X',ylab='Variable Y',las=1)

## Let's add a grid and a title:
grid(col=grey(.4))
title('Example of dens2Dplot')
}
\keyword{hplot}
