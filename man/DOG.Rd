\name{DOG}
\alias{DOG}
\title{Create a Derivative Of Gaussian Wavelet}
\description{
A DOG wavelet is defined as
\deqn{\Psi(t) = D^m(e^{\frac{-t^2}{2}})}{
      \Psi(t) = D^m(exp(-t^2/2))}
the m-th derivative of a standard Gaussian function (a Gaussian function with mean
0 and variance 1)
}
\usage{
DOG(m = 2L)
}
\arguments{
  \item{m}{
  positive integer; the order / / derivative number of the DOG wavelet
}
}
\details{
A DOG wavelet is a real-valued wavelet, and as such is better for isolating
peaks or discontinuities in a signal

A DOG wavelet of order 2 is also called the Mexican hat / / Ricker wavelet.
}
\value{
A function of class c("DOG", "Wavelet", "function")
}
\references{
\href{https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf}{Torrence and Compo 1998}
}
\examples{
plot(DOG(), -5, 5, n = 1001,
    xlab = "t", ylab = "\u03a8(t)", main = "Mexican hat wavelet")


plot(DOG(6), -5, 5, n = 1001,
    xlab = "t", ylab = "\u03a8(t)", main = "DOG wavelet, order 6")
}
