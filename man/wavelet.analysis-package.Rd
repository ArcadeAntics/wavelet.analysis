\name{wavelet.analysis-package}
\alias{wavelet.analysis-package}
\alias{wavelet.analysis}
\docType{package}
\title{
\packageTitle{wavelet.analysis}
}
\description{
\packageDescription{wavelet.analysis}
}
\author{
\packageAuthor{wavelet.analysis}

Maintainer: \packageMaintainer{wavelet.analysis}
}
\references{
\href{https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf}{Torrence and Compo 1998}
}

\keyword{ package }
\seealso{
\href{https://en.wikipedia.org/wiki/Morlet_wavelet}{Correction for Gabor (complex Morlet) wavelet}
}
\examples{
# There are three types of wavelets in wavelet.analysis: Gabor (complex Morlet),
# Paul, and DOG (Derivative Of Gaussian)


plot(Gabor(), -3, 3, n = 1001,
    xlab = "time", ylab = "", main = "frequency (\u03c3) = 6")
plot(Paul(), -3, 3, n = 1001,
    xlab = "time", ylab = "", main = "order (m) = 4")
plot(DOG(), -4, 4, n = 1001,
    xlab = "time", ylab = "", main = "derivative (m) = 2")


# Use WaveletTransform to calculate the wavelet transform of a data set


value <- WaveletTransform(rnorm(100), frequency = 32)


print(value)


}
