\name{Wavelet-class}
\Rdversion{1.1}
\docType{class}
\alias{Wavelet-class}

\title{Class "Wavelet"}
\description{
A class to hold a generalized wavelet and its common attributes
}
\details{
Objects of primary class Wavelet should never be created. A wavelet has far too many specific requirements for this general class to be practical. This class is only a template for other wavelets that share these common attributes. Objects can be created by calls to the functions \link{Gabor}, \link{Paul}, and \link{DOG}. All of these functions create Wavelet class objects, however, they all inherit from Wavelet, none are primary class Wavelet.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{\link{function}}. The mother wavelet.
A mother wavelet is defined as a function of time f(t) that is absolutely and square
integrable with zero mean and square norm one.}
    \item{\code{normal}:}{Object of class \code{"numeric_or_complex"}. Normalizes the mother wavelet.}
    \item{\code{Fourier Transform}:}{Object of class \code{"FourierTransform"}. The Fourier transform of the mother wavelet.}
    \item{\code{e-folding time}:}{Object of class \code{\link{numeric}}. The decay factor of the mother wavelet.}
    \item{\code{Fourier wavelength}:}{Object of class \code{\link{numeric}} ~~ }
    \item{\code{reconstruction factor}:}{Object of class \code{\link{numeric}} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{function}"}, from data part.

Class \code{"\linkS4class{OptionalFunction}"}, by class "function", distance 2.

Class \code{"\linkS4class{PossibleMethod}"}, by class "function", distance 2.
}
% \section{Methods}{
% \code{\link{}}
% }
\keyword{classes}
