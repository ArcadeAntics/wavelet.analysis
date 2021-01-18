

NormalAttrSymbol <- "normal"
FourierTransformAttrSymbol <- "fourierTransform"
EFoldingTimeAttrSymbol <- "eFoldingTime"
FourierWavelengthAttrSymbol <- "fourierWavelength"
ReconstructionFactorAttrSymbol <- "reconstructionFactor"


FourierTransformClassSymbol <- "FourierTransform"
WaveletClassSymbol <- "Wavelet"
DOGWaveletClassSymbol <- "DOG"
GaborWaveletClassSymbol <- "Gabor"
PaulWaveletClassSymbol <- "Paul"










setClass(
    Class     = FourierTransformClassSymbol,
    contains  = c("function", "oldClass"),
    slots     = `names<-`("numeric_or_complex", NormalAttrSymbol),
    prototype = prototype(.S3Class = c(FourierTransformClassSymbol, "function")))


setValidity(FourierTransformClassSymbol, function (object)
{
    if (!identical(formals(object), pairlist(omega = quote(expr = ))))
        "object must be a function of a single parameter \"omega\""
    else if (!is.call(body(object)))
        "object body must be class \"call\""
    else if (!tryCatch(is.numeric_or_complex(object(0)), error = FALSE))
        "object must produce numeric / / complex output"
    else validNormal(Normal(object))
})


# setOldClass(c(FourierTransformClassSymbol, "function"), S4Class = FourierTransformClassSymbol)










setClass(
    Class     = WaveletClassSymbol,
    contains  = c("function", "oldClass"),
    slots     = `names<-`(
        c("numeric_or_complex", FourierTransformClassSymbol, "numeric"             , "numeric"                  , "numeric"                     ),
        c(NormalAttrSymbol    , FourierTransformAttrSymbol , EFoldingTimeAttrSymbol, FourierWavelengthAttrSymbol, ReconstructionFactorAttrSymbol)),
    prototype = prototype(.S3Class = c(WaveletClassSymbol, "function")))


setValidity(WaveletClassSymbol, function (object)
{
    if (!identical(formals(object), pairlist(t = quote(expr = ))))
        return("object must be a function of a single parameter \"t\"")
    if (!is.call(body(object)))
        return("object body must be class \"call\"")
    if (!tryCatch(is.numeric_or_complex(object(0)), error = FALSE))
        return("object must produce numeric / / complex output")

    x <- validNormal(Normal(object))
    if (!is.null(x))
        return(x)

    validObject(FourierTransform(object))

    x <- validEFoldingTime(EFoldingTime(object))
    if (!is.null(x))
        return(x)

    x <- validFourierWavelength(FourierWavelength(object))
    if (!is.null(x))
        return(x)

    validReconstructionFactor(ReconstructionFactor(object))
})


# setOldClass(c(WaveletClassSymbol, "function"), S4Class = WaveletClassSymbol)










setClass(
    Class     = DOGWaveletClassSymbol,
    contains  = WaveletClassSymbol,
    slots     = c(derivative = "integer"),
    prototype = prototype(.S3Class = c(DOGWaveletClassSymbol, WaveletClassSymbol, "function")))


setValidity(DOGWaveletClassSymbol, function (object)
{
    m <- attr(object, "derivative")
    if (length(m) != 1L)
        paste0("slot \"derivative\" must be of length 1")
    else if (is.na(m) || m <= 0L)
        "slot \"derivative\" must be finite and positive"
    else if (!tryCatch(gamma(m + 1/2), warning = function(w) FALSE))
        "slot \"derivative\" out of range in \"gammafn\""
})


# setOldClass(c(DOGWaveletClassSymbol, "function"), S4Class = DOGWaveletClassSymbol)










setClass(
    Class     = GaborWaveletClassSymbol,
    contains  = WaveletClassSymbol,
    slots     = c(frequency = "numeric"),
    prototype = prototype(.S3Class = c(GaborWaveletClassSymbol, WaveletClassSymbol, "function")))


setValidity(GaborWaveletClassSymbol, function (object)
{
    sigma <- attr(object, "frequency")
    if (length(sigma) != 1L)
        paste0("slot \"frequency\" must be of length 1")
    else if (!is.finite(sigma) || sigma == 0)
        "slot \"frequency\" must be finite and non-zero"
})


# setOldClass(c(GaborWaveletClassSymbol, "function"), S4Class = GaborWaveletClassSymbol)










setClass(
    Class     = PaulWaveletClassSymbol,
    contains  = WaveletClassSymbol,
    slots     = c(order = "integer"),
    prototype = prototype(.S3Class = c(PaulWaveletClassSymbol, WaveletClassSymbol, "function")))


setValidity(PaulWaveletClassSymbol, function (object)
{
    m <- attr(object, "order")
    if (length(m) != 1L)
        paste0("slot \"order\" must be of length 1")
    else if (is.na(m) || m <= 0L)
        "slot \"order\" must be finite and positive"
    else if (!tryCatch(factorial(2 * m), warning = function(w) FALSE))
        "slot \"order\" out of range in \"gammafn\""
})


# setOldClass(c(PaulWaveletClassSymbol, "function"), S4Class = PaulWaveletClassSymbol)










setClass(
    Class     = "WaveletTransform",
    contains  = c("matrix", "oldClass"),
    slots     = c(series = "numeric_or_complex", wavelet = WaveletClassSymbol,
        frequency = "numeric", deltaj = "numeric"),
    prototype = prototype(.S3Class = c("WaveletTransform", "matrix")))


setValidity("WaveletTransform", function (object)
{
    if (!is.numeric_or_complex(object))
        return("object must be a numeric or complex matrix")
    N <- nrow(object)
    if (N < 3L)
        return("object must have at least 3 rows")
    n <- length(series(object))
    if (n != N)
        return(paste0("length of slot \"series\" and number of rows must be equal (",
            n, ", ", N, ")"))
    validObject(wavelet(object))
    frequency <- frequency(object)
    if (length(frequency) != 1L)
        return("slot \"frequency\" must be of length 1")
    else if (!is.finite(frequency) || frequency <= 0)
        return("slot \"frequency\" must be finite and positive")
    deltaj <- deltaj(object)
    if (length(deltaj) != 1L)
        return("slot \"deltaj\" must be of length 1")
    else if (!is.finite(deltaj) || deltaj <= 0)
        return("slot \"deltaj\" must be finite and positive")
    j <- J(N = N, deltaj = deltaj)
    nc <- ncol(object)
    if (nc != j + 1)
        return(sprintf("object should have %d column%s, has %d column%s",
            j + 1, if (j) "s" else "", nc, if (nc != 1L) "s" else ""))
})


# setOldClass(c("WaveletTransform", "matrix"), S4Class = "WaveletTransform")










setClass(
    Class     = "CompatibleWaveletTransforms",
    contains  = c("list", "oldClass"),
    prototype = prototype(.S3Class = c("CompatibleWaveletTransforms", "list")))


setValidity("CompatibleWaveletTransforms", function (object)
{
    if (!length(object))
        return(TRUE)
    if (!all(vapply(object, is.WaveletTransform, NA)))
        stop("each element of 'object' should be or extend class \"WaveletTransform\"")
    lapply(object, validObject)
    nrows <- unique(vapply(object, nrow, 0L))
    if (length(nrows) != 1L)
        return(paste0("differing number of rows: ", paste0(nrows, collapse = ", ")))
    ncols <- unique(vapply(object, ncol, 0L))
    if (length(ncols) != 1L)
        return(paste0("differing number of columns: ", paste0(ncols, collapse = ", ")))
    frequencies <- unique(vapply(object, frequency, 0))
    if (length(frequencies) != 1L)
        return(paste0("differing 'deltat' values: ", paste0(signif(1/frequencies, getOption("digits")), collapse = ", ")))
    wavelets <- unique(lapply(lapply(object, wavelet), deparse))
    if (length(wavelets) != 1L)
        return("differing wavelets when deparsed")
    envirs <- lapply(lapply(object, wavelet), environment)
    if (length(envirs) != 1L) {
        obs <- unique(lapply(envirs, as.list, all.names = TRUE, sorted = TRUE))
        if (length(obs) != 1L)
            return("differing objects within environments of wavelets")
    }
    deltajs <- unique(vapply(object, deltaj, 0))
    if (length(deltajs) != 1L)
        return(paste0("differing 'deltaj' values: ", paste0(signif(deltajs, getOption("digits")), collapse = ", ")))
})


# setOldClass(c("CompatibleWaveletTransforms", "list"), S4Class = "CompatibleWaveletTransforms")










setClass(
    Class    = "WaveletCovariance",
    contains = c("matrix", "oldClass"),
    slots    = c(wavelet = WaveletClassSymbol, frequency = "numeric",
        deltaj = "numeric"),
    prototype = prototype(.S3Class = c("WaveletCovariance", "matrix")))


setValidity("WaveletCovariance", function (object)
{
    if (!is.numeric_or_complex(object))
        return("object must be a numeric or complex matrix")
    N <- nrow(object)
    if (N < 3L)
        return("object must have at least 3 rows")
    validObject(wavelet(object))
    frequency <- frequency(object)
    if (length(frequency) != 1L)
        return("slot \"frequency\" must be of length 1")
    else if (!is.finite(frequency) || frequency <= 0)
        return("slot \"frequency\" must be finite and positive")
    deltaj <- deltaj(object)
    if (length(deltaj) != 1L)
        return("slot \"deltaj\" must be of length 1")
    else if (!is.finite(deltaj) || deltaj <= 0)
        return("slot \"deltaj\" must be finite and positive")
    j <- J(N = N, deltaj = deltaj(object))
    nc <- ncol(object)
    if (ncol(object) != j + 1)
        return(sprintf("object should have %d column%s, has %d column%s",
            j + 1, if (j) "s" else "", nc, if (nc != 1L) "s" else ""))
})


# setOldClass(c("WaveletCovariance", "matrix"), S4Class = "WaveletCovariance")
