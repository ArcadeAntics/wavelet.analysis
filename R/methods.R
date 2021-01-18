

sqrtgamma_max <- local({
    value <- 2
    increment <- 1
    fun <- function(x) sqrt(gamma(1 + x%%1)) * prod(sqrt((1 + x%%1):(x - 1)))
    while (value != value + increment) {
        if (is.finite(fun(value + increment)))
            value <- value + increment
        else increment <- increment/10
    }
    value
})





.coi <- function() NULL
formals(.coi) <- `names<-`(rep(list(quote(expr = )),
    4L), c("N", "deltat", FourierWavelengthAttrSymbol, EFoldingTimeAttrSymbol))
body(.coi) <- substitute({
    N <- (N - 2)/2
    c(1e-05, seq_len(ceiling(N)), seq.int(floor(N), 1, length.out = floor(N)),
        1e-05) * (deltat * FourierWavelengthAttrSymbol * EFoldingTimeAttrSymbol)
}, list(FourierWavelengthAttrSymbol = as.symbol(FourierWavelengthAttrSymbol),
    EFoldingTimeAttrSymbol = as.symbol(EFoldingTimeAttrSymbol)))





.CompatibleWaveletTransforms <- function (xx, cl = c("CompatibleWaveletTransforms", "list"))
{
    class(xx) <- cl
    xx
}


.period <- function (N, deltat, deltaj)
2 * deltat * 2^(0:J(N = N, deltaj = deltaj) * deltaj)





.scale <- function ()
NULL


formals(.scale) <- `names<-`(rep(list(quote(expr = )),
    2L), c("period", FourierWavelengthAttrSymbol))


body(.scale) <- substitute(period/sub1, list(sub1 = as.symbol(FourierWavelengthAttrSymbol)))





.WaveletTransform <- function (xx, cl = c("WaveletTransform", "matrix"))
{
    class(xx) <- cl
    xx
}


`[.CompatibleWaveletTransforms` <- function (x, i, ..., drop = FALSE)
{
    if (missing(i))
        return(x)
    value <- NextMethod("[")
    if (any(vapply(X = value, FUN = "is.null", FUN.VALUE = NA)))
        stop("subscript out of bounds")
    if (drop && length(value) == 1L)
        .subset2(value, 1L)
    else .CompatibleWaveletTransforms(value)
}


all.equal.Wavelet <- function (target, current, ...)
{
    x <- all.equal.default(target, current)
    if (!isTRUE(x))
        x
    else {
        x <- lapply(list(environment(target), environment(current)),
            as.list, all.names = TRUE, sorted = TRUE)
        if (identical(x[[1L]], x[[2L]]))
            TRUE
        else "objects within environments of target, current do not match"
    }
}


as.function.Wavelet <- function (x, ...)
unclass(x)


as.ts.WaveletTransform <- function (x, ...)
stats::ts(series(x), start = 0, frequency = frequency(x))





as.Wavelet <- function(x) NULL
body(as.Wavelet) <- substitute({
    if (!is.Wavelet(x))
        sub1
    cl <- oldClass(x)
    i <- match(WaveletClassSymbol, cl)
    if (i > 1L)
        class(x) <- cl[-(1:(i - 1L))]
    x
}, list(sub1 = call("stop", call("gettextf", sprintf("cannot coerce class %%s to a %s",
    WaveletClassSymbol), quote(sQuote(deparse(class(x))[1L])), domain = NA)),
    WaveletClassSymbol = WaveletClassSymbol))





CalculateAngularFrequency <- function (N, deltat = 2 * pi/N)
2 * pi/(N * deltat) * c(0:floor(N/2), seq.int(-floor((N - 1)/2),
    -1, length.out = ceiling(N/2 - 1)))


CalculateBound <- function (f, limit)
{
    if (abs(f(limit)) != 0)
        return(limit + sign(limit))
    value <- which(f(limit:0) == 0)
    value <- which(value == seq_along(value))
    limit + (1 - value[length(value)]) * sign(limit)
}


CalculateBound2 <- function (f, ...)
CalculateBound(function(x) f(x)^2, ...)


# CalculateFourierTransform <- function(f, omegas = seq(-20, 20, 0.01)) NULL
# body(CalculateFourierTransform) <- substitute({
#     safe_exp <- function(x) {
#         value <- exp(x)
#         value[is.na(value)] <- 0
#         value
#     }
#     i <- 1i
#     a <- function(t, omega) Re(f(t) * safe_exp(-i * omega * t))
#     b <- function(t, omega) Im(f(t) * safe_exp(-i * omega * t))
#     a <- vapply(omegas, function(omega) integrate(a, -Inf, Inf,
#         omega = omega)$value, 0)
#     b <- vapply(omegas, function(omega) integrate(b, -Inf, Inf,
#         omega = omega)$value, 0)
#     envir <- new.env(parent = baseenv())
#     envir$omegas <- omegas
#     envir$maxomega <- max(omegas)
#     envir$minomega <- min(omegas)
#     envir$y <- if (all(b == 0))
#         a
#     else a + b * i
#     value <- function(omega) {
#         omega
#         i <- vapply(omega, function(x) if (is.na(x))
#             NA_integer_
#         else which.min(abs(omegas - x)), 0L)
#         value <- y[i]
#         value[is.infinite(omega) | omega > maxomega | omega <
#             minomega] <- 0
#         value
#     }
#     environment(value) <- envir
#     Normal(value) <- CalculateNormal(value, use.stats = FALSE)
#     class(value) <- FourierTransformClass
#     return(value)
# }, list(FourierTransformClass = S3Class(new(FourierTransformClassSymbol))))


CalculateNormal <- function (f, use.stats = TRUE)
{
    if (use.stats) {
        1/sqrt(stats::integrate(function(x) abs(f(x))^2, -Inf, Inf)$value)
    }
    else {
        lower <- CalculateBound2(f, -100L)
        upper <- CalculateBound2(f, 100L)
        delta <- 0.001
        x <- seq(lower, upper, delta)
        y <- abs(f(x))^2
        area <- sum(y[-length(y)] + y[-1L])/2 * delta
        1/sqrt(area)
    }
}





CalculateReconstructionFactor <- function (wavelet, use.stats = TRUE)
{
    FT <- FourierTransform(wavelet)
    if (use.stats) {
        k <- sqrt(2 * pi) * Normal(FT)
        stats::integrate(function(x) abs(k * FT(x))^2/abs(x),
            -Inf, Inf)$value/1.3696
    }
    else {
        lower <- CalculateBound2(FT, -10000)
        upper <- CalculateBound2(FT, 10000)
        delta <- min((upper - lower)/1e6, 0.001)
        omega <- seq.int(lower, upper, delta)
        omega <- omega[omega != 0]
        sum(abs(NormalizeFourierTransform(FT, omega))^2/abs(omega))/1.3696 * delta
    }
}





CalculateReconstructionFactor2 <- function(wavelet, N = 4096, deltaj = 0.125) NULL
body(CalculateReconstructionFactor2) <- substitute({
    omega <- CalculateAngularFrequency(N, 1)
    s <- `.scale() sub`
    FT <- FourierTransform(wavelet)
    deltaj/Re(Normal(wavelet) * wavelet(0)) * sum(Re(vapply(X = s,
        FUN = function(sj) Conj(mean(NormalizeFourierTransform(FT,
            omega, sj))), FUN.VALUE = NA_complex_))/sqrt(s))
}, list(`.scale() sub` = `names<-`(quote(.scale(
          .period(N = N, deltat = deltat, deltaj = deltaj), FourierWavelength(wavelet))),
    c("", "period"                                        , FourierWavelengthAttrSymbol))))





coi <- function (x)
NULL


body(coi) <- `names<-`(call(".coi", quote(nrow(x)),
    quote(deltat(x)), quote(FourierWavelength(wavelet(x))), quote(EFoldingTime(wavelet(x)))),
    c("", "N", "deltat", FourierWavelengthAttrSymbol,
        EFoldingTimeAttrSymbol))





deltaj <- function (x)
attr(x, "deltaj")


deltat.WaveletCovariance <- function (x, ...)
1/attr(x, "frequency")


deltat.WaveletTransform <- deltat.WaveletCovariance





DOG <- function(m = 2L) NULL
body(DOG) <- substitute({
    m <- as.integer(m)
    m <- aslength1(m)
    if (is.na(m) || m <= 0L)
        stop("argument 'm' must be a positive number")
    FT <- function(omega) NULL
    body(FT) <- substitute({
        value <- `omega^m sub` * exp(-omega^2/2)
        value[is.na(value) & !is.na(omega)] <- 0
        value
    }, list(`omega^m sub` = if (m == 1L) as.symbol("omega") else substitute(omega^m, list(m = as.numeric(m)))))
    Normal(FT) <- as.numeric_or_complex(-1i^m/sqrtgamma(m + 0.5))
    value <- function(t) NULL
    body(value) <- substitute({
        value <- `-hermite(m) sub` * exp(-t^2/2)
        value[is.na(value) & !is.na(t)] <- 0
        value
    }, list(`-hermite(m) sub` = as.body(-hermite(m), xname = "t")))
    attr(value, "derivative") <- m
    Normal(value) <- 1/sqrtgamma(m + 0.5)
    FourierTransform(value) <- FT
    EFoldingTime(value) <- sqrt(2)
    FourierWavelength(value) <- 2 * pi/sqrt(m + 0.5)
    ReconstructionFactor(value) <- 2 * sqrt(pi) * beta(m, 0.5)/1.3696
    class(value) <- DOGWaveletClass
    return(value)
}, list(DOGWaveletClass = S3Class(new(DOGWaveletClassSymbol))))





EFoldingTime <- function(x) NULL
body(EFoldingTime) <- substitute(attr(x, EFoldingTimeAttrSymbol),
    list(EFoldingTimeAttrSymbol = EFoldingTimeAttrSymbol))


`EFoldingTime<-` <- function(x, value) NULL
body(`EFoldingTime<-`) <- substitute(`attr<-`(x, EFoldingTimeAttrSymbol, value),
    list(EFoldingTimeAttrSymbol = EFoldingTimeAttrSymbol))





FourierTransform <- function(x) NULL
body(FourierTransform) <- substitute(attr(x, FourierTransformAttrSymbol),
    list(FourierTransformAttrSymbol = FourierTransformAttrSymbol))


`FourierTransform<-` <- function(x, value) NULL
body(`FourierTransform<-`) <- substitute({
    class(value) <- FourierTransformClass
    `attr<-`(x, FourierTransformAttrSymbol, value)
}, list(FourierTransformClass = S3Class(new(FourierTransformClassSymbol)),
    FourierTransformAttrSymbol = FourierTransformAttrSymbol))





FourierWavelength <- function(x) NULL
body(FourierWavelength) <- substitute(attr(x, FourierWavelengthAttrSymbol),
    list(FourierWavelengthAttrSymbol = FourierWavelengthAttrSymbol))


`FourierWavelength<-` <- function(x, value) NULL
body(`FourierWavelength<-`) <- substitute(`attr<-`(x, FourierWavelengthAttrSymbol, value),
    list(FourierWavelengthAttrSymbol = FourierWavelengthAttrSymbol))





frequency.WaveletCovariance <- frequency.WaveletTransform <- function (x, ...)
attr(x, "frequency")





Gabor <- function(sigma = 6) NULL
body(Gabor) <- substitute({
    sigma <- as.numeric(sigma)
    sigma <- aslength1(sigma)
    if (!is.finite(sigma) || sigma == 0)
      stop("argument 'sigma' must be finite and non-zero")
    isigma <- 1i * sigma
    kappa <- exp(-sigma^2/2)
    FT <- function(omega) exp(-(sigma - omega)^2/2) - kappa * exp(-omega^2/2)
    Normal(FT) <- pi^(-1/4)/sqrt(1 + exp(-sigma^2) - 2 * exp(-3/4 * sigma^2))
    value <- function(t) {
        value <- exp(-t^2/2) * (exp(isigma * t) - kappa)
        value[is.na(value) & !is.na(t)] <- 0
        value
    }
    attr(value, "frequency") <- sigma
    Normal(value) <- Normal(FT)
    FourierTransform(value) <- FT
    EFoldingTime(value) <- sqrt(2)
    FourierWavelength(value) <- 4 * pi/(sigma + sqrt(2 + sigma^2))
    ReconstructionFactor(value) <- CalculateReconstructionFactor(value)
    class(value) <- GaborWaveletClass
    return(value)
}, list(GaborWaveletClass = S3Class(new(GaborWaveletClassSymbol))))





is.CompatibleWaveletTransforms <- function (x)
inherits(x, "CompatibleWaveletTransforms")





is.Wavelet <- function(x) NULL
body(is.Wavelet) <- substitute(inherits(x, WaveletClassSymbol),
    list(WaveletClassSymbol = WaveletClassSymbol))





is.WaveletCovariance <- function (x)
inherits(x, "WaveletCovariance")


is.WaveletTransform <- function (x)
inherits(x, "WaveletTransform")


isCompatibleWaveletTransforms <- function (x, y)
{
    if (!is.WaveletTransform(x))
        stop("argument 'x' should be or extend class \"WaveletTransform\"")
    if (!is.WaveletTransform(y))
        stop("argument 'y' should be or extend class \"WaveletTransform\"")
    methods::validObject(x)
    methods::validObject(y)
    if (nrow(x) != nrow(y))
        stop("differing number of rows: ", nrow(x), ", ", nrow(y))
    if (ncol(x) != ncol(y))
        stop("differing number of columns: ", ncol(x), ", ",
            ncol(y))
    if (deltat(x) != deltat(y))
        stop("differing 'deltat' values: ", signif(deltat(x),
            getOption("digits")), ", ", signif(deltat(y), getOption("digits")))
    if (length(unique(list(deparse(wavelet(x)), deparse(wavelet(y))))) !=
        1L)
        stop("differing wavelets when deparsed")
    if (length(unique(list(as.list(environment(wavelet(x)), all.names = TRUE,
        sorted = TRUE), as.list(environment(wavelet(y)), all.names = TRUE,
        sorted = TRUE)))) != 1L)
        stop("differing objects within environments of wavelets")
    if (deltaj(x) != deltaj(y))
        stop("differing 'deltaj' values: ", signif(deltaj(x),
            getOption("digits")), ", ", signif(deltaj(y), getOption("digits")))
    TRUE
}


J <- function (N, deltaj)
as.integer(1/deltaj * log2((N - 1)/2))





setWaveletSubClass <- function(Class, slots = list()) NULL
body(setWaveletSubClass) <- substitute({
    if (!is.vector(Class, "character"))
        stop("argument 'Class' must be class \"character\"")
    Class <- aslength1(Class)
    if (is.na(Class))
        stop("argument 'Class' must not be NA")
    if (grepl("\\s", Class))
        stop("argument 'Class' should not contain white-space characters")
    if (!nzchar(Class))
        stop("attempt to use zero-length variable name")
    subclasses <- methods::getClassDef(WaveletClassSymbol)@subclasses
    matched <- match(tolower(Class), tolower(names(subclasses)),
        nomatch = 0L)
    if (matched) {
        if (subclasses[[matched]]@package == .packageName)
            stop("cannot overwrite current definition of Wavelet subclass \"",
                names(subclasses)[[matched]], "\"")
        else warning("overwriting current definition of Wavelet subclass \"",
            names(subclasses)[[matched]], "\"")
    }
    where <- topenv(parent.frame())
    methods::setClass(
        Class = Class,
        contains = WaveletClassSymbol,
        slots = slots,
        prototype = methods::prototype(.S3Class = c(Class, WaveletClassSymbol, "function")),
        where = where)
    methods::getClassDef(Class, where = where)
}, list(WaveletClassSymbol = WaveletClassSymbol, .packageName = .packageName))





Normal <- function(x) NULL
body(Normal) <- substitute(attr(x, NormalAttrSymbol), list(NormalAttrSymbol = NormalAttrSymbol))


`Normal<-` <- function(x, value) NULL
body(`Normal<-`) <- substitute(`attr<-`(x, NormalAttrSymbol, as.numeric_or_complex(value)),
    list(NormalAttrSymbol = NormalAttrSymbol))





NormalizeFourierTransform <- function() NULL
formals(NormalizeFourierTransform) <- local({
    value <- list(quote(expr = ), quote(expr = ), 1, 1)
    names(value) <- c(FourierTransformAttrSymbol, "omegak", "s" ,"deltat")
    value
})
body(NormalizeFourierTransform) <- substitute(sqrt(2 * pi * s/deltat) *
    Normal(FourierTransformAttrSymbol) * FourierTransformAttrSymbol(s *
    omegak), list(FourierTransformAttrSymbol = as.symbol(FourierTransformAttrSymbol)))





NormalizeWavelet <- function (wavelet, t, s = 1, deltat = 1)
sqrt(deltat/s) * Normal(wavelet) * wavelet(t * deltat/s)





Paul <- function(m = 4L) NULL
body(Paul) <- substitute({
    m <- as.integer(m)
    m <- aslength1(m)
    if (is.na(m) || m <= 0L)
        stop("argument 'm' must be a positive number")
    else if (m > max)
        stop("argument 'm' is too large a number")
    i <- 1i
    FT <- function(omega) NULL
    body(FT) <- substitute({
        value <- `omega^m sub` * exp(-omega)
        value[is.na(value) & !is.na(omega) | omega < 0] <- 0
        value
    }, list(`omega^m sub` = if (m == 1L) as.symbol("omega") else substitute(omega^m,
        list(m = as.numeric(m)))))
    Normal(FT) <- 2^m/sqrt(m)/sqrtgamma(2 * m)
    value <- function(t) NULL
    body(value) <- substitute((1 - i * t)^-`m + 1 sub`, list(`m + 1 sub` = m + 1))
    attr(value, "order") <- m
    Normal(value) <- 2^m * 1i^m * factorial(m)/sqrt(pi)/sqrtfactorial(2 * m)
    FourierTransform(value) <- FT
    EFoldingTime(value) <- 1/sqrt(2)
    FourierWavelength(value) <- 4 * pi/(2 * m + 1)
    ReconstructionFactor(value) <- 2 * pi/(m * 1.3696)
    class(value) <- PaulWaveletClass
    return(value)
}, list(PaulWaveletClass = S3Class(new(PaulWaveletClassSymbol)),
    max = as.integer(sqrtgamma_max/2)))





period <- function (x)
.period(N = nrow(x), deltat = deltat(x), deltaj = deltaj(x))


plot.FourierTransform <- function (x, y = 0, to = 1, from = y, n = 101, add = FALSE, type = "l",
    xlab = NULL, ylab = NULL, log = NULL, xlim = NULL, ylim = NULL,
    col = c("orange2", "steelblue2"), bg = graphics::par("bg"),
    pch = graphics::par("pch"), cex = graphics::par("cex"), lty = graphics::par("lty"),
    lwd = graphics::par("lwd"), ...)
{
    fun <- x
    if (grDevices::dev.cur() == 1L && !isFALSE(add)) {
        warning("'add' will be ignored as there is no existing plot")
        add <- FALSE
    }
    addF <- isFALSE(add)
    if (is.null(xlab)) {
        xlab <- if (is.Wavelet(x))
            as.symbol("t")
        else as.symbol("omega")
    }
    if (is.null(ylab)) {
        ylab <- if (is.Wavelet(x))
            call("*", as.symbol("Psi"), call("(", xlab))
        else call("*", call("hat", as.symbol("Psi")), call("(", xlab))
    }
    if (is.null(from) || is.null(to)) {
        xl <- if (!is.null(xlim))
            xlim
        else if (!addF) {
            pu <- graphics::par("usr")[1L:2L]
            if (graphics::par("xaxs") == "r")
                pu <- grDevices::extendrange(pu, f = -1/27)
            if (graphics::par("xlog"))
                10^pu
            else pu
        }
        else c(0, 1)
        if (is.null(from))
            from <- xl[1L]
        if (is.null(to))
            to <- xl[2L]
    }
    lg <- if (length(log))
        log
    else if (!addF && graphics::par("xlog"))
        "x"
    else ""
    if (length(lg) == 0)
      lg <- ""
    if (grepl("x", lg, fixed = TRUE)) {
        if (from <= 0 || to <= 0)
            stop("'from' and 'to' must be > 0 with log=\"x\"")
        x <- exp(seq.int(log(from), log(to), length.out = n))
    }
    else x <- seq.int(from, to, length.out = n)
    y <- fun(x)
    if (!is.null(Normal(fun)))
        y <- y * Normal(fun)
    if (is.null(ylim))
        ylim <- range(Re(y), Im(y), na.rm = TRUE)
    type <- rep(type, length.out = 2L)
    col  <- rep(col , length.out = 2L)
    bg   <- rep(bg  , length.out = 2L)
    pch  <- rep(pch , length.out = 2L)
    cex  <- rep(cex , length.out = 2L)
    lty  <- rep(lty , length.out = 2L)
    lwd  <- rep(lwd , length.out = 2L)
    if (is.numeric(y)) {
        if (isTRUE(add))
            graphics::lines(x = x, y = y, type = type[[1L]],
                col = col[[1L]], bg = bg[[1L]], pch = pch[[1L]],
                cex = cex[[1L]], lty = lty[[1L]], lwd = lwd[[1L]],
                ...)
        else graphics::plot(x = x, y = y, type = type[[1L]],
            col = col[[1L]], bg = bg[[1L]], pch = pch[[1L]],
            cex = cex[[1L]], lty = lty[[1L]], lwd = lwd[[1L]],
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
            log = lg, ...)
    }
    else {
        if (isTRUE(add))
            graphics::lines(x = x, y = Im(y), type = type[[2L]],
                col = col[[2L]], bg = bg[[2L]], pch = pch[[2L]],
                cex = cex[[2L]], lty = lty[[2L]], lwd = lwd[[2L]],
                ...)
        else graphics::plot(x = x, y = Im(y), type = type[[2L]],
            col = col[[2L]], bg = bg[[2L]], pch = pch[[2L]],
            cex = cex[[2L]], lty = lty[[2L]], lwd = lwd[[2L]],
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
            log = lg, ...)
        graphics::lines(x = x, y = Re(y), type = type[[1L]],
            col = col[[1L]], bg = bg[[1L]], pch = pch[[1L]],
            cex = cex[[1L]], lty = lty[[1L]], lwd = lwd[[1L]],
            ...)
    }
    invisible(list(x = x, y = y))
}


plot.Wavelet <- plot.FourierTransform


legend.dimensions <- function (expr, local = TRUE)
{
    # `legend.dimensions` returns the width and height of the legend's box produced by 'expr'
    #
    # `expr` is a call to function `legend`
    # `local` controls the environment in which `expr` is evaluated
    #     TRUE means that `expr` will be evaluated in the same
    #         environment as the call to `legend.dimensions`
    #     FALSE means that `expr` will be evaluated in the global environment
    #     can also be an environment to evaluate `expr`
    #
    # return value is a list with components w, h
    #     giving width and height of the legend's box IN INCHES
    #
    # `legend.dimensions` is useful for knowing how wide a legend's box is BEFORE PLOTTING
    # the function `legend` relies on a plot already existing before calling.
    # if a plot is already produced, `legend.dimensions` is not useful
    #
    # I've used `legend.dimensions` as seen here:
    # expr <- expression(legend(  # the call that produces the desired legend
    #     x = par("usr")[2L], y = par("usr")[4L],
    #     legend = letters[1:5], fill = 1:5,
    #     xpd = TRUE))
    # x <- legend.dimensions(expr)  # capture the width and height of the resultant legend
    # y <- par("mai")  # the size of the margins, in inches
    # y[4L] <- x$w + 0.1  # make the right side margin exactly wide enough to fit this legend, plus 0.1 inches
    # par(mai = y)  # set the new margin sizes
    # plot(1:5)  # produce the desired plot
    # eval(expr)  # add the legend

    envir <- if (isTRUE(local))
        parent.frame()
    else if (isFALSE(local))
        globalenv()
    else if (is.environment(local))
        local
    else stop("'local' must be TRUE, FALSE or an environment")
    if (!is.expression(expr))
        stop("invalid 'expr', must be class \"expression\"")
    expr <- aslength1(expr)
    opar <- graphics::par("cex", "family", "font", "lheight", mar = rep(0, 4L), xlog = FALSE, ylog = FALSE)  # capture the graphical parameters that may affect the legend size
    on.exit(graphics::par(opar))
    filename <- tempfile()
    on.exit(suppressWarnings(file.remove(filename)),  # before exiting the function, remove the file
        add = TRUE, after = FALSE)
    old.device <- grDevices::dev.cur()
    on.exit(grDevices::dev.set(old.device), add = TRUE, after = FALSE)
    grDevices::png(filename, width = 480, height = 480, res = 48)
    num.device <- grDevices::dev.cur()
    on.exit(grDevices::dev.off(num.device), add = TRUE, after = FALSE)  # shut down the plotting device, BEFORE attempting to remove the file
    graphics::par(opar)
    graphics::plot.default(
        xlim = c(0, 10), ylim = c(0, 10), xaxs = "i", yaxs = "i",
        x = NA_real_, y = NA_real_,
        axes = FALSE, frame.plot = FALSE, ann = FALSE
    )
    value <- eval(expr, envir)$rect
    list(w = value$w/graphics::xinch(warn.log = FALSE), h = value$h/graphics::yinch(warn.log = FALSE))
}


plot.WaveletCovariance <- function (x, y = c("time", "frequency", "2d"), type = "l", main = NULL, sub = NULL, xlab = NULL,
    ylab = NULL, ann = graphics::par("ann"), axes = TRUE, frame.plot = axes, ..., col = graphics::par("col"), legend.n = 5L,
    digits = max(3L, getOption("digits") - 3L), filter.time = NULL, filter.frequency = NULL)
{
    how <- match.arg(y)
    switch(how, `2d` = {
        if (is.null(xlab))
            xlab <- "time"
        if (is.null(ylab))
            ylab <- "period"
        if (missing(col) || is.null(col))
            col <- grDevices::heat.colors(12)
        legend.legend <- signif(seq.int(min(x, na.rm = TRUE),
            max(x, na.rm = TRUE), length.out = legend.n), digits)
        legend.command <- expression(graphics::legend(
            x = graphics::par("usr")[2L], y = graphics::par("usr")[4L],
            legend = legend.legend,
            fill = col[seq.int(1, length(col), length.out = length(legend.legend))],
            bty = "n", xpd = TRUE))
        ld <- legend.dimensions(legend.command)
        omai <- graphics::par("mai")
        on.exit(graphics::par(mai = omai))
        mai <- omai
        mai[4L] <- ld$w + 0.1
        graphics::par(mai = mai)
        z <- x
        m <- nrow(z)
        n <- ncol(z)
        x <- seq_len(m)
        y <- seq_len(n)
        graphics::image(x, y, z, xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab,
            col = col, ...)
        graphics::box()
        graphics::axis(1L, x[seq.int(1L, m, length.out = 5L)], signif(time(z),
            2L)[seq.int(1L, m, length.out = 5L)], ...)
        graphics::axis(2L, y[seq.int(1L, n, length.out = 5L)], signif(period(z),
            2L)[seq.int(1L, n, length.out = 5L)], ...)
        eval(legend.command)
    }, time = {
        if (is.null(xlab))
            xlab <- "time"
        if (is.null(ylab))
            ylab <- "spectrum"
        plot(time(x), rowSums(x)/deltat(x), xlab = xlab, ylab = ylab,
            col = col, type = type, ...)
    }, frequency = {
        if (is.null(xlab))
            xlab <- "frequency"
        if (is.null(ylab))
            ylab <- "spectrum"
        plot(period(x), colSums(x)/deltaj(x), xlab = xlab, ylab = ylab,
            col = col, type = type, ...)
    })
}


plot.WaveletTransform <- function (x, y = NULL, ..., xlab = "time", ylab = "series")
plot(as.ts(x), y = NULL, ..., xlab = xlab, ylab = ylab)


print.CompatibleWaveletTransforms <- function (x, ...)
{
    cat(sprintf("<List of %d compatible wavelet power spectra>\n\n",
        length(x)))
    print(unclass(x), ...)
    invisible(x)
}


printWaveletAttr <- function (x, digits = getOption("digits"))
{
    value <- list(`Wavelet normal` = Normal(x), `Fourier transform normal` = Normal(FourierTransform(x)),
        `e-folding time` = EFoldingTime(x), `Fourier wavelength` = FourierWavelength(x),
        `Reconstruction factor` = ReconstructionFactor(x))
    value <- lapply(X = value, FUN = "signif", digits = digits)
    cat(paste0(format(names(value)), " : ", value), sep = "\n")
    invisible(x)
}


assign(paste0("print.", DOGWaveletClassSymbol),
    function (x, digits = getOption("digits"), verbose = getOption("verbose"), ...)
{
    cat("DOG (Derivative Of Gaussian) wavelet of order ", attr(x,
        "derivative"), "\n", sep = "")
    if (isTRUE(verbose)) {
        cat("\n",
            "function (t)\n",
            "-`m-th hermite polynomial`(t) * exp(-t^2/2)\n",
            # paste0(deparse(body(x)[[2L]][[3L]]), "\n"),
            "\n",
            "Fourier transform:\n",
            "function (\u03c9)\n",
            "\u03c9^m * exp(-\u03c9^2/2)\n",
            "\n", sep = "")
        printWaveletAttr(x, digits = digits)
    }
    invisible(x)
})


assign(paste0("print.", GaborWaveletClassSymbol),
    function (x, digits = getOption("digits"), verbose = getOption("verbose"), ...)
{
    cat("Gabor (complex Morlet) wavelet of frequency ", as.character(signif(attr(x,
        "frequency"), digits = digits)), "\n", sep = "")
    if (isTRUE(verbose)) {
        cat("\n",
            "function (t)\n",
            "exp(-t^2/2) * (exp(i * \u03c3 * t) - exp(-\u03c3^2/2))\n",
            "\n",
            "Fourier transform:\n",
            "function (\u03c9)\n",
            "exp(-(\u03c3 - \u03c9)^2/2) - exp(-\u03c3^2/2) * exp(-\u03c9^2/2)\n",
            "\n", sep = "")
        printWaveletAttr(x, digits = digits)
    }
    invisible(x)
})


assign(paste0("print.", PaulWaveletClassSymbol),
    function (x, digits = getOption("digits"), verbose = getOption("verbose"), ...)
{
    cat("Paul wavelet of order ", attr(x, "order"), "\n", sep = "")
    if (isTRUE(verbose)) {
        cat("\n",
            "function (t)\n",
            "(1 - i * t)^-(m + 1)\n",
            "\n",
            "Fourier transform:\n",
            "function (\u03c9)\n",
            "Heaviside(\u03c9) * \u03c9^m * exp(-\u03c9)\n",
            "\n", sep = "")
        printWaveletAttr(x, digits = digits)
    }
    invisible(x)
})


printWaveletCovariance_or_WaveletTransform <- function (Class)
{
    value <- function(x, digits = getOption("digits"), ...) NULL
    body(value) <- substitute({
        cat(paste0("<", paste0(dim(x), collapse = " x "), sub))
        cat("\n")
        cat("\u03b4j : ", as.character(signif(deltaj(x), digits)), "\n", sep = "")
        cat("\u03b4t : ", as.character(signif(deltat(x), digits)), "\n", sep = "")
        cat("wavelet : ")
        print(wavelet(x), digits = digits, ...)
        invisible(x)
    }, list(sub = switch(Class, WaveletCovariance = " wavelet covariance matrix>\n",
        WaveletTransform = " wavelet power spectrum>\n")))
    value
}


print.WaveletCovariance <- printWaveletCovariance_or_WaveletTransform("WaveletCovariance")


print.WaveletTransform <- printWaveletCovariance_or_WaveletTransform("WaveletTransform")





ReconstructionFactor <- function(x) NULL
body(ReconstructionFactor) <- substitute(attr(x, ReconstructionFactorAttrSymbol),
    list(ReconstructionFactorAttrSymbol = ReconstructionFactorAttrSymbol))


`ReconstructionFactor<-` <- function(x, value) NULL
body(`ReconstructionFactor<-`) <- substitute({
    if (is.null(value))
        `attr<-`(x, ReconstructionFactorAttrSymbol, CalculateReconstructionFactor(x))
    else `attr<-`(x, ReconstructionFactorAttrSymbol, value)
}, list(ReconstructionFactorAttrSymbol = ReconstructionFactorAttrSymbol))





scale.WaveletCovariance <- function(x, ...) NULL
body(scale.WaveletCovariance) <- substituteDirect(local({
    value <- list(quote(period(x)), quote(FourierWavelength(wavelet(x))))
    names(value) <- c("period", FourierWavelengthAttrSymbol)
    as.call(c(list(as.symbol(".scale")), value))
}), list(FourierWavelengthAttrSymbol = FourierWavelengthAttrSymbol))


scale.WaveletTransform <- scale.WaveletCovariance





series <- function (x)
attr(x, "series")


Spectrum <- function (x, use.coi = FALSE)
{
    if (use.coi) {
        coi <- coi(x)
        period <- period(x)
        attributes(x) <- list(dim = dim(x))
        for (i in seq_along(period)) x[coi < period[[i]], i] <- NA
    }
    else attributes(x) <- list(dim = dim(x))
    x
}


sqrtfactorial <- function (x)
sqrtgamma(x + 1)





sqrtgamma <- function(x) NULL
body(sqrtgamma) <- substitute({
    value <- numeric(length(x))
    i <- is.na(x)
    value[i] <- x[i]
    j <- !i & x <= 2
    value[j] <- sqrt(gamma(x[j]))
    i <- i | j
    j <- !i & x > sqrtgamma_max
    value[j] <- Inf
    i <- !(i | j)
    x <- x[i]
    value[i] <- unlist(.mapply(function(y, z) {
        sqrt(gamma(y)) * prod(sqrt(y:z))
    }, list(y = 1 + x%%1, z = x - 1), NULL), recursive = FALSE)
    value
}, list(sqrtgamma_max = sqrtgamma_max))





time.WaveletCovariance <- function (x, offset = 0, ...)
{
    n <- nrow(x)
    xtsp <- c(0, (n - 1L)/attr(x, "frequency"), attr(x, "frequency"))
    y <- seq.int(xtsp[1L], xtsp[2L], length.out = n) + offset/xtsp[3L]
    stats::tsp(y) <- xtsp
    attr(y, "class") <- "ts"
    y
}


time.WaveletTransform <- time.WaveletCovariance


validAttr <- function (AttrSymbol)
{
    value <- function(x) NULL
    body(value) <- substitute({
        if (!is.numeric(x))
            sub1
        else if (length(x) != 1L)
            sub2
        else if (!is.finite(x) || x <= 0)
            sub3
    }, list(sub1 = sprintf("slot \"%s\" must be numeric", AttrSymbol),
        sub2 = sprintf("slot \"%s\" must be of length 1", AttrSymbol),
        sub3 = sprintf("slot \"%s\" must be finite and positive",
            AttrSymbol)))
    environment(value) <- parent.frame()
    value
}


validEFoldingTime <- validAttr(EFoldingTimeAttrSymbol)


validFourierWavelength <- validAttr(FourierWavelengthAttrSymbol)





validNormal <- function(x) NULL
body(validNormal) <- substitute({
    if (!is.numeric_or_complex(x))
        sub1
    else if (length(x) != 1L)
        sub2
    else if (!is.finite(x) || x == 0)
        sub3
}, list(
    sub1 = sprintf("slot \"%s\" must be numeric or complex", NormalAttrSymbol),
    sub2 = sprintf("slot \"%s\" must be of length 1", NormalAttrSymbol),
    sub3 = sprintf("slot \"%s\" must be finite and non-zero", NormalAttrSymbol)))





validReconstructionFactor <- validAttr(ReconstructionFactorAttrSymbol)


wavelet <- function (x)
attr(x, "wavelet")


WaveletCovariance <- function (x, ..., fun = "Re", use.coi = FALSE, simplify = TRUE)
{
    if (is.WaveletCovariance(x))
        return(sum(Spectrum(x, use.coi), na.rm = TRUE))
    if (!is.WaveletTransform(x) && !is.CompatibleWaveletTransforms(x)) {
        x <- WaveletTransform(x, ...)
        if (is.WaveletTransform(x)) {
            Wx <- x
            Wy <- x
        }
    }
    else if (is.WaveletTransform(x)) {
        Wx <- x
        if (missing(...)) {
            Wy <- Wx
        }
        else {
            Wy <- ...elt(1L)
            isCompatibleWaveletTransforms(Wx, Wy)
        }
    }
    if (is.CompatibleWaveletTransforms(x)) {
        Wx <- x[[1L]]
        Wy <- x[[2L]]
    }
    fun <- match.fun(fun)
    value <- fun(Wx * Conj(Wy))
    value <- sweep(x = value, MARGIN = 2L, STATS = scale(Wx), FUN = "/")
    value <- deltat(Wx) * deltaj(Wx)/(ReconstructionFactor(wavelet(Wx)) * nrow(Wx)) * value
    if (simplify)
        return(sum(Spectrum(value, use.coi), na.rm = TRUE))
    attrib <- attributes(value)
    attrib[c("series", "class")] <- list(NULL, c("WaveletCovariance", "matrix"))
    value <- Spectrum(value, use.coi)
    attributes(value) <- attrib
    return(value)
}


WaveletCovariancePlot <- function (x, ..., digits = max(3L, getOption("digits") - 3L),
    against = c("time", "period"), index, value, xlab = against,
    ylab = "spectrum")
{
    against <- match.arg(against)
    if (!is.WaveletCovariance(x))
        x <- WaveletCovariance(x, ..., SIMPLIFY = FALSE)
    if (against == "time") {
        index <- if (missing(index) && missing(value))
            seq_len(ncol(x))
        else if (!missing(value)) {
            period <- period(x)
            if (length(value) == 1L)
                which.min(abs(period - value))
            else if (length(value) == 2L)
                which(period >= value[1L] & period <= value[2L])
        }
        else if (length(index) == 2L)
            index[1L]:index[2L]
        data_ <- data.frame(x = time(x), y = (if (length(index) ==
            1L)
            x[, index]
        else rowSums(x[, index]))/deltat(x))
        maxdata_ <- signif(unlist(data_[which.max(data_$y), ]),
            digits)
        main <- sprintf("Maximum = %s, Time = %s, Flux = %s",
            maxdata_[2L], maxdata_[1L], signif(sum(data_$y) *
                deltat(x), digits))
    }
    else if (against == "period") {
        index <- if (missing(index) && missing(value))
            seq_len(nrow(x))
        else if (!missing(value)) {
            time <- time(x)
            if (length(value) == 1L)
                which.min(abs(time - value))
            else if (length(value) == 2L)
                which(time >= value[1L] & time <= value[2L])
        }
        else if (length(index) == 2L)
            index[1L]:index[2L]
        data_ <- data.frame(x = period(x), y = (if (length(index) ==
            1L)
            x[index, ]
        else colSums(x[index, ]))/deltaj(x))
        maxdata_ <- signif(unlist(data_[which.max(data_$y), ]),
            digits)
        main <- sprintf("Maximum = %s, Period = %s, Flux = %s",
            maxdata_[2L], maxdata_[1L], signif(sum(data_$y) *
                deltaj(x), digits))
    }
    plot(data_, type = "l", xlab = xlab, ylab = ylab, main = main)
    invisible(data_)
}





WaveletTransform <- function(x, wavelet = Gabor(), deltaj = 0.125, frequency = 1,
    deltat = 1, methodWaveletTransform = c("discrete", "continuous")) NULL
body(WaveletTransform) <- substitute({
    DWT <- function(x) {
        devx <- x - mean.default(x, na.rm = TRUE)
        devx <- ZeroPadding(devx, NPad)
        X <- stats::fft(devx)/NPad
        W <- matrix(vapply(s, function(sj) stats::fft(X * Conj(NormalizeFourierTransform(sub1,
            omega, sj, deltat)), inverse = TRUE)[seq_along(x)],
            complex(N)), nrow = N)
        attributes(W) <- list(dim = dim(W), series = x, wavelet = wavelet,
            frequency = 1/deltat, deltaj = deltaj)
        .WaveletTransform(W)
    }
    CWT <- function(x) {
        devx <- x - mean.default(x, na.rm = TRUE)
        W <- t(vapply(0:(N - 1L), function(n) colSums(x * conjPsiLookup(n)),
            complex(length(s))))
        attributes(W) <- list(dim = dim(W), series = x, wavelet = wavelet,
            frequency = 1/deltat, deltaj = deltaj)
        .WaveletTransform(W)
    }
    methodWaveletTransform <- match.arg(methodWaveletTransform)
    if (!is.Wavelet(wavelet))
        stop("'wavelet' must be or extend class \"Wavelet\"")
    if (!is.numeric_or_complex(x)) {
        if (is.logical(x)) {
            attrib <- list(dim = dim(x), dimnames = dimnames(x),
                names = names(x), tsp = stats::tsp(x))
            x <- as.numeric(x)
            attributes(x) <- attrib
        }
        else if (!is.data.frame(x))
            stop("'x' should be numeric or complex or a data.frame")
        else if (!all(vapply(x, is.numeric_or_complex, NA)))
            stop("all columns of 'x' must be numeric or complex")
    }
    N <- NROW(x)
    if (N < 3L)
        stop("'x' must contain at least 3 observations")
    if (inherits(x, "ts"))
        deltat <- stats::deltat(x)
    else {
        if (missing(frequency)) {
            if (!missing(deltat)) {
                if (!is.numeric(deltat) || length(deltat) !=
                  1L || !is.finite(deltat) || deltat <= 0)
                  stop("invalid 'deltat' value")
                frequency <- 1/deltat
            }
        }
        else if (!is.numeric(frequency) || length(frequency) !=
            1L || !is.finite(frequency) || frequency <= 0)
            stop("invalid 'frequency' value")
        if (frequency > 1 && abs(frequency - round(frequency)) <
            getOption("ts.eps"))
            frequency <- round(frequency)
        deltat <- 1/frequency
    }
    dl <- length(dim(x))
    x <- if (dl < 2L)
        c(x)
    else if (dl == 2L)
        as.data.frame(x, optional = TRUE)
    else stop("'x' cannot have more than two dimensions")
    s <- sub2
    switch(methodWaveletTransform, discrete = {
        sub1 <- FourierTransform(wavelet)
        NPad <- 2^ceiling(log2(N))
        omega <- CalculateAngularFrequency(NPad, deltat)
        if (is.null(dim(x)))
            DWT(x)
        else .CompatibleWaveletTransforms(lapply(x, DWT))
    }, continuous = {
        `n' - n` <- (1L - N):(N - 1L)
        conjPsi <- Conj(vapply(s, function(sj) NormalizeWavelet(wavelet,
            `n' - n`, sj, deltat), complex(length(`n' - n`))))
        conjPsiLookup <- function(n) conjPsi[(N - n):(N * 2 -
            1 - n), ]
        if (is.null(dim(x)))
            CWT(x)
        else .CompatibleWaveletTransforms(lapply(x, CWT))
    })
}, list(sub1 = as.symbol(FourierTransformAttrSymbol), sub2 = `names<-`(quote(.scale(
          .period(N = N, deltat = deltat, deltaj = deltaj), FourierWavelength(wavelet))),
    c("", "period"                                        , FourierWavelengthAttrSymbol))))





ZeroPadding <- function (x, length.out)
c(x, rep(0, length.out - length(x)))


remove(NormalAttrSymbol, FourierTransformAttrSymbol, EFoldingTimeAttrSymbol,
    FourierWavelengthAttrSymbol, ReconstructionFactorAttrSymbol,
    FourierTransformClassSymbol, WaveletClassSymbol, DOGWaveletClassSymbol,
    GaborWaveletClassSymbol, PaulWaveletClassSymbol, printWaveletCovariance_or_WaveletTransform,
    sqrtgamma_max, validAttr)





if (FALSE) {
    (function() {
        odir <- getwd()
        on.exit(setwd(odir))
        setwd("~")
        pkg <- "wavelet.analysis"
        system(sprintf("R CMD INSTALL --no-multiarch --with-keep.source %s", pkg))
        system(sprintf("R CMD build %s", pkg))
        system(sprintf("R CMD check --as-cran %s_%s.tar.gz", pkg, utils::packageVersion(pkg)))
    })()
}
