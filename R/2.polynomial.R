

hermite <- function (n)
{
    k <- 1
    for (i in seq.int(0L, length.out = n)) {
        k <- c(0, k) - c(k[-1L] * seq_len(i), 0, 0)
    }
    if (any(!is.finite(k)))
        stop(paste0(strwrap("invalid 'n', one or more hermite polynomial coefficients are not finite",
            exdent = 4), collapse = "\n"))
    if (any(abs(k) > .Machine$double.base^.Machine$double.digits - 1))
        warning(paste0(strwrap("exceedingly large 'n', one or more hermite polynomial coefficients exceed the largest integer valued floating-point number",
            exdent = 4), collapse = "\n"))
    k
}


as.body <- function (x, xname = "x", ...)
{
    if (!is.call(xname))
        xname <- as.symbol(xname)
    x <- as.numeric_or_complex(x)
    k <- as.double(seq.int(0L, along.with = x))
    i <- is.na(x) | x  # remove 0 from coefficients
    x <- x[i]
    k <- k[i]
    if (!length(x))
        return(substitute(rep(0, length(xname)), list(xname = xname)))
    if (length(x) == 1L) {
        num <- if (k == 1)
            xname
        else call("^", xname, k)
        body <- if (is.na(x))
            call("*", x, num)
        else if (x == 1)
            num
        else if (x == -1)
            call("-", num)
        else call("*", x, num)
    }
    else {
        nums <- lapply(X = k, FUN = function(y) if (y == 1)
            xname
        else if (y == 0)
            NULL
        else call("^", xname, y))
        .sign <- function(x) {
            if (is.numeric(x)) {
                value <- sign(x)
                value[is.na(value) | value == 0] <- 1
            }
            else {
                value <- sign(Re(x))
                i <- is.na(value) | value == 0
                value2 <- sign(Im(x[i]))
                value2[is.na(value2) | value2 == 0] <- 1
                value[i] <- value2
            }
            value
        }
        signs <- .sign(x[-1L])
        parts <- .mapply(function(e1, e2) {
            if (is.null(e2))
                as.numeric_or_complex(e1)
            else if (!is.na(e1) && e1 == 1)
                e2
            else call("*", as.numeric_or_complex(e1), e2)
        }, list(x * c(1, signs), nums), NULL)
        i <- signs == -1
        signs[i] <- "-"
        signs[!i] <- "+"
        body <- parts[[1L]]
        for (i in seq_along(signs)) {
            body <- call(signs[i], body, parts[[i + 1L]])
        }
    }
    body
}