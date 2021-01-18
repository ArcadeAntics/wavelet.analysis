

setClassUnion(
    name    = "numeric_or_complex",
    members = c("numeric", "complex"))





as.numeric_or_complex <- function (x)
UseMethod("as.numeric_or_complex")


as.numeric_or_complex.character <- function (x)
{
    x <- as.complex(x)
    if (all(is.na(x) | Im(x) == 0))
        as.numeric(x)
    else x
}


as.numeric_or_complex.complex <- function (x)
if (all(is.na(x) | Im(x) == 0)) as.numeric(x) else as.complex(x)


as.numeric_or_complex.default <- function (x)
as.numeric(x)


as.numeric_or_complex.ts <- function (x)
if (is.complex(x)) as.numeric_or_complex.complex(x) else as.numeric(x)





is.numeric_or_complex <- function (x)
is.numeric(x) || is.complex(x)





aslength1 <- function (x)
{
    len <- length(x)
    if (len == 1L) {
        x
    }
    else if (len > 1L) {
        warning(gettextf("first element used of '%s' argument",
            deparse(substitute(x), nlines = 1L)[1L], domain = NA))
        x[1L]
    }
    else stop(gettextf("'%s' must be of length 1", domain = NA,
        deparse(substitute(x), nlines = 1L)[1L]))
}
