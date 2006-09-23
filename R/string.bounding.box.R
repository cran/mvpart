"string.bounding.box" <-
function(s)
{
    s2 <- strsplit(s, "\n")
    rows <- sapply(s2, length)
    columns <- sapply(s2, function(x) max(nchar(x)))
    list(columns=columns, rows=rows)
}

