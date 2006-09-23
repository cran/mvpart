"rpart.dist" <-
function (y, offset, parms, wt) 
{
    if (!is.null(offset)) 
        y <- y - offset
        ny <-  1
        list(y = y, parms = 0, numresp = ny, numy = ny, summary = function(yval, 
        dev, wt, ylevel, digits) {
        paste("  Means=", apply(formatg(yval, digits - 3), 1, 
            paste, collapse = ",", sep = ""), ", Summed MSE=", 
            formatg(dev/wt, digits), sep = "")
    }, text = function(yval, dev, wt, ylevel, digits, n, use.n) {
        if (use.n) {
            paste(formatg(dev, digits), " : n=", 
                n, sep = "")
        } 
          else {
            paste(formatg(dev, digits))
        }
    }, bar = function(yval2) {
            yval2
    })
}

