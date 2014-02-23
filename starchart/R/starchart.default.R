starchart.default <-
function(x, y,
                             splittable = 10,
                             fitted = FALSE, ...) {
    x = as.data.frame(x)
    y = as.numeric(y)
    
    tree = starchartBuildTree(x, y, splittable = splittable)
    
    # (Optionally) compute a few additional fields for some generic functions.
    if (fitted) {
        tree$fitted.values = predict.starchart(tree, data = x)
        tree$residuals = y - tree$fitted.values
    }
    tree$call = match.call()
    
    class(tree) = "starchart"
    return(tree)
}
