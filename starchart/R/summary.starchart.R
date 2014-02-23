summary.starchart <-
function(object, ...) {
    result = list(call = object$call,
                  nodes = length(object$nodes),
                  leaves = sum(object$nodes == ""))
    if ("residuals" %in% names(object))
        result$residuals = object$residuals
    class(result) = "summary.starchart"
    return(result)
}
