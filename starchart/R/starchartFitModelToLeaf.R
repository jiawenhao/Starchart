starchartFitModelToLeaf <-
function(x, y) {
    data = data.frame(x, response = y)
    # NOTE: Change the lm model below to customize leaf node models.
    lm(response ~ 1, data)
}
