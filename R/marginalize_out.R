#' @export
marginalize_out = function(a, index){

  indices = 1:length(dim(a))


  apply(a, indices[-index], mean)
}
