#' Pone el origen de un vector en otro lado,
#' moviendo los elementos despues del borde al principio
#' y los de antes del borde quedan al final
#' @export
recenter.vector <- function(v, i){
  # c(
  #   v[(i+1):length(v)],
  #   v[1:i]
  # )
  c(
    v[(1:length(v)) >= i],
    v[(1:length(v)) < i]
  )
}
