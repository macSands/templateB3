#' Hello
#'
#' Prints a greeting message.
#'
#' @param name A character string; defaults to "World".
#' @return Prints a greeting message.
#' @export
hello <- function(name = "World") {
  message(sprintf("Hello, %s!", name))
}
