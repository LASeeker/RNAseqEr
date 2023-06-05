#' colour_palette
#' @description sets colour palette to my default
#' @import ggsci RColorBrewer
#' @return returns a palette of colours that can be used for plots
#' @export
#'
#' @examples
#' library(RColorBrewer)
#' library(ggsci)
#' my_cols <- colour_palette()
colour_palette<- function(){
  cols_new <- c(brewer.pal(12, "Paired"),
                  brewer.pal(12, "Set3"),
                  brewer.pal(8, "Dark2"))
  cols_new <- cols_new[-14]
  cols_new <- cols_new[-12]
  cols_new <- cols_new[-9]
  cols_new <- cols_new[-20]
  cols_new <- cols_new[-10]
  mypal <- c(pal_npg("nrc", alpha = 0.7)(10),
             pal_tron("legacy", alpha = 0.7)(7),
             pal_lancet("lanonc", alpha = 0.7)(9),
             pal_simpsons(palette = c("springfield"), alpha = 0.7)(16),
             pal_rickandmorty(palette = c("schwifty"), alpha = 0.7)(6),
             pal_futurama(palette = c("planetexpress"), alpha = 0.7)(5),
             pal_startrek(palette = c("uniform"), alpha = 0.7)(5))
  mycoloursP<- c(cols_new, mypal)
  return(mycoloursP)
}

