#' Create new directory if not already present
#'
#' @param dir_path path where directory is being created
#'
#' @return nothing
#' @export
#'
#' @examples
#' new_dir()
new_dir <- function(dir_path = paste0(getwd(), "outs")){
  if (dir.exists(dir_path) == FALSE) {
    dir.create(dir_path, recursive = TRUE)
    print("New directory created.")
  }else{
    print("Directory already exists")
  }
}
