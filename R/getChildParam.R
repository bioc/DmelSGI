getChildParam <- function(name, default = c()) {
  if (length(name) > 1) {
    warning("Only first element of name will be used.")
    name = name[1]
  }
  
  param = default
  if (exists("opts_chunk")) {
    op = opts_chunk$get()
    if (name %in% names(op)) {
      param = op[[name]]
    }
  }
  param
}

getBaseDir <- function(default = ".") {
  basedir = default
  if (exists("opts_knit")) {
    op = opts_knit$get()
    if ("output.dir" %in% names(op)) {
      basedir = op[["output.dir"]]
    }
  }
  basedir
}

