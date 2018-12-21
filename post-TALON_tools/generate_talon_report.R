main <- function() {

    load_packages()



}


load_packages <- function() {

    message("Loading required packages...")

    tryCatch({ suppressPackageStartupMessages(library("DBI"))
               suppressPackageStartupMessages(library("ggplot2"))
               suppressPackageStartupMessages(library("knitr"))
               suppressPackageStartupMessages(library("optparse"))
               suppressPackageStartupMessages(library("RSQLite"))
               suppressPackageStartupMessages(library("readr"))
               suppressPackageStartupMessages(library("rmarkdown"))
               suppressPackageStartupMessages(library("tinytex"))
 
              }, error = function(cond) {
                   message("One of the required packages is missing.")
                   message("Here's the original error message:")
                   message(cond)
              })

    return
}

main()
