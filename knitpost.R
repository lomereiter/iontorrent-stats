KnitPost <- function(input, base.url = "/") {
    require(knitr)
    opts_knit$set(base.url = base.url)
    fig.path <- paste0("figures/", sub(".Rmd$", "", basename(input)), "/")
    opts_chunk$set(fig.path = fig.path)
    opts_chunk$set(fig.cap = "center")
    opts_chunk$set(fig.width = 10)
    opts_chunk$set(tidy=FALSE) # no thanks, I'll better tidy up my code manually
    opts_chunk$set(highlight=TRUE)
    render_jekyll()
    knit(input, envir = parent.frame())
}
