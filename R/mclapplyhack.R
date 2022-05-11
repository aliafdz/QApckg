#' @title Execution of parallel:mclapply() on Windows machines
#'
#' @description Mimics forking on Windows machines just like it is done with
#'   \code{\link{mclapply}} in Mac or Linux.
#'
#' @param mc.cores Number of cores to use for parallelization.
#'
#' @note The function code was extracted from \code{post-10-mclapply-hack.R} written by Nathan VanHoudnos, and was only
#'   added in this package for parallelization using Windows operating system. See the original code
#'   in <https://github.com/nathanvan/mcmc-in-irt/blob/master/post-10-mclapply-hack.R>
#'
#' @author Nathan VanHoudnos
#' @import parallel
#' @import ltools
#' @export
#' @seealso \code{\link{mclapply}}

## Define the hack
mclapply.hack <- function(..., mc.cores=NULL) {
  ## Create a cluster
  if( is.null(mc.cores) ) {
    size.of.list <- length(list(...)[[1]])
    mc.cores <- min(size.of.list, detectCores())
  }
  ## N.B. setting outfile to blank redirects output to
  ##      the master console, as is the default with
  ##      mclapply() on Linux / Mac
  cl <- makeCluster( mc.cores, outfile="" )

  ## Find out the names of the loaded packages
  loaded.package.names <- c(
    ## Base packages
    sessionInfo()$basePkgs,
    ## Additional packages
    names( sessionInfo()$otherPkgs ))

  tryCatch( {

    ## Copy over all of the objects within scope to
    ## all clusters.
    this.env <- environment()
    while( identical( this.env, globalenv() ) == FALSE ) {
      clusterExport(cl,
                    ls(all.names=TRUE, env=this.env),
                    envir=this.env)
      this.env <- parent.env(environment())
    }
    clusterExport(cl,
                  ls(all.names=TRUE, env=globalenv()),
                  envir=globalenv())

    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parLapply( cl, 1:length(cl), function(xx){
      lapply(loaded.package.names, function(yy) {
        require(yy , character.only=TRUE)})
    })

    ## Run the lapply in parallel
    return( parLapply( cl, ...) )
  }, finally = {
    ## Stop the cluster
    stopCluster(cl)
  })
}
