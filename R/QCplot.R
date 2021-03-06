#' @title QC plot by read position
#'
#' @description Draws a scatter plot with Phred score values by read position.
#'
#' @note This function is only defined for correct execution of \code{\link{PoolQCbyPos}}
#'  function from the same package.
#'
#' @param fvnm1,fvnm2 Matrix or array containing Phred score values by read position for
#'   0.05, 0.25, 0.5, 0.75 and 0.95 quantiles.
#'   If a QC by position plot for R1 and R2 files is required, both arguments are needed, one
#'   for each file.
#' @param snm Character indicating the name of the pool of evaluated reads. Only required when
#'   both \code{fvnm1} and \code{fvnm2} are provided.
#' @param SW Logical indicating if the plot should include quality profile by SW (Sliding Window).
#'   If \code{TRUE}, Phred Scores are computed as moving averages for windows of 10 base length
#'   slid along the sequence.
#' @param FL Logical indicating if the first argument corresponds to the scores of FLASH extended reads.
#'   If \code{TRUE}, only \code{fvnm1} argument is required.
#'
#' @return A scatter plot with desired values will be shown in the active plots window.
#'
#' @details \code{fvnm1} and \code{fvnm2} arguments are obtanied from the \code{\link{QCscores}}
#'   function with argument \code{byPos=TRUE}.
#'
#' @seealso \code{\link{PoolQCbyPos}}, \code{\link{QCscores}}
#' @examples
#' flashDir <- "./flash"
#' flashfiles <- list.files(flashDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' lst1 <- QCscores(file.path(flashDir,flashfiles[1]),byPos=T)
#' fvnm1 <- lst1$fvnq
#' QCplot(fvnm1,FL=TRUE) # QC plot for FLASH fastq file
#' QCplot(fvnm1,SW=TRUE,FL=TRUE) # QC plot by SW for FLASH fastq file
#' @export
#' @author Alicia Aranda

# Gràfic de QC per posició i SW dels reads originals abans i després de flash
# Les línies dels quantils indiquen que com a màxim el % indicat dels reads tenen aquell Phred score
QCplot <- function(fvnm1,fvnm2,snm,SW=FALSE,FL=FALSE){

  if(missing(fvnm2)&missing(snm)&FL==TRUE){
    # El primer argument seran les dades provinent de FLASH
    fvnm <- fvnm1
    # Avalua si volem realitzar el gràfic de FLASH per sliding window i el calcula
    if(SW){
      nc <- ncol(fvnm)
      fvnm <- sapply(1:(nc-10),function(iwin)
        rowMeans(fvnm[,iwin:(iwin+9),drop=FALSE]))
    }
    # Gràfic de QC per posició i SW dels reads després de filtrar per FLASH
    # Nº de columnes
    nc <- ncol(fvnm)
    # Transforma en llista els resultats de quantils de phred score entre total de reads
    bxp.dt <- list(stats=fvnm,names=1:nc)
    # Dibuixa el gràfic i les línies de cada quantil
    plot(1:nc,fvnm[2,],type="l",col="darkgreen",ylim=c(0,40),xaxt="n",
         xlab="Position",ylab="Phred score")
    lines(1:nc,fvnm[3,],type="l",lwd=2,col="darkgreen")
    lines(1:nc,fvnm[4,],type="l",col="darkgreen")
    lines(1:nc,fvnm[1,],type="l",lty=3,col="darkgreen")
    lines(1:nc,fvnm[5,],type="l",lty=3,col="darkgreen")
    axis(1,at=seq(10,nc,10),lab=seq(10,nc,10),las=2,cex.axis=0.8)
    abline(h=c(20,23,30,33),lty=4,col="gray")
    legend("bottomleft",lty=c(3,1,1,1,3),lwd=c(1,1,2,1,1),cex=0.8,
           legend=rev(c(0.05,0.25,0.5,0.75,0.95)),title="quantiles")
    legend("bottom",lwd=2,col="darkgreen",legend="Flash reads",horiz=TRUE)
  }
  else {
    # Columnes de la taula (fracció dels quantils de phred score entre total de reads)
    if(SW){
      nc <- ncol(fvnm1) # Abans de fer SW
      # Aplica SW: per cada bloc de 10 columnes fa la mitjana del phred score de les files
      fvnm1 <- sapply(1:(nc-10),function(iwin)
        rowMeans(fvnm1[,iwin:(iwin+9),drop=FALSE]))
      # Mateix per R2
      nc <- ncol(fvnm2)
      fvnm2 <- sapply(1:(nc-10),function(iwin)
        rowMeans(fvnm2[,iwin:(iwin+9),drop=FALSE]))
    }
    nc <- ncol(fvnm1)
    # Dibuixa el gràfic amb el phred score per posició amb les línies dels quantils per R1
    plot(1:nc,fvnm1[2,],type="l",col="navy",ylim=c(0,40),xaxt="n",
         xlab="",ylab="Phred score")
    lines(1:nc,fvnm1[3,],type="l",lwd=2,col="navy")
    lines(1:nc,fvnm1[4,],type="l",col="navy")
    lines(1:nc,fvnm1[1,],type="l",lty=3,col="navy")
    lines(1:nc,fvnm1[5,],type="l",lty=3,col="navy")
    axis(1,at=seq(10,nc,10),lab=seq(10,nc,10),las=2,cex.axis=0.8)
    # El mateix per R2
    nc <- ncol(fvnm2)
    lines(1:nc,fvnm2[2,],type="l",col="maroon")
    lines(1:nc,fvnm2[3,],type="l",lwd=2,col="maroon")
    lines(1:nc,fvnm2[4,],type="l",col="maroon")
    lines(1:nc,fvnm2[1,],type="l",lty=3,col="maroon")
    lines(1:nc,fvnm2[5,],type="l",lty=3,col="maroon")
    # Línies grises per alguns valors phred score que puguin ser d'utilitat
    abline(h=c(20,23,30,33),lty=4,col="gray")
    legend("bottomleft",lty=c(3,1,1,1,3),lwd=c(1,1,2,1,1),cex=0.8,
           legend=rev(c(0.05,0.25,0.5,0.75,0.95)),title="quantiles")
    legend("bottom",lwd=2,col=c("navy","maroon"),legend=c("R1","R2"),horiz=TRUE)
    # Títols en funció de si es SW o no
    if(SW)
      title(main=paste("Quality profile by SW (size 10,step 1):",snm))
    else
      title(main=paste("Quality profile by position:",snm))
  }}
