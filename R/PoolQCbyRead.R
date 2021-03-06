#' @title Evaluate QC by read
#'
#' @description This function evaluates fastq files after the execution of the FLASH program to extend
#'   paired-end reads, and returns QC by read plots in pdf format. The results of this function are
#'   important for defining the maximum fraction of bases below Q30 allowed in reads, which will
#'   be used in \code{\link{FiltbyQ30}} function.
#'
#'
#' @param flashfiles Vector including the paths of FLASH processed files, with fastq extension.
#' @param samples Data frame with relevant information to identify the samples of the sequencing experiment, including
#'   \code{Patient.ID, MID, Primer.ID, Region, RefSeq.ID}, and \code{Pool.Nm} columns.
#' @param primers Data frame with information about the \emph{primers} used in the experiment, including
#'   \code{Ampl.Nm, Region, Primer.FW, Primer.RV, FW.pos, RV.pos, FW.tpos, RV.tpos, Aa.ipos},
#'     and \code{Aa.lpos} columns.
#' @param ncores Number of cores to use for parallelization with \code{\link{mclapply.hack}}.
#'
#' @importFrom foreach foreach
#' @import parallel
#' @return After execution a message will appear in console, indicating that the following
#' report files have been generated (and saved in a reports folder):
#' \enumerate{
#'   \item{\code{PoolQCbyRead_PoolName.pdf}: This file is generated for each pool used in the
#'     experiment, after extracting its name from \code{samples} data frame. The pdf includes
#'     includes a representation of bases below Q30 (in nº of reads and percentage) by read.}
#'   \item{\code{PoolReadLengths.pdf}: Includes one plot for each pool representing the read length distribution.}}
#'
#' @export
#' @seealso \code{\link{R1R2toFLASH}}, \code{\link{QCscores}}
#' @examples
#' flashDir <- "./flash"
#' repDir <- "./reports"
#' # Save the file names with complete path
#' flashfiles <- list.files(flashDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' # Get data
#' samples <- read.table("./data/samples.csv", sep="\t", header=T,
#'                      colClasses="character",stringsAsFactors=F)
#' primers <- read.table("./data/primers.csv", sep="\t", header=T,
#'                       stringsAsFactors=F)
#' PoolQCbyRead(flashfiles,samples,primers)
#' @author Alicia Aranda



PoolQCbyRead <- function(flashfiles,samples,primers,ncores=1) {

  # Si la ruta on es troben els fitxers flash no està ben especificada, intenta buscar la
  # ruta correcta a partir del directori de treball
  # Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
  if(length(flashfiles)==0) {
    flashfiles <- list.files(paste(getwd(),"/flash",sep=''))
    if(length(flashfiles)==0) {
      stop("Couldn't find FLASH result files, please indicate correct path.\n")
    }}

  # Si cap dels fitxers indicats a la carpeta flash no existeix, atura l'execució i
  # mostra un missatge d'error
  if(any(!file.exists(flashfiles))){
    stop(paste(flashfiles[!file.exists(flashfiles)],"does not exist.\n"))
  }

  # Si la ruta on es troben els fitxers data no està ben especificada, atura l'execució
  # i mostra un missatge d'error
  if(length(samples)==0|length(primers)==0) {
    stop("Please check data folder files, something is missing.\n")}

  # Si no existeix la carpeta on es guarden els fitxers resultants de la funció,
  # es genera automàticament a la carpeta de treball
  if(!dir.exists("./reports")) {
    dir.create("./reports")}
    repDir <- "./reports"

  ### Fitxers que resulten de Flash
  flnms <- basename(flashfiles)
  # Guarda del fitxer només el nom del pool seguit de S1 o S2
  snms <- sub("_flash\\.fastq$","",flnms)
  # Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
  parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
  if(is.vector(parts))
    parts <- matrix(parts,nrow=1)
  colnames(parts) <- c("PatID","SmplID")

  ### Amb 'mclapply()' aplica la funció 'QCscores()' del paquet sobre tots els fitxers
  # de la carpeta flash per calcular les puntuacions per read (argument byRead)
  flashdata <-mclapply.hack(c(1:length(snms)),function(i){
    QCscores(flashfiles[i],byRead=T)},mc.cores=ncores)
  names(flashdata) <- flnms

######## POOLQCBYREAD

  ### Loop sobre pools
  foreach(i=1:length(snms)) %do% {
    # Guarda els resultats de la funció 'QCscores()' sobre FLASH que corresponen al pool avaluat
    lst1 <- flashdata[[i]]

    # Genera el pdf on aniran els gràfics
    pdf.flnm <- paste("PoolQCbyRead_",parts[i,"PatID"],".pdf",sep="")
    pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=10.5)
    par(mfrow=c(2,1))
    # Després d'aplicar la funció, de la llista resultant agafa la matriu de les bases<Q30, i agafa només
    # aquelles per sota del quantil a 0.99
    nl30 <- lst1$all.nl30[lst1$all.nl30<quantile(lst1$all.nl30,p=0.99)]

    # Histograma amb el subset anterior
    hdt <- hist(nl30,breaks=30,main="",xlab="# bases below Q30")
    # Aplica els 3 quantils indicats a les bases<Q30
    qs <- quantile(lst1$all.nl30,p=c(0.5,0.8,0.95))
    # Dibuixa una fletxa blava en l'eix X per indicar el rang del quantil 0.95
    arrows(x0=0,x1=qs[3],y0=0,code=3,angle=90,len=0.1,col="blue",lwd=2)
    # Dibuixa una fletxa rosa en l'eix X per indicar el rang del quantil 0.80
    arrows(x0=0,x1=qs[2],y0=0,code=3,angle=90,len=0.1,col="maroon",lwd=2)
    # Dibuixa una línia vertical blava indicant el quantil 0.5
    abline(v=qs[1],lty=4,col="blue",lwd=2)
    q95 <- paste("q95% ",qs[3])
    q80 <- paste("q80% ",qs[2])
    med <- paste("q50% ",qs[1])
    # hist$mids= the n cell midpoints; hist$counts= n integers; for each cell, the number of x[] inside.
    text(x=max(hdt$mids)*0.95,y=max(hdt$counts)*0.9,adj=1,cex=0.8,
         paste(q95,q80,med,sep="\n"))
    title(main=parts[i,"PatID"],line=1)

    # Guarda les bases<30 entre el nº de cicles -> fracció de les bases per read
    all.fnl30 <- lst1$all.fnl30

    # D'aquesta matriu de divisió agafa els valors per sota del quantil 0.99
    fnl30 <- all.fnl30[all.fnl30<quantile(all.fnl30,p=0.99)]
    # Nou histograma però amb les fraccions de bases per read
    hdt <- hist(fnl30,breaks=30,main="",
                xlab="fraction of bases below Q30 by read")
    qs <- round(quantile(all.fnl30,p=c(0.5,0.8,0.95)),3)
    arrows(x0=0,x1=qs[3],y0=0,code=3,angle=90,len=0.1,col="blue",lwd=2)
    arrows(x0=0,x1=qs[2],y0=0,code=3,angle=90,len=0.1,col="maroon",lwd=2)
    abline(v=qs[1],lty=4,col="blue",lwd=2)
    q95 <- paste("q95%",qs[3])
    q80 <- paste("q80%",qs[2])
    med <- paste("q50%",qs[1])
    # Valors amb 3 decimals
    vals <- sprintf("%5.3f",c(qs[3],qs[2],qs[1]))
    vals <- paste(c("q95%","q80%","q50%"),vals,collapse="\n")
    text(x=max(hdt$mids)*0.98,y=max(hdt$counts)*0.95,adj=1,cex=0.8,vals)
    # Línies verticals grises en les posicions x indicades
    abline(v=c(0.01,0.02,0.05),lty=4,col="gray")

    # Calcula el sumatori (en %) de reads que disposen de menys de 1, 2 o 5% de les seves bases per sota de Q30
    # Si Pr(f<=0.05)=82%, per exemple, vol dir que si eliminem els reads amb més del 5% de les bases<Q30 encara
    # ens quedem amb un 82% dels reads totals.
    p1pct <- sum(all.fnl30<=0.01)/length(all.fnl30)*100
    p2pct <- sum(all.fnl30<=0.02)/length(all.fnl30)*100
    p5pct <- sum(all.fnl30<=0.05)/length(all.fnl30)*100
    vals <- sprintf("%5.1f",c(p1pct,p2pct,p5pct))
    txt <- paste("Pr(f<=",c(1,2,5),"%)  ",vals,"%",sep="")
    txt <- paste(txt,collapse="\n")
    text(x=max(hdt$mids)*0.98,y=max(hdt$counts)*0.75,adj=1,cex=0.8,txt)
    title(main=parts[i,"PatID"],line=1)
    dev.off()
  }

##### LENPEAKSBYPOOL
  # Genera el fitxer pdf indicat en el directori de reports
  pdf(file.path(repDir,"PoolReadLengths.pdf"),paper="a4",width=6,height=10)
  par(mfrow=c(2,1),mar=c(5,4,4,2.5)+0.1)
  par(mfrow=c(2,1))

  foreach(i=1:length(snms)) %do% {

    # Guarda els resultats de la funció 'QCscores()' sobre FLASH que corresponen al pool avaluat
    lst1 <- flashdata[[i]]

    # Taula per indicar la freqüència de cada nº de cicles extrets en vln
    lnfrq <- table(lst1$all.ln)
    # Defineix els valors del nº de cicles (longitud dels reads)
    x <- as.integer(names(lnfrq))
    # Defineix les freqüències
    y <- as.vector(lnfrq)
    # Gràfic de la longitud dels reads i la seva freq
    plot(x,y,type="h",xlab="Read length",ylab="Frequency")
    title(main=paste(snms[i],"- read lengths"),line=2)
    par(new=T)
    # Calcula la freq acumulada de la longitud dels reads per dibuixar la línea al gràfic
    plot(x,cumsum(y/sum(y)),type="l",ylim=c(0,1),
         axes=F,xlab=NA,ylab=NA,col="blue")
    # L'eix y de la línia de freq és de color blau amb intervals de 0.1
    axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2,cex.axis=0.8)
    grid()
    # Guarda en una variable les longituds de read amb freq acumulada major al 5%
    idx <- which( y/sum(y) >= 0.05 )
    x[idx]
    # Mostra un altre títol de gràfic indicant on es troben els pics de longitud de seq
    title(main=paste("Peaks at",paste(x[idx],collapse=", ")),
          cex.main=0.8,line=0.8)
  }

  dev.off()
  return(cat('The generated files are in the reports folder.\n'))}
