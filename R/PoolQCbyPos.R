#' @title Evaluate QC by position
#'
#' @description This function evaluates fastq files before and after the execution of FLASH program to extend
#'   paired-end reads, and returns QC by position plots in pdf format.
#' @description It can be applied also after filtering FLASH fastq files by Phred Score.
#'
#'
#' @param runfiles Vector including the paths of Illumina MiSeq Raw Data files, often with fastq.gz extension.
#'   If the function is applied for filtered fastq files, this argument must be NA or missing.
#' @param flashfiles Vector including the paths of FLASH processed/filtered files, with fastq extension.
#' @param samples Data frame with relevant information about the samples of the sequencing experiment, including
#'   \code{Patient.ID, MID, Primer.ID, Region, RefSeq.ID}, and \code{Pool.Nm} columns.
#' @param primers Data frame with information about the \emph{primers} used in the experiment, including
#'   \code{Ampl.Nm, Region, Primer.FW, Primer.RV, FW.pos, RV.pos, FW.tpos, RV.tpos, Aa.ipos},
#'     and \code{Aa.lpos} columns.
#' @param ncores Number of cores to use for parallelization with \code{\link{mclapply.hack}}.
#'
#' @return After execution, a pdf file for each pool used in the experiment will be saved in a
#'   reports folder (if it is not previously defined, the function will create this folder),
#'   and a message indicating that the files are generated will appear in console.
#'
#' @export
#' @import parallel
#' @examples
#' runDir <- "./run"
#' flashDir <- "./flash"
#' repDir <- "./reports"
#' # Save the file names with complete path
#' runfiles <- list.files(runDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' flashfiles <- list.files(flashDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' # Get data
#' samples <- read.table("./data/samples.csv", sep="\t", header=T,
#'                      colClasses="character",stringsAsFactors=F)
#' primers <- read.table("./data/primers.csv", sep="\t", header=T,
#'                       stringsAsFactors=F)
#' PoolQCbyPos(flashfiles,samples,primers,runfiles)
#' @author Alicia Aranda


PoolQCbyPos <- function(flashfiles,samples,primers,runfiles,ncores=1) {

# Si aquesta funció s'aplica sobre les dades filtrades per Q30, és a dir no
# s'inclou l'argument runfiles, aquests es defineixen com NA
# D'aquesta manera es pot avaluar d'on venen les dades
if(missing(runfiles)){
  runfiles <- NA
}

# Si s'inclouen els arxius de run però la ruta on es troben no està ben especificada,
# intenta buscar la ruta correcta a partir del directori de treball
# Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
if(!is.na(runfiles[1])&length(runfiles)==0) {
  runfiles <- list.files(paste(getwd(),"/run",sep=''),recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
  if(length(runfiles)==0) {
    stop("Couldn't find any Raw Data files, please indicate correct path.\n")
  }
}
  # Si cap dels fitxers indicats a la carpeta run no existeix, atura l'execució i
  # mostra un missatge d'error
  if(!is.na(runfiles[1])){
    if(any(!file.exists(runfiles))){
    stop(paste(runfiles[!file.exists(runfiles)],"does not exist.\n"))
  }}

# Si la ruta on es troben els fitxers flash no està ben especificada, intenta buscar la
# ruta correcta a partir del directori de treball
# Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
if(length(flashfiles)==0) {
  flashfiles <- list.files(paste(getwd(),"/flash",sep=''),recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
  if(length(flashfiles)==0) {
    stop("Couldn't find FLASH result files, please indicate correct path.\n")
  }
}
  # Si cap dels fitxers indicats a la carpeta flash no existeix, atura l'execució i
  # mostra un missatge d'error
  if(any(!file.exists(flashfiles))){
    stop(paste(flashfiles[!file.exists(flashfiles)],"does not exist.\n"))
  }

# Si la ruta on es troben els fitxers data no està ben especificada, atura l'execució
# i mostra un missatge d'error
if(length(samples)==0|length(primers)==0) {
  stop("Please check data folder files, something is missing.\n")
}

  # Si no existeix la carpeta on es guarden els fitxers resultants de la funció,
  # es genera automàticament a la carpeta de treball
  if(!dir.exists("./reports")) {
    dir.create("./reports")}
    repDir <- "./reports"


###  Guarda la llista de pools en el fitxer samples
pools <- unique(samples$Pool.Nm)

# Amb la funció 'sapply()' aplica, sobre els pools que tenim, una funció que busca
# a quins primers correspon el pool i calcula la longitud màxima de l'amplicó (restant les
# posicions dels primers)
pln <- sapply(pools,function(p)
{ # Guarda els items de la taula de mostres que corresponen al pool avaluat
  idx <- which(samples$Pool.Nm==p)
  # Si no hi ha cap item que correspongui retorna 0
  if(length(idx)==0) return(0)
  # samples$Primer.ID[idx] indica les posicions de l'amplicó del pool avaluat
  # Indica a quin item de la taula de primers corresponen les posicions del pool avaluat
  idx <- which(primers$Ampl.Nm %in% samples$Primer.ID[idx])
  # Calcula el màxim de la resta entre la posició del primer reverse de l'amplicó avaluat i la posició del forward
  max(primers$RV.pos[idx]-primers$FW.pos[idx]+1)})
# Guarda el nom dels fitxers de flash sense la ruta completa
ffiles <- basename(flashfiles)

if(!is.na(runfiles[1])){
### Fitxers que resulten de Flash
# Guarda del fitxer només el nom del pool seguit de S1 o S2
snms <- sub("_flash\\.fastq$","",ffiles)

# Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

### Suma a la longitud de l'amplicó la longitud del MID i M13
pln <- pln + 2*(20+10)

### Fitxers R1 i R2 originals, de la carpeta run
# Guarda el nom dels fitxers R1 i R2 buscant quins contenen el terme desitjar en els noms
# dels fitxers run
rfiles <- basename(runfiles)
R1.flnms <- rfiles[which(grepl('R1',runfiles))]
R2.flnms <- rfiles[which(grepl('R2',runfiles))]

### Sincronitza la longitud dels amplicons dels pools amb el nom d'aquests
pln <- pln[parts[,"PatID"]]

### Amb 'mclapply()' per Windows aplica la funció 'QCscores()' del paquet sobre tots els fitxers
# de la carpeta run per calcular les puntuacions per posició (argument byPos)
rundata <-mclapply.hack(c(1:length(rfiles)),function(i){
  QCscores(runfiles[i],byPos=T)},mc.cores=ncores)
# Assigna cada resultat al seu nom de fitxer corresponent
names(rundata) <- rfiles

### Amb 'mclapply()' aplica la funció 'QCscores()' del paquet sobre tots els fitxers
# de la carpeta flash per calcular les puntuacions per posició (argument byPos)
flashdata <-mclapply.hack(c(1:length(snms)),function(i){
  QCscores(flashfiles[i],ln=pln[i],byPos=T)},mc.cores=ncores) # pln per guardar com argument ln la longitud de l'amplicó!
names(flashdata) <- ffiles
}
else{ # Si no hi ha fitxers run perquè estem aplicant sobre FLASH filtrat per Q30
  ### Fitxers que resulten de Flash després de filtrar per Q30
  # Guarda el nom dels pools treient el sufix dels fitxers de flashFiltDir
  snms <- sub("_flashFilt\\.fastq$","",ffiles)
  # Sincronitza la longitud dels amplicons dels pools amb el nom d'aquests
  pln <- pln[snms]
  # Aplica la funció QCscores per calcular la qualitat per posició i fer els gràfics corresponents
  flashdata <-mclapply.hack(c(1:length(snms)),function(i){
    QCscores(flashfiles[i],ln=pln[i],byPos=T)},mc.cores=ncores) # pln per guardar com argument ln la longitud de l'amplicó!
  names(flashdata) <- ffiles
}

### Loop sobre pools
for(i in 1:length(snms))
{ if(!is.na(runfiles[1])){
  # Genera el pdf PoolQCbyPos del pool avaluat
  pdf.flnm <- paste("PoolQCbyPos",parts[i,1],"pdf",sep=".")
  # Genera la ruta on es guardarà el pdf (reports)
  pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10.5,height=6.5)
  par(mfrow=c(2,1),mar=c(3,4,1.5,2)+0.1)

  # Guarda els resultats de la funció 'QCscores()' que corresponen al pool avaluat
  # i diferenciant els fitxers R1 i R2
  # Això és pel gràfic d'abans de flash, avalua la qualitat dels reads inicials
  lst1 <- rundata[[which(names(rundata) %in% R1.flnms[i])]]
  lst2 <- rundata[[which(names(rundata) %in% R2.flnms[i])]]

  # Guarda la taula amb quantils de phred score per posició entre total de reads
  # per les dades de R1 i R2
  fvnm1 <- lst1$fvnq
  fvnm2 <- lst2$fvnq

  # Genera el gràfic QCbyPos per les dades crues de la carpeta run
  QCplot(fvnm1,fvnm2,snms[i])

  # Guarda els resultats de la funció 'QCscores()' sobre FLASH que corresponen al pool avaluat
  lst <- flashdata[[i]]

  # Matriu amb quantils de phred score per posició entre total de reads per les dades FLASH
  fvnm <- lst$fvnq
  # Genera el gràfic QCbyPos dels resultats de FLASH
  QCplot(fvnm,FL=TRUE)

  # Gràfic SW (sliding window, per regions) pels fitxers originals de run
  QCplot(fvnm1,fvnm2,snms[i],SW=TRUE)

  # Gràfic SW pels resultats de flash
  QCplot(fvnm,SW=TRUE,FL=TRUE)

  # Gràfic de la longitud dels reads abans i després de flash
  par(mfrow=c(1,2))
  # Guarda la llista dels quantils de longitud entre el total de reads per R1,R2 i flash
  stats=cbind(lst1$fvnl,lst2$fvnl,lst$fvnl)
  colnames(stats) <- c("R1","R2","Flash")
  # 'bxp()' serveix per realitzar boxplots basats en sumatoris
  bxp.dt <- list(stats=stats,names=colnames(stats))
  bxp(bxp.dt,pars=list(boxfilll="lavender"),border="navy",
      ylab="Read length",las=2,ylim=c(0,max(stats)))
  title(main="Read length distributions")
}
  else if(is.na(runfiles)){
    # Genera el pdf indicat pel pool avaluat
    pdf.flnm <- paste("PoolFiltQCbyPos",snms[i],"pdf",sep=".")
    # Genera la ruta on es guardarà el pdf (reports)
    pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10.5,height=6.5)
    par(mfrow=c(2,1),mar=c(3,4,1.5,2)+0.1)

    # Guarda els resultats de la funció 'QCscores()' sobre FLASH que corresponen al pool avaluat
    lst <- flashdata[[i]]

    # Matriu amb quantils de phred score per posició entre total de reads per les dades FLASH
    fvnm <- lst$fvnq

    # Genera el gràfic QCbyPos dels resultats de FLASH
    QCplot(fvnm,FL=TRUE)

    # Gràfic SW pels resultats de flash
    QCplot(fvnm,SW=TRUE,FL=TRUE)

    # Canvia la disposició pel gràfic de longitud dels reads
    par(mfrow=c(1,2))}

  # Gràfic de longitud de reads després de flash
  # Guarda la longitud dels reads i la seva freqüència (calculada amb 'table()')
  lnfrq <- table(lst$all.ln)
  # Guarda les longituds
  x <- as.integer(names(lnfrq))
  # Guarda les freqüències
  y <- as.vector(lnfrq)
  # Gràfic de tipus histograma amb línies verticals
  plot(x,y,type="h",xlab="Read length",ylab="Frequency")
  title(main="Read lengths (Flash reads)")
  par(new=T)
  # Càlcul de la línia de freqüència acumulada per representar-la al gràfic
  plot(x,cumsum(y/sum(y)),type="l",ylim=c(0,1),
       axes=F,xlab=NA,ylab=NA,col="blue")
  # L'eix y de la línia de freq és de color blau amb intervals de 0.1
  axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2)
  grid()

  dev.off()
}
return(cat('The generated files are in the reports folder.\n'))}
