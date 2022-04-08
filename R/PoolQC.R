#' PoolQC
#'
#' @details Files indicated in \code{runfiles} and \code{flashfiles} arguments must be located in different folders,
#'   named runDir and flashDir respectively.
#' @param runfiles Vector including the names of Illumina MiSeq Raw Data files, often with fastq.gz extension
#' @param flahsfiles Vector including the names of FLASH processed files, with fastq extension
#' @param samples Data frame
#' @param primers Data frame
#'
#' @export

#runfiles <- list.files(runDir)
#flashfiles <- list.files(flashDir)

PoolQC <- function(runfiles,flashfiles,samples,primers) {

# Si la ruta on es troben els fitxers run no està ben especificada, intenta buscar la
# ruta correcta a partir del directori de treball
# Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
if(length(runfiles)==0) {
  runfiles <- list.files(paste(getwd(),"/run",sep=''))
  if(length(runfiles)==0) {
    stop("Couldn't find any Raw Data file, please indicate correct path")
  }
}

# Si la ruta on es troben els fitxers flash no està ben especificada, intenta buscar la
# ruta correcta a partir del directori de treball
# Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
if(length(flashfiles)==0) {
  flashfiles <- list.files(paste(getwd(),"/flash",sep=''))
  if(length(flashfiles)==0) {
    stop("Couldn't find FLASH result files, please indicate correct path")
  }
}

# Si la ruta on es troben els fitxers data no està ben especificada, atura l'execució
# i mostra un missatge d'error
if(length(samples)==0||length(primers)==0) {
  stop("Please check data folder files, something is missing")
}

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

### Suma a la longitud de l'amplicó la longitud del MID i M13
pln <- pln + 2*(20+10)

### Fitxers que resulten de Flash
flnms <- flashfiles
# Guarda del fitxer només el nom del pool seguit de S1 o S2
snms <- sub("_flash\\.fastq$","",flashfiles)
# Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

### Fitxers R1 i R2 originals, de la carpeta run
# Guarda el nom dels fitxers R1 i R2 buscant quins contenen el terme desitjar en els noms
# dels fitxers run
R1.flnms <- runfiles[which(grepl('R1',runfiles))]
R2.flnms <- runfiles[which(grepl('R2',runfiles))]

### Sincronitza la longitud dels amplicons dels pools amb el nom d'aquests
pln <- pln[parts[,"PatID"]]

### Amb 'sapply()' aplica la funció 'fn.fastq()' del paquet sobre tots els fitxers
# de la carpeta run
rundata <- sapply(c(1:length(runfiles)),function(i){
  fn.fastq(file.path(runDir,runfiles)[i])})
# Assigna cada resultat al seu nom de fitxer corresponent
colnames(rundata) <- runfiles

### Amb 'sapply()' aplica la funció 'fn.fastq()' del paquet sobre tots els fitxers
# de la carpeta flash
flashdata <- sapply(c(1:length(snms)),function(i){
  fn.fastq(file.path(flashDir,flashfiles[i]),pln[i])}) # pln per guardar com argument ln la longitud de l'amplicó!
colnames(flashdata) <- flashfiles

### Loop sobre pools
for(i in 1:length(snms))
{ # Genera el pdf PoolQCbyPos del pool avaluat
  pdf.flnm <- paste("PoolQCbyPos",parts[i,1],"pdf",sep=".")
  # Genera la ruta on es guardarà el pdf (reports)
  pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10.5,height=6.5)
  par(mfrow=c(2,1),mar=c(3,4,1.5,2)+0.1)


  # Guarda els resultats de la funció 'fn.fastq()' que corresponen al pool avaluat
  # i diferenciant els fitxers R1 i R2
  # Això és pel gràfic d'abans de flash, avalua la qualitat dels reads inicials
  lst1 <- rundata[,which(colnames(rundata) %in% R1.flnms[i])]
  lst2 <- rundata[,which(colnames(rundata) %in% R2.flnms[i])]

  # Guarda la taula amb quantils de phred score per posició entre total de reads
  # per les dades de R1 i R2
  fvnm1 <- lst1$fvnq
  fvnm2 <- lst2$fvnq

  # Genera el gràfic QCbyPos per les dades crues de la carpeta run
  plot.QC(fvnm1,fvnm2,snms[i])

  # Guarda els resultats de la funció 'fn.fastq()' sobre FLASH que corresponen al pool avaluat
  lst <- flashdata[,i]

  # Matriu amb quantils de phred score per posició entre total de reads per les dades FLASH
  fvnm <- lst$fvnq
  # Genera el gràfic QCbyPos dels resultats de FLASH
  plot.QC(fvnm,FL=TRUE)

  # Gràfic SW (sliding window, per regions) pels fitxers originals de run
  plot.QC(fvnm1,fvnm2,snms[i],SW=TRUE)

  # Gràfic SW pels resultats de flash
  plot.QC(fvnm,SW=TRUE,FL=TRUE)

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
  axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2)
  grid()

  dev.off()
}
return(cat('QC by position plots done successfully'))}
