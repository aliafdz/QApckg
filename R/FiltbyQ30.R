#' @title Filter haplotypes by Q30
#'
#' @description This function applies \code{\link{FastqStreamer}} over a fastq file and removes
#'   all reads that have a defined fraction of bases below Q30. The remaining reads will be saved
#'   in a new fastq file.
#'
#' @param max.pct The maximum percentage of bases below Q30 allowed in reads (by default,5\%).
#' @param flashfiles Vector including the paths of files that are going to be processed,
#'   with fastq extension.
#' @param flashres Table of results obtained after the execution of \code{\link{R1R2toFLASH}}
#'   function.
#'
#' @details This function is designed to be applied after \code{\link{R1R2toFLASH}} function from
#'   the same package. If \code{flashres} is not specified but FLASH extension was previously
#'   done, the function will try to load the FLASH results table from the reports folder.
#'
#' @return A \code{\link{data.frame}} object containing FLASH and Filtering results.
#' @return After the execution, a fastq file with remaining reads for each pool will be saved in
#'   a new flashFilt folder (if it is not previously created). Additionaly, two report files will be
#'   generated in a reports folder:
#'   \enumerate{
#'   \item{\code{FiltQ30.barplot.pdf}: Includes a first Bar plot representing raw reads,
#'     extended reads by FLASH and filtered reads, and a second Bar plot with
#'     the yield by process for each pool.}
#'   \item{\code{FiltQ30_report.txt}: Includes the same data returned by the function.}
#' }
#'   The results table obtained includes two new columns with respect to FLASH results
#'     table, named FiltQ30 (number of filtered reads) and Raw (total sequencing reads).
#'
#' @importFrom RColorBrewer brewer.pal
#' @import ShortRead
#' @seealso \code{\link{R1R2toFLASH}}, \code{\link{FastqStreamer}}
#' @export
#' @examples
#' runDir <- "./run"
#' runfiles <- list.files(runDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' flash <- "./FLASH/flash.exe"
#' flashDir <- "./flash"
#' flashfiles <- list.files(flashDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' flashres <- R1R2toFLASH(runfiles,flash)
#' filtres <- FiltbyQ30(max.pct=0.05,flashfiles,flashres)
#' @author Alicia Aranda

FiltbyQ30 <- function(max.pct=0.05,flashfiles,flashres){
  # La taula resulta d'aplicar la funció R1R2toFLASH
  # Si no existeix la variable o no s'ha incorporat, llegeix la taula del fitxer
  # .txt que podria estar guardat a la carpeta de reports
  if(!exists('flashres')|missing(flashres)|length(flashres)==0){
    tryCatch(
      expr = {
        flashres <- read.table("./reports/FLASH_report.txt",skip=8)},
      # Si no pot trobar la taula resultant de FLASH, atura l'execució i genera un
      # missatge d'error
      error = function(e){
        message("Error: Check if you have executed the previous functions needed.\n")
        stop(print(e))
      })
  }
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

  # Si no existeixen les carpetes on es guarden els fitxers resultants de la funció,
  # es generen automàticament a la carpeta de treball
  if(!dir.exists("./flashFilt")) {
    dir.create("./flashFilt")}
    flashFiltDir <- "./flashFilt"

  if(!dir.exists("./reports")) {
    dir.create("./reports")}
    repDir <- "./reports"

### Fitxers que resulten de Flash
# Guarda del fitxer només el nom del pool seguit de S1 o S2
snms <- sub("_flash\\.fastq$","",basename(flashfiles))

# Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

# Genera dos fitxers amb extensió fastq (1 per pool) que resultaran d'aquest script
oflnms <- paste(parts[,"PatID"],"flashFilt.fastq",sep="_")
# Indica la ruta dels fitxers output en la carpeta flashFilt
oflnms <- file.path(flashFiltDir,oflnms)

### Loop sobre pools
# Vector que inclou tants 0s com nº de pools
freads <- integer(length(snms))

for(i in 1:length(snms))
{
  # Si existeix el fitxers flashFilt a la carpeta, l'elimina (ja que s'ha de generar ara)
  if(file.exists(oflnms[i]))
    file.remove(oflnms[i])

  ### Aplica streamer (iteració) en el fitxer fastq avaluat
  strm <- FastqStreamer(flashfiles[i])

  # Aquesta variable és per la funció 'writeFastq()' i permet generar el fitxer fastq
  # de nou
  mode <- "w"
  # Defineix les variables inicialment a 0
  in.rds <- filt.rds <- 0
  ### Carrega el fitxer fastq per chuncks
  while(length(sqq <- yield(strm)))
  { ### Actualitza el nº de reads a cada iteració
    # Aquest nº ha de ser igual al nº de reads units per FLASH del pool avaluat
    in.rds <- in.rds+length(sqq)
    ### Calcula els Phred scores (Codi ASCII-33) amb la funció 'quality()'
    phrsc <- as(quality(sqq),"matrix")
    ### A la matriu amb els Phred scores aplica el sumatori de les bases amb puntuació menor a 30 (per sota de Q30)
    nl30 <- apply(phrsc,1,function(x) sum(x<30,na.rm=TRUE))

    ### Divideix les bases<30 entre el nº de cicles (longitud) -> fracció de les bases per read
    fnl30 <- nl30/width(sqq)

    ### Aplica el filtre per eliminar els reads amb més del 5% de les bases<Q30
    sqq <- sqq[fnl30<=max.pct]
    # Actualitza el nº de reads després de filtrar per Q30
    filt.rds <- filt.rds+length(sqq)

    # Genera el fitxer amb extensió fastq que es guardarà a flashFilt
    writeFastq(sqq,oflnms[i],mode,compress=TRUE)
    # Canvia la variable mode perquè a la propera iteració es guardin les noves
    # seqüències filtrades en el mateix fastq sense generar un de nou
    mode <- "a"
  }
  close(strm)
  # Guarda per cada pool el nº de reads després del filtrat
  freads[i] <- filt.rds
}

# Afegeix al fitxer que tenia de FLASH els valors de nº de reads després de filtrat Q30
flashres$FiltQ30 <- freads
# També afegeix al fitxer el total de reads (resultat de sumar els que van fer extensió a FLASH i els que no)
flashres$Raw <- flashres$Extended+flashres$NoExtd

# Guarda els resultats de la taula actualitzada en un fitxer .txt
sink(file.path(repDir,"FiltQ30_report.txt"))
print(flashres)
sink()

### PROVA per veure si es podrien juntar els fitxers de FLASH i de filtrat
# Per mantenir els resultats del pipeline original, de moment aquest codi no s'executa.
#file.rename(file.path(repDir,"FLASH_report.txt"),file.path(repDir,"FLASH+FiltQ30_report.txt"))
#sink(file.path(repDir,"FLASH+FiltQ30_report.txt"),append=TRUE)
#cat("\nFiltering haplotypes by Q30:")
#cat("\n    Accept max of ", max.pct*100, "% bases below Q30 by read\n",sep='')
#capture.output(flashres,file=file.path(repDir,"FLASH+FiltQ30_report.txt"),append=TRUE)
#sink()


# Genera una matriu amb les dades de la taula de resultats: total de reads inicials,
# extensió per FLASH i filtrats per Q30
res <- data.matrix(flashres[,c("Raw","Extended","FiltQ30")])

# Genera el fitxer pdf on es generaran els gràfics
pdf(file.path(repDir,"FiltQ30.barplot.pdf"),paper="a4",width=6.6,height=10)
par(mfrow=c(2,1))

pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
# Defineix el límit de l'eix y a partir del valor màxim de la taula res
ymx <- max(res)*1.15
# Genera un gràfic de barres amb la trasposada de la taula res -> nº de reads per pas
bp <- barplot(t(res),col=pal2[1:3],border=pal1[1:3],beside=TRUE,
              xaxt="n",ylim=c(0,ymx))
axis(1,at=colMeans(bp),rownames(res),cex.axis=0.8,las=2)
legend("top",horiz=TRUE,fill=pal2,cex=0.8,legend=colnames(res))
# Calcula el % de reads respecte el total per veure el rendiment
resy <- round(res/res[,1]*100,1)
ymx <- 115
# Gràfic de barres que mostra el % de reads per pas de filtrat
bp <- barplot(t(resy),col=pal2[1:3],border=pal1[1:3],beside=TRUE,
              xaxt="n",ylim=c(0,ymx))
axis(1,at=colMeans(bp),rownames(res),cex.axis=0.8,las=2)
legend("top",horiz=TRUE,fill=pal2,cex=0.8,legend=colnames(res))
# Afegeix l'etiqueta dels % calculats al gràfic
y <- min(resy)/2
text(x=as.vector(bp),y=y,lab=as.vector(t(resy)),cex=0.5,font=2,srt=90,
     col=pal1[1:3])
dev.off()

return(flashres)
}

