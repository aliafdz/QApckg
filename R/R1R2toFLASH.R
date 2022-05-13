#' @title Run FLASH to extend paired-end reads and generate report graphs.
#' @author Alicia Aranda
#'
#' @description This function applies \code{\link{executeFLASH}} over R1 and R2 reads
#'   for multiple sample pools and returns a file with the extended reads for each pool.
#'
#'
#' @param runfiles Vector including the paths of Illumina MiSeq Raw Data files, often with fastq.gz extension.
#' @param flash File path of the FLASH executable.
#' @param min.ov Minimum overlap (in nt) between R1 and R2 reads.
#' @param max.ov Maximum overlap (in nt) between R1 and R2 reads.
#' @param err.lv Mismatch fraction accepted in overlapping.
#'
#' @return The function returns a \code{\link{data.frame}} object containing FLASH results
#'   for sequenced regions.
#'
#' @return After the execution, a fastq file with extended reads for each pool will be saved in
#'   a new folder named flash. Additionaly, two files will be saved in a reports folder:
#'   \enumerate{
#'   \item{\code{FLASH_barplot.pdf}: Bar plots representing extended vs not extended reads
#'     and the yield of the process for each pool.}
#'   \item{\code{FLASH_report.txt}: Includes the data returned by the function with used FLASH parameters.}
#' }
#' @import grDevices
#' @import graphics
#' @import utils
#' @import foreach
#' @import ShortRead
#' @import Biostrings
#' @importFrom RColorBrewer brewer.pal
#' @seealso \code{\link{executeFLASH}}
#' @export
#' @examples
#' runDir <- "./run"
#' flash <- "./FLASH/flash.exe"
#' # Save the file names with complete path
#' runfiles<-list.files(runDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' min.ov <- 20
#' max.ov <- 300
#' err.lv <- 0.1
#' flashres <- R1R2toFLASH(runfiles,flash,min.ov,max.ov,err.lv)

R1R2toFLASH <- function(runfiles,flash,min.ov=20,max.ov=300,err.lv=0.10)
{
# Si la ruta on es troben els fitxers run no està ben especificada, intenta buscar la
# ruta correcta a partir del directori de treball
  # Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
  if(length(runfiles)==0) {
    runfiles <- list.files(paste(getwd(),"/run",sep=''),recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
    if(length(runfiles)==0) {
      stop("Couldn't find any Raw Data file, please indicate the correct path.\n")
    }
  }
  # Si cap dels fitxers indicats a la carpeta run no existeix, atura l'execució i
  # mostra un missatge d'error
  if(any(!file.exists(runfiles))){
    stop(paste(runfiles[!file.exists(runfiles)],"does not exist.\n"))
  }

  # Afegim un test per comprovar que l'executable FLASH existeix
  if(!file.exists(flash)){
    stop("Couldn't find the FLASH executable, please indicate the correct path.\n")
  }

  # Si no existeixen les carpetes on es guarden els fitxers resultants de la funció,
  # es generen automàticament a la carpeta de treball
  if(!dir.exists("./flash")) {
    dir.create("./flash")}
    flashDir <- "./flash"

  if(!dir.exists("./reports")) {
  dir.create("./reports")}
  repDir <- "./reports"

# Guarda amb la funció 'basename()' el nom dels arxius de la carpeta run, sense la ruta completa
flnms <- basename(runfiles)

# La funció sub() permet substituir un patró pel que indiquem com 2n argument
# En aquest cas, les variables runfiles i snms son idèntiques (de moment)
snms <- sub("\\.fastq$","",flnms)

# Taula on es disposa els nom de cada fitxer fastq en la primera columna i després es
# separa cadascun del termes del nom del fitxer en successives columnes
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
# Elimina les columnes 3 i 5 de la taula (que indiquen L001 i l'extensió)
parts <- parts[,-c(3,5)]
# Assigna els noms a les columnes: ID pacient, ID mostra/pool, read
colnames(parts) <- c("PatID","SmplID","Read")

# Dels 4 fitxers que hi havia, guarda en 2 variables els que corresponen a R1 i R2
R1.flnms <- runfiles[parts[,3]=="R1"] # parts[,3] = columna 3 (read) de la taula parts
R2.flnms <- runfiles[parts[,3]=="R2"]

# Guarda dos noms de fitxer corresponents a R1: pool, ID de mostra i extensió flash.fastq
# Com hi ha 2 fitxers amb R1, vol dir que avaluem 2 regions o pools
out.flnms <- paste(parts[parts[,3]=="R1",1],parts[parts[,3]=="R1",2],
                   "flash.fastq",sep="_")

# Genera la ruta dels fitxers definits abans a la carpeta flash
out.flnms <- file.path(flashDir,out.flnms)

# Guarda de la taula només els R1
parts <- parts[parts[,3]=="R1",,drop=FALSE]

# Defineix les opcions necessàries per a l'execució de FLASH
# Aquesta variable es guarda com a global, per tal de poder accedir a ella des de la funció
# 'executeFLASH()' i que no aparegui error
flash.opts <- paste("-m",min.ov,"-M",max.ov,"-x",err.lv)

## Itera sobre el total de pools (nº de fitxers que es generaran) i aplica la funció 'executeFLASH()'
# del paquet, que permet realitzar l'extensió dels reads R1 i R2 i guardar el nº de reads units (extended)
# i no units (no extended)
# La funció 'foreach()' funciona com un bucle for però de forma més ràpida, i permet executar la funció
# definida a partir de l'argument .export
# Parteix dels fitxers R1 i R2 de cadascun dels pools, així com el nom del fitxer fastq
# resultant que es guardarà a la carpeta flash
flashres <- foreach(i=1:length(out.flnms),.export='executeFLASH',.packages='ShortRead') %do%
  executeFLASH(R1.flnms[i],R2.flnms[i],flash,flash.opts,out.flnms[i])

# Construeix una matriu 2x2 que inclogui els resultats d'aplicar FLASH per cada pool
res <- matrix(unlist(flashres),nrow=length(out.flnms),ncol=2,byrow=TRUE) # length(out.flnms)= 2 en aquest cas, que son els pools
# Assigna com a nom de fila el pool (regió de VHB) i com a columnes els reads segons si s'ha donat o no extensió
rownames(res) <- paste(parts[,1],parts[,2],sep="_")
colnames(res) <- c("Extended","NoExtd")

# Guarda un dataframe amb les dades resultants del FLASH (taula res), afegint la
# columna Yield calculada dividint els reads extended entre el total *100
df.res <- data.frame(res,Yield=round(res[,1]/(res[,1]+res[,2])*100,1))

# Guarda el fitxer de report de FLASH en format .txt
txt.flnm <- file.path(repDir,"FLASH_report.txt")
# Comandes per omplir el fitxer txt generat. Recull els paràmetres indicats al FLASH
# i la taula df.res que conté els resultats obtinguts
sink(txt.flnm)
cat("\nExtending Illumina reads by FLASH\n")
cat("\nFLASH parameters:")
cat("\n    Minimum overlap:",min.ov)
cat("\n    Maximum overlap:",max.ov)
cat("\n        Error level:",err.lv,"\n\n")
# 'capture.output()' permet mostrar la taula de resultats tal i com es mostra en l'output de R
capture.output(df.res,file=txt.flnm,append=TRUE)
sink() # Tanca el fitxer

# Genera el pdf que contindrà el gràfic barplot dels resultats FLASH
pdf.flnm <- file.path(repDir,"FLASH_barplot.pdf")
pdf(pdf.flnm,paper="a4",width=6,height=10)
par(mfrow=c(2,1))

pal=brewer.pal(8,"Dark2") # Crea la paleta de colors
M <- data.matrix(df.res[,1:2]) # Guarda les columnes 1 i 2 de la taula de resultats
ymx <- max(M)*1.15 # Defineix límit superior del gràfic
barplot(t(M),beside=TRUE,las=2,col=pal[1:2],ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal[1:2],cex=0.8,
       legend=c("Extended","Not extended"))
title(main="FLASH results on paired-ends",line=2.5)
title(main="Yield in number of reads",cex.main=1,line=1)

# Genera un altre barplot on s'inclouen els resultats del Yield de FLASH
bp <- barplot(df.res$Yield,col="Lavender",ylim=c(0,100),ylab="Yield",
              names.arg=rownames(df.res),las=2)
text(bp,10,paste(df.res$Yield,"%",sep=""),srt=90,col="navy")
title(main="Yield in percentage",cex.main=1,line=1)

dev.off()

return(result <- df.res)
}
