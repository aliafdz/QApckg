#' R1R2toFLASH
#'
#' A function that runs FLASH program to extend paired-end reads and generates some report graphics.
#'
#' @details Files indicated in \code{runfiles} must be located in a folder named runDir.
#'   Also, a reports folder must be created in the project environment, whose path will be
#'   named as repDir.
#' @param runfiles Character indicating which files are going to be processed, often with fastq.gz extension
#' @param flash Folder path containing FLASH executable
#' @param min.len Minimum length to consider a sequence
#' @param min.ov Minimum overlap (in nt) between R1 and R2
#' @param max.ov Maximum overlap (in nt) between R1 and R2
#' @param err.lv Mismatch fraction accepted in overlapping
#' @param chunck.sz Chunck size to be used by \code{\link{FastqStreamer}} function
#' @return The function returns a \code{\link{data.frame}} object containing FLASH results
#'   for your sequenced regions, but also two report files:
#'   \enumerate{
#'   \item \code{FLASH_barplot.pdf}: Bar plots representing extended vs not extended reads
#'     and the yield of the process for each pool.
#'   \item \code{FLASH_report.txt}: Includes the data returned by the function with FLASH parameters used.
#' }
#' @import grDevices
#' @import graphics
#' @import utils
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' runDir <- "C:/run"
#' flash <- "C:/FLASH/flash.exe"
#' runfiles <- list.files(runDir)
#' R1R2toFLASH(runfiles,flash,min.len=200,min.ov=20,max.ov=300,err.lv=0.10,chunck.sz=1.e6)

R1R2toFLASH <- function(runfiles,flash,min.len=200,min.ov=20,max.ov=300,err.lv=0.10,chunck.sz=1.e6)
{
# Si la ruta on es troben els fitxers run no està ben especificada, intenta buscar la
# ruta correcta a partir del directori de treball
  # Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
  if(length(runfiles)==0) {
    runfiles <- list.files(paste(getwd(),"/run",sep=''))
    if(length(runfiles)==0) {
      stop("Couldn't find any Raw Data file, please indicate correct path")
    }
  }

  # Definim la variable per a la iteració dels chuncks com a global, per tal
  # de poder accedir a ella des de funcions posteriors del pipeline sense
  # requerir la seva definició múltiples cops
  chunck.sz <<- chunck.sz

# La funció sub() permet substituir un patró pel que indiquem com 2n argument
# En aquest cas, les variables runfiles i snms son idèntiques (de moment)
snms <- sub("\\.fastq$","",runfiles)

# Taula on es disposa els nom de cada fitxer fastq en la primera columna i després es
# separa cadascun del termes del nom del fitxer en successives columnes
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
# Elimina les columnes 3 i 5 de la taula (que indiquen L001 i l'extensió)
parts <- parts[,-c(3,5)]
# Assigna els noms a les columnes: ID pacient, ID mostra/pool, read
colnames(parts) <- c("PatID","SmplID","Read")

# Guarda a la variable la ruta dels fitxers que es troben a la carpeta run
flnms <- file.path(runDir,runfiles)

# Dels 4 fitxers que hi havia, guarda en 2 variables els que corresponen a R1 i R2
R1.flnms <- flnms[parts[,3]=="R1"] # parts[,3] = columna 3 (read) de la taula parts
R2.flnms <- flnms[parts[,3]=="R2"]

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
flash.opts <<- paste("-m",min.ov,"-M",max.ov,"-x",err.lv)

## Itera sobre el total de pools (nº de fitxers que es generaran) i aplica la funció 'executeFLASH()'
# del paquet, que permet realitzar l'extensió dels reads R1 i R2 i guardar el nº de reads units (extended)
# i no units (no extended)
# La funció 'foreach()' funciona com un bucle for però de forma més ràpida, i permet executar la funció
# definida a partir de l'argument .export
# Parteix dels fitxers R1 i R2 de cadascun dels pools, així com el nom del fitxer fastq
# resultant que es guardarà a la carpeta flash
flashres <- foreach(i=1:length(out.flnms),.export='executeFLASH',.packages='ShortRead') %do%
  executeFLASH(R1.flnms[i],R2.flnms[i],out.flnms[i],chunck.sz)

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
#cat(capture.output(df.res), sep = '\n')
capture.output(flashres,file=txt.flnm,append=TRUE)
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

return(df.res)
}
