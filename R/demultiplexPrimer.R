#' @title Trim specific primer sequences
#'
#' @description Demultiplex reads by identifying primer sequences within windows of expected positions in the sequenced reads.
#' It is important to note that MID and primer sequences will be trimmed from reads after the identification of primers,
#'   but amplicon length is not predetermined.
#'
#' @details After demultiplexing reads by MID with \code{\link{demultiplexMID}} function, primer sequences are identified
#'   in both strands. First, forward strands are recognized by searching FW primer sequence in 5' end and the
#'   reverse complement of RV primer sequence in 3' end. Then, reverse strands are recognized by searching RV
#'   primer sequence in 5' end and FW primer sequence in 3' end, obtaining the reverse complement of all reads
#'   identified as reverse strands. So, both strands are obtained in a way that facilitates their intersection.
#'
#'
#' @param splitfiles Vector including the paths of demultiplexed files by MID, with fna extension.
#' @param samples Data frame with relevant information about the samples of the sequencing experiment, including
#'   \code{Patient.ID, MID, Primer.ID, Region, RefSeq.ID}, and \code{Pool.Nm} columns.
#' @param primers Data frame with information about the \emph{primers} used in the experiment, including
#'   \code{Ampl.Nm, Region, Primer.FW, Primer.RV, FW.pos, RV.pos, FW.tpos, RV.tpos, Aa.ipos},
#'     and \code{Aa.lpos} columns.
#' @param prmm Number of mismatches allowed between the primers and read sequences.
#' @param min.len Minimum length desired for haplotypes. Any sequence below this length will be discarted.
#' @param target.st,target.end Initial and end positions between which primer sequences will be searched.
#'
#' @return A list containing the following:
#'   \item{fileTable}{A table with relevant data of each FASTA file generated in execution,
#'     including their associated strand, mean read length, total reads and total haplotypes obtained.}
#'   \item{poolTable}{A table with the number of total trimmed reads and the yield of the process by pool.}
#'
#'   After execution, a FASTA file for each combination of strand, MID and pool will be saved in a newly
#'   created trim folder, including its associated reads.
#'   Additionaly, some report files will be generated in a reports folder:
#'   \enumerate{
#'   \item{\code{AmpliconLengthsRprt.txt}: Includes the amplicon lengths of both strands
#'     for each sample (with their corresponding MID identifier).}
#'   \item{\code{AmpliconLengthsPlot.pdf}: Includes a barplot for each sample representing the amplicon
#'     lengths of both strands.}
#'   \item{\code{SplitByPrimersOnFlash.txt}: Includes a table of reads identified by primer, total reads identified by patient
#'     and the yield by pool.}
#'   \item{\code{SplitByPrimersOnFlash.pdf,SplitByPrimersOnFlash-hz.pdf}: Includes some plots representing primer matches
#'     by patient (in nº of reads) and the coverage of forward/reverse matches by pool.}
#'   \item{\code{SplittedReadsFileTable.txt}: A file containing the same information as \code{fileTable}.}}
#'
#' @import ShortRead
#' @import Biostrings
#' @import QSutils
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @seealso \code{\link{demultiplexMID}},\code{\link{primermatch}}
#' @examples
#' # Set parameters
#' prmm <- 3
#' min.len <- 180
#' # The expected window for primer sequences will depend on the presence of
#' # adapters, MID sequences and/or M13 primer.
#' target.st <- 1
#' target.end <- 100
#' splitDir <- "./splits"
#' # Save the file names with complete path
#' splitfiles <- list.files(splitDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' # Get data
#' samples <- read.table("./data/samples.csv", sep="\t", header=T,
#'                       colClasses="character",stringsAsFactors=F)
#' mids <- read.table("./data/mids.csv", sep="\t", header=T,
#'                    stringsAsFactors=F)
#' pm.res <- demultiplexPrimer(splitfiles,samples,primers,prmm,min.len,target.st,target.end)
#' @author Alicia Aranda


demultiplexPrimer <- function(splitfiles,samples,primers,prmm=3,min.len=180,target.st=1,target.end=100){

  # Si la ruta on es troben els fitxers de la carpeta splits no està ben especificada, intenta buscar la
  # ruta correcta a partir del directori de treball
  # Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
  if(length(splitfiles)==0) {
    splitfiles <- list.files(paste(getwd(),"/splits",sep=''))
    if(length(splitfiles)==0) {
      stop("Please indicate correct path for demultiplexed reads by MID.\n")
    }
  }

  # Si la ruta on es troben els fitxers data no està ben especificada, atura l'execució
  # i mostra un missatge d'error
  if(length(samples)==0|length(primers)==0) {
    stop("Please check data folder files, something is missing.\n")
  }

  # Si cap dels fitxers indicats a la carpeta splits no existeix, atura l'execució i
  # mostra un missatge d'error
  if(any(!file.exists(splitfiles))){
    stop(paste(splitfiles[!file.exists(splitfiles)],"does not exist.\n"))
  }

  # Si no existeixen les carpetes on es guarden els fitxers resultants de la funció,
  # es generen automàticament a la carpeta de treball
  if(!dir.exists("./trim")) {
    dir.create("./trim")}
  trimDir <- "./trim"

  if(!dir.exists("./reports")) {
    dir.create("./reports")}
  repDir <- "./reports"

### Els primers específics emprats en l'amplificació (que afegeixen la cua M13 per seqüenciar) s'han d'eliminar
# dels amplicons a estudiar, ja que a l'emprar el mateix primer per totes les mostres es perden tots els possibles
# mismatches que determinen variabilitat genètica.
# Per tant, es calculen les posicions de tall pels primers específics FW i RV. A la posició 5' del primer FW se li
# suma la seva longitud per determinar el punt de tall 5' de l'amplicó.
# A la posició 5' del primer RV se li resta la longitud per determinar el punt de tall 3' de l'amplicó.
primers$FW.tpos <- primers$FW.pos+nchar(primers$Primer.FW)
primers$RV.tpos <- primers$RV.pos-nchar(primers$Primer.RV)

### Inicialitzacions
# Calcula el nº de mostres (files del fitxer samples) i fa el doble
Ns <- nrow(samples)*2
# Genera una matriu amb el doble de files que el total de mostres i 4 columnes
# Asigna el nom a les columnes de la matriu on s'aniran afegint els resultats
pr.res <- matrix(0,nrow=Ns,ncol=4,dimnames=list(NULL,c("Tot.reads","matches","shorts","fn.reads")))

# Genera un data frame amb el doble de files que el total de mostres i 9 columnes amb els noms dels arguments
# 'character()' aplicat sobre un nombre enter retorna tants caràcters buits com el nº indiqui
# 'integer()' sobre un nombre enter retorna tants 0s com el nº indiqui
FlTbl <- data.frame(File.Name=character(Ns),Pat.ID=character(Ns),
                    Ampl.Nm=character(Ns),Pr.ID=integer(Ns),
                    Str=character(Ns),Pos=integer(Ns),Len=integer(Ns),
                    Reads=integer(Ns),Hpls=integer(Ns),
                    stringsAsFactors=FALSE)

# Guarda el nom dels pools
pools <- unique(samples$Pool.Nm)

# Guarda un data frame amb 3 columnes indicant el total de reads assignats a algun MID, el total de reads assignats
# als primers (cadenes forward o reverse) i el rendiment calculat, tot això per cada regió avaluada
PoolTbl <- data.frame(MIDReads=numeric(length(pools)),PrimerReads=0,Pct=0,row.names=pools)

# Genera els fitxers resultants en format .txt i .pdf que es guardaran a la carpeta reports
sink(file.path(repDir,"AmpliconLengthsRprt.txt"))
pdf(file.path(repDir,"AmpliconLengthsPlot.pdf"),paper="a4",
    width="6",height=10.5)
par(mfrow=c(2,1))

### Loop sobre el total de pools emprats (extrets del fitxer samples)
k <- 0
for(i in 1:length(pools))
{ # Identifica les mostres que corresponen al pool avaluat
  idx <- which(samples$Pool.Nm==pools[i])

  # Guarda el noms dels fitxers de les mostres (MID) corresponents al pool avaluat
  flnms <- paste(pools[i],".MID",samples$MID[idx],".fna",sep="")

  # Guarda el pool que s'està avaluant
  pool<-pools[i]
  # La funció 'environment()' permet assignar totes les variables locals d'aquesta funció
  # a la subfunció 'primermatch()' per tal d'evitar errors.
  environment(primermatch) <- environment()
  # Aplica la funció 'primermatch()' del mateix paquet sobre cadascuna de les mostres
  # del pool avaluat
  foreach(j=1:length(idx)) %do% {primermatch(j,idx,flnms,pool)}
}
# Tanca els fitxers .txt i .pdf generats
sink()
dev.off()

# Calcula el rendiment d'aquest pas, dividint els reads que s'han pogut assignar a alguna de les cadenes
# entre el total de reads de cada pool o regió (sumatori de tots els MIDS de cada pool)
PoolTbl$Pct <-  PoolTbl$PrimerReads/PoolTbl$MIDReads*100

### Sincronització d'estructures
# Guarda la columna de la taula de resultats amb els identificadors dels primers majors a 0
# en cas de trobar-ne algun que sigui 0
if(any(FlTbl$Pr.ID==0)) {
  fl <- FlTbl$Pr.ID>0
  # Guarda totes les entrades de la taula que compleixin la condició anterior (en aquest cas totes)
  FlTbl <- FlTbl[fl,]
  # Guarda també les entrades de la taula amb els nº de reads que compleixen la primera condició
  pr.res <- pr.res[fl,]
}
# # Guarda tots els identificadors dels pacients concatenas amb la regió del HBV avaluada
# anms <- paste(FlTbl$Pat.ID,FlTbl$Ampl.Nm,sep=".")
# # Assigna a la taula de mostres (variable samples) els noms de les files, que corresponen de nou als
# # identificadors dels pacients amb la regió avaluada (tot i que en aquest cas extrau les dades de la taula de mostres)
# rownames(samples) <- paste(samples$Patient.ID,samples$Primer.ID,sep=".")

### Gràfics de resultats
# Generació de les paletes de colors
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")

# Guarda en dues variables diferents les entrades de la taula de resultats corresponents a les cadenes forward i reverse
fw.idx <- which(FlTbl$Str=="fw")
rv.idx <- which(FlTbl$Str=="rv")
# Genera un data frame amb dades de les dues taules de resultats que inclou, només per les cadenes forward:
# ID dels pacients, regió amplificada, total de reads de la mostra, reads eliminats per longitud curta,
# reads associats a cadascuna de les dues cadenes (per separat) i sumatori del total de reads en el MID avaluat (mostra) que s'han assignat
mres <- data.frame(PatID=FlTbl$Pat.ID[fw.idx],
                     PrimerID=FlTbl$Ampl.Nm[fw.idx],
                     Treads=pr.res[fw.idx,1],
                     Shorts=pr.res[fw.idx,3]+pr.res[rv.idx,3],
                     FW.match=pr.res[fw.idx,4],
                     RV.match=pr.res[rv.idx,4],
                     Fn.reads=pr.res[fw.idx,4]+pr.res[rv.idx,4],
                     stringsAsFactors=FALSE)

# Calcula el sumatori dels reads assignats a cadena forward o reverse en funció dels pacients, és a dir,
# suma els reads assignats per a les dues regions del HBV avaluades
# 'tapply()' calcula el sumatori (sum) del nº reads segons els pacients
T.reads <- apply(mres[,5:6],2, function(x)
  tapply(x,mres$PatID,sum))

# Aquest condicional només s'aplica si a la taula només hi ha un sol pacient
if(length(unique(mres$PatID))==1)
{ # En aquest cas, es genera una matriu d'una sola fila amb els valors de reads assignats a les dues cadenes,
  # indicant el pacient i les columnes de la taula T.reads
  x <- matrix(T.reads,nrow=1)
  rownames(x) <- mres$PatID[1]
  colnames(x) <- names(T.reads)
  T.reads <- x
}

# Genera el fitxer .pdf que es guardarà a la carpeta reports
pdf(file.path(repDir,"SplitByPrimersOnFlash.pdf"),paper="a4",width=6,height=11)
par(mfrow=c(2,1),mar=c(7,4,4,2)+0.1)

# Defineix el límit de l'eix Y del gràfic a partir del sumatori de reads assignats per pacient (no segons el pool)
ymx <- max(rowSums(T.reads))*1.2
# Gràfic de barres amb la trasposada de la taula T.reads, per representar els reads assignats a la cadena forward
# i reverse per cada pacient, independentment de la regió del HBV avaluada
# En aquest gràfic es representa una única barra per pacient i d'indica amb colors diferents les assignacions up i down
barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

# Genera un altre gràfic amb les mateixes dades, però aquest cop es representen dues barres per pacient, una per
# cada cadena forward o reverse, i diferenciant les dues regions del HBV avaluades
res.mat <- mres[,5:6]
ymx <- max(res.mat)*1.2
# També es defineixen els noms de l'eix X com el ID dels pacients i la regió HBV avaluada
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
dev.off()


# Genera el fitxer .txt que es guardarà a la carpeta reports
sink(file.path(repDir,"SplitByPrimersOnFlash.txt"))
# Afegeix la taula on es registren, per cada MID avaluat, els reads assignats a forward i reverse
cat("\nTable of reads identified by primer\n\n")
print(mres)
cat("\n")
# També afegeix la taula on s'indica la longitud dels reads i els haplotips detectats per cadena (treient el nom del fitxer)
print(FlTbl[,-1])
# Inclou la taula on es mostren els reads assignats per pacient, sumant les dues cadenes
cat("\nTotal reads identified by patient\n\n")
print(T.reads)
# Mostra el rendiment d'aquest pas: els reads totals per pool (tots els MIDS de cada pool) i els que s'han assignat
cat("\nYield by pool\n\n")
print(PoolTbl)
sink()

### Guarda les taules de resultats que retornarà la funció
result <- list(fileTable=FlTbl,poolTable=PoolTbl)

# Guarda un altre fitxer .txt amb les dades incloses a cada fitxer .fna generat a la carpeta trim
sink(file.path(repDir,"SplittedReadsFileTable.txt"))
print(FlTbl)
sink()

### Genera els gràfics de barres anteriors però en un full A4 horitzontal
pdf(file.path(repDir,"SplitByPrimersOnFlash-hz.pdf"),paper="a4r",width=10,height=6)
par(mar=c(7,4,4,2)+0.1)

ymx <- max(res.mat)*1.2
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))

ymx <- max(rowSums(T.reads))*1.2
bp <- barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx),xaxt="n")
axis(side=1,at=bp,rownames(T.reads),cex.axis=0.6,las=2)
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

### També genera dos gràfics boxplot per representar els reads assignats a la cadena forward i reverse per pool
par(mfrow=c(1,2))
# Defineix el límit màxim de l'eix Y
ymx <- max(c(res.mat[,1],res.mat[,2]))
# Defineix els noms de les files de la taula de primers, amb la regió avaluada
rownames(primers) <- primers$Ampl.Nm
# Guarda el nom de la regió avaluada de totes les entrades classificades a la cadena forward
reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="fw"],"Region"]
# Genera el boxplot dels reads assignats a la cadena forward en funció de la regió o pool
boxplot(res.mat[,1]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,1],
       pch="+",cex=0.8)
title(main="Primers identified on forward reads")

# Guarda el nom de la regió avaluada de totes les entrades classificades a la cadena reverse
reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="rv"],"Region"]
# Genera el boxplot dels reads assignats a la cadena reverse en funció de la regió o pool
boxplot(res.mat[,2]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,2],
       pch="+",cex=0.8)
title(main="Primers identified on reverse reads")
dev.off()

# Borra les taules guardades a l'entorn global per evitar que augmenti la memòria d'execució, ja que estan
# guardades a la variable results.
rm(FlTbl,PoolTbl,pr.res,envir = globalenv())
return(result)
}
