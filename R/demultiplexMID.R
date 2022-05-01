#' demultiplexMID
#'
#' Demultiplex reads by identifying MID sequences within windows of expected positions in the sequenced reads.
#'   MIDs are 10 base-length oligonucleotides that allow the identification of samples from different patients
#'   or origins.
#'
#' It is important to note that MID sequences will be not trimmed from reads, they are only identified for
#'   association with each sample.
#'
#' @details Files indicated in \code{flashffiles} must be located in a directory named flashFiltDir, and
#'   a directory for resulting FASTA files with demultiplexed reads must be named as splitDir.
#'   Also, a reports folder must be created in the project environment, whose path will be
#'   named as repDir.
#'
#' @param flashffiles Vector including the names of FLASH filtered files, with fastq extension.
#' @param samples Data frame with relevant information about the samples of the sequencing experiment, including
#'   \code{Patient.ID, MID, Primer.ID, Region, RefSeq.ID}, and \code{Pool.Nm} columns.
#' @param primers Data frame with information about the \emph{primers} used in the experiment, including
#'   \code{Ampl.Nm, Region, Primer.FW, Primer.RV, FW.pos, RV.pos, FW.tpos, RV.tpos, Aa.ipos},
#'     and \code{Aa.lpos} columns.
#' @param mids Data frame with the association between MID identifiers and their sequences.
#' @param maxdif Number of mismatches allowed between MID and read sequences.
#' @param mid.start Expected start position for MID in sequence.
#' @param mid.end Expected end position for MID in sequence.
#'
#' @return A list containing the following:
#'   \item{nreads}{A table with the number of reads
#'     identified for each MID.}
#'   \item{by.pools}{A table with the coverage of reads demultiplexed by pool.}
#'   After execution, a FASTA file for each combination of MID and pool will be saved in the splits folder,
#'   including its associated reads. Additionaly, two report files will be generated:
#'   \enumerate{
#'   \item \code{SplidByMIDs.barplots.pdf}: Includes a first barplot representing \code{nreads} data values,
#'     and a second plot with the \code{by.pools} data values.
#'   \item \code{SplidByMIDs.Rprt.txt}: Includes the same data tables returned by the function.}
#'
#'
#' @import ShortRead
#' @importFrom RColorBrewer brewer.pal
#' @seealso \code{\link{FiltbyQ30}}
#' @export
#' @examples
#' maxdif <- 1
#' mid.start <- 1
#' mid.end <- 40
#' flashFiltDir <- "./flashFilt"
#' flashffiles <- list.files(flashFiltDir)
#' samples <- read.table("./data/samples.csv", sep="\t", header=T,
#'                       colClasses="character",stringsAsFactors=F)
#' mids <- read.table("./data/mids.csv", sep="\t", header=T,
#'                    stringsAsFactors=F)
#'
#' demultiplexMID(flashffiles,samples,mids,maxdif,mid.start,mid.end)
#'

demultiplexMID <- function(flashffiles,samples,mids,maxdif=1,mid.start=1,mid.end=40){

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

  # Guarda la columna Pools (regions 5' o preS1). 'unique' per eliminar elements duplicats (només hi haurà 2)
  pools <- unique(samples$Pool.Nm)

###  Inicialitzacions
# 'paste' per concatenar strings separats per punt
# Així s'identifica cada MID (mostra) amb la regió avaluada (pool)
nms <- unique(paste(samples$Pool.Nm,samples$MID,sep="."))

# Genera un data frame amb 3 columnes que indiquen el nom del pool, el MID associat i el nº de reads,
# que inicialment serà 0
nreads <- data.frame(Pool.Nm=samples$Pool.Nm,MID=samples$MID,Reads=0,stringsAsFactors=FALSE)
# Indica com a nom de les files de la taula l'indicador nms (del tipus Pool.MID)
rownames(nreads) <- nms

# Guarda un data frame que inclorà la quantificació de reads totals, els que no s'han assignat (MID0)
# i els que sí s'han assignat
by.pools <- data.frame(TotReads=integer(length(pools)),NoMID=integer(length(pools)),row.names=pools)

###  Identificar MIDs en pool
# Amb 'sapply()' guarda quins elements de la columna Pool coincideixen amb cadascun
# dels pools avaluats (quines mostres son de cada pool)
idxs <- sapply(c(1:length(pools)),function(x) which(samples$Pool.Nm==pools[x]))
# Associa el nom dels pools amb la taula anterior en funció de la classe de l'objecte obtingut
if(class(idxs)[1]=='matrix') {colnames(idxs)<-pools} else names(idxs)<-pools

# Guarda els MIDS que corresponen a cadascun dels pools
mid.sets <- tapply(idxs,samples$Pool.Nm,function(x) unique(samples$MID[x]))

flashfdata <- sapply(c(1:length(pools)),function(i)
  FastqStreamer(file.path(flashFiltDir,flashffiles[i])))

###  Loop sobre pools en samples
for(i in 1:length(pools)) {

  ###  Identificar fastq del pool
  # Busca dins dels arxius de la carpeta flashDir el fitxer per cada pool (regió VHB).
  # 'grep' serveix per buscar un patró dins d'un vector de caracters, i retorna l'index del vector que encaixa amb la cerca
  ip <- grep(pools[i],flashffiles)

  # Si no troba l'arxiu de la iteració que està fent o aquest no té MIDS associats passa a la següent iteració.
  if(length(ip)==0|length(idxs[,ip])==0) next

  reads <- 0 # nº inicial de reads a 0
  app.flag <- m0.app.flag <- FALSE

  strm <- flashfdata[[i]]
  ###  Carrega el fitxer fastq per chuncks
  # Itera sobre tots els blocs del fastq que inclouen: id, seq, +, qualitats
  while(length(sqq <- yield(strm)))
  { seqs <- sread(sqq) # Guarda la seq de cada conjunt o chunck

  ###  Actualitza el nº total de reads
  reads <- reads + length(seqs)
  ###  Es descarten aquells reads de longitud menor a 150
  seqs <- seqs[ width(seqs) > 150 ]

  ###  Loop sobre els MIDs del pool avaluat
  for(j in mid.sets[[i]]){
    # Guarda en una variable la seq inclosa en les posicions indicades.
    # En aquest cas, es busca el MID entre les posicions 1-40 (start-end), definides
    # a l'arxiu de paràmetres.
    sbsq <- subseq(seqs,start=mid.start,end=mid.end)
    k <- which(mids$MID.ID==j) # Guarda els MIDs que coincideixen amb l'avaluat.
    pr.up <- mids$MID.Seq[k] # Guarda la seq del MID que ha coincidit.

    # Recompte de coincidències entre la seq 1-40 extreta i la del MID que ha
    # coincidit abans, amb màxim 1 diferència de mismacth.
    up.matches <- vcountPattern(pattern=pr.up,subject=sbsq,
                                max.mismatch=maxdif,fixed=TRUE)

    # Només guarda els que tinguin més d'una coincidència.
    flags <- up.matches>=1
    if(sum(flags)){ # La funció 'sum' assegura que tinguem un valor major a 1

      # Guarda en un arxiu .fna la seq dels reads que s'han associat al MID
      # que s'està avaluant, amb nom `pool.nºMID.fna`, en el directori splits
      KK <- which(nms==paste(pools[i],j,sep="."))[1]
      nreads$Reads[KK] <-  nreads$Reads[KK]+sum(flags)
      up.flnm <- paste(pools[i],".MID",j,".fna",
                       sep="")
      writeXStringSet(seqs[flags],file.path(splitDir,up.flnm),
                      append=app.flag)
      # Actualitza les seqs dels reads que no han fet match en aquesta iteració
      seqs <- seqs[!flags]}
  }

  by.pools$TotReads[i] <- reads # Sumatori dels reads obtinguts en el pool avaluat

  if(length(seqs)){ # En cas que s'hagin iterat tots els MIDs i quedin reads sense assignar:
    # Guarda un altre arxiu .fna amb els reads que no s'han assignat a MID i ho anomena MID0
    mid0.flnm <- paste(pools[i],"MID0.fna",sep=".")
    writeXStringSet(seqs,file.path(splitDir,mid0.flnm),
                    append=m0.app.flag)
    # Canvia aquest paràmetre perquè les properes iteracions les noves seqs s'afegeixin al final del fitxer
    m0.app.flag <- TRUE

    # Recompte els reads no assignats a cap MID per aquell pool
    by.pools$NoMID[i] <- by.pools$NoMID[i]+length(seqs)
  }
  # Canvia aquest paràmetre perquè les properes iteracions les noves seqs s'afegeixin al final del fitxer
  app.flag <- TRUE
  }
  # Tanca el fastq avaluat
  close(strm)
}

# Guarda una nova columna a la taula de resultats on s'inclou el nº total de reads assignats a un MID
# per pool
by.pools$MIDReads <- tapply(nreads$Reads,factor(nreads$Pool.Nm,levels=pools),sum)
# Borra els noms de les files de la taula nreads i guarda el resultat retornat per la funció
rownames(nreads) <- NULL
result <- list(nreads=nreads,by.pools=by.pools)

# Genera un fitxer .txt que es guardarà a la carpeta reports, on es guadarà en format taula els resultats
# representats als gràfics
sink(file.path(repDir,"SplidByMIDs.Rprt.txt"))
cat("\nCoverage by Pool and MIDs\n\n")
print(nreads)
cat("\n")
print(by.pools)
sink()

# Genera el pdf indicar que es guardarà a la carpeta reports
pdf.flnm <- file.path(repDir,"SplidByMIDs.barplots.pdf")
pdf(pdf.flnm,paper="a4",width=5.5,height=10)
par(mfrow=c(2,1),mar=c(7.5,4,4,2)+0.1)

pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
# Genera un gràfic de barres amb les dades del nº de reads per MID
bp <- barplot(nreads$Reads,col="lavender",border="navy")
axis(1,at=bp,las=2,rownames(nreads),cex.axis=0.8)
title(main="Coverage by Pool and MID",line=1.5)

# Defineix el límit de l'eix y a partir del màxim de reads totals obtinguts
ymx <- max(data.matrix(by.pools))*1.2
# Genera un altre gràfic de barres amb les dades del nº de reads per pool (assignats a MID i no assignats)
bp <- barplot(t(data.matrix(by.pools)),beside=TRUE,ylim=c(0,ymx),
              col=pal2[1:3],border=pal1[1:3])
legend("top",horiz=TRUE,fill=pal2[1:3],legend=colnames(by.pools),cex=0.8)
title(main="Coverage by Pool",line=1.5)

dev.off()

return(result)
}
