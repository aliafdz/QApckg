#' @title Check for insertions and deletions in samples
#'
#' @description Generates insertion/deletion plots for each evaluated sample
#'   from consensus haplotypes' sequences.
#'
#' @param machfiles Vector including the paths of files generated by \code{\link{ConsHaplotypes}} function,
#'   with fna extension.
#' @param pm.res The list returned by \code{\link{demultiplexPrimer}}, including \code{fileTable}
#'   and \code{poolTable} data frames.
#'
#' @return After execution, a file named \code{GapsBarPlots.pdf} will be saved in the reports folder,
#'   including the plots generated for all samples. In each plot, insertions are represented by red lines
#'   and deletions are represented by blue lines.
#'
#' @note This function is designed to be applied at the end of the quality assessment analysis and requires
#'   the previous execution of \code{\link{demultiplexPrimer}} and \code{\link{ConsHaplotypes}} functions
#'   from the same package.
#'
#' @import stringr
#' @import Biostrings
#' @export
#' @examples
#' ## Execute demultiplexing by primer
#' splitDir <- "./splits"
#' splitfiles <- list.files(splitDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' pm.res <- demultiplexPrimer(splitfiles,samples,primers)
#'
#' ## Obtain consensus haplotypes (default parameters)
#' trimDir <- "./trim"
#' trimfiles <- list.files(trimDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' int.res <- ConsHaplotypes(trimfiles, pm.res, thr, min.seq.len)
#'
#' ## Apply function
#' mach.Dir <- "./MACH"
#' machfiles <- list.files(mach.Dir,recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
#' PlotInDels(machfiles,pm.res)

PlotInDels <- function(machfiles,pm.res){

  # Si la taula pm.res no s'ha inclòs o no existeix, retorna un missatge d'error
  if(missing(pm.res)|!exists('pm.res')|length(pm.res)==0) {
    stop("The list obtained by demultiplexPrimer function is needed.\n")
  }
  if(!is.list(pm.res)){
    stop("pm.res argument must be a list.\n")
  }

  # Si la ruta on es troben els fitxers de la carpeta splits no està ben especificada, intenta buscar la
  # ruta correcta a partir del directori de treball
  # Si tot i així no troba fitxers, atura l'execució i mostra un missatge d'error
  if(length(machfiles)==0) {
    machfiles <- list.files(paste(getwd(),"/MACH",sep=''))
    if(length(machfiles)==0) {
      stop("Please indicate correct path for intersected reads.\n")
    }
  }

  # Filtra dels fitxers de la carpeta MACH només aquells resultants de la intersecció
  machfiles <- machfiles[grep("MACHpl02",machfiles)]

  # Si cap dels fitxers indicats a la carpeta splits no existeix, atura l'execució i
  # mostra un missatge d'error
  if(any(!file.exists(machfiles))){
    stop(paste(machfiles[!file.exists(machfiles)],"does not exist.\n"))
  }

  # Si la carpeta reports no s'ha creat a l'entorn de treball, la genera per emmagatzemar els
  # gràfics de resultats.
  if(!dir.exists("./reports")) {
    dir.create("./reports")}
  repDir <- "./reports"

# En aquest cas només fa falta la taula FlTbl del demultiplexat per primers
FlTbl <- pm.res$fileTable

# Retorna els indexs de la taula FlTbl derivats d'ordenar les mostres en funció de l'ID del
# pacient, de l'amplicó (regió) avaluat i la cadena forward o reverse
o <- order(FlTbl$Pat.ID,FlTbl$Ampl.Nm,FlTbl$Str)
# Reordena les entrades de la taula en funció dels indexs anteriors, de manera que s'agrupen
# les entrades corresponents al mateix pacient i en segon ordre a la mateixa regió avaluada
FlTbl <- FlTbl[o,]


### Generació dels noms dels fitxers fasta
# Guarda els indexs de la taula FlTbl segons la cadena asignada als reads
idx.fw <- which(FlTbl$Str=="fw")
idx.rv <- which(FlTbl$Str=="rv")

# Guarda els noms dels fitxers fasta inclosos a la carpeta MACH generats després
# de fer la intersecció entre cadenes (MACHpl02).
flnms <- paste("MACHpl02",FlTbl$Pat.ID[idx.fw],
               FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
flnms <- machfiles

# Guarda les concatenacions dels ID dels pacients amb la regió amplificada units per guionet
pnms <- paste(FlTbl$Pat.ID[idx.fw],FlTbl$Ampl.Nm[idx.fw],sep=" - ")

## Genera el fitxer .pdf on es guardaran els gràfics de resultats d'aquest codi
# Es generara un gràfic per cada mostra avaluada
pdf(file=file.path(repDir,"GapsBarPlots.pdf"),paper="a4",width=7,
    height=10)
par(mfrow=c(4,1))
par(mar=c(4.5, 4, 3, 4) + 0.1)

## Bucle sobre totes les mostres (2 per pacient), que coincideix amb el nº de fitxers MACHpl02
foreach(i=1:length(flnms)) %do% {
  ## Aplica la funció del paquet QSutils per llegir el fitxer fasta de la mostra avaluada
  lst <- ReadAmplSeqs(flnms[i])
  # Guarda totes les seqüències dels haplotips del fitxer avaluat
  seqs <- as.character(lst$hseqs)
  # Guarda la llista amb el nº de reads de tots els haplotips
  nr <- lst$nr

  # Guarda l'index de l'haplotip amb màxim nº de reads
  imast <- which.max(nr)
  # Guarda la seqüència de l'haplotip amb màxim nº de reads (màster)
  rsq <- seqs[imast]
  # Guarda la longitud de la seq de l'haplotip màster
  xmx <- nchar(rsq)

  ## Separa la seqüència de l'haplotip màster per nucleòtids
  rnt <- strsplit(rsq,split="")[[1]]
  # Guarda una variable buida amb tants 0 com nucleòtids hi ha a la seq
  tp <- integer(length(rnt))
  # Indica quins elements de la seq corresponen a insercions (símbol "-")
  ins.pos <- rnt=="-"
  # A la variable buida d'abans, indica amb un 1 les posicions amb insercions
  tp[ins.pos] <- 1

  ## Ara separa per nucleòtids les seqs de tots els haplotips de la mostra avaluada
  ntmat <- t(sapply(seqs,function(s) strsplit(s,split="")[[1]]))
  # Indica, sobre cada element de les seqs, quins corresponen a delecions ("-")
  del.pos <- apply(ntmat,2,function(vnt) any(vnt=="-"))
  # A la variable d'abans, indica amb un 2 les posicions on no hi ha insercions
  # però sí delecions
  tp[!ins.pos & del.pos] <- 2
  # Suma una unitat a tots els elemnts de la variable tp, de manera que un 1 indica
  # on hi ha nucleòtids, 2 on hi ha insercions i 3 on hi ha delecions
  tp <- tp+1

  # Guarda una altra variable buida amb tants 0 com nucleòtids hi ha a la seq màster
  nb <- integer(length(rnt))
  # En cas de que es trobin insercions en alguna de les seqs,
  if( sum(ins.pos) )
  { # A les posicions on s'hagin detectat insercions, en la variable buida nb
    # es guarda el nº de reads dels haplotips on s'ha donat la inserció (i per tant no hi ha gap)
    nb[ins.pos] <- apply(ntmat[,ins.pos,drop=FALSE],2,function(vnt)
      sum(nr[vnt!="-"]))
  }
  # En cas de que es trobin delecions en alguna de les seqs,
  if( sum(del.pos) )
  { # A les posicions on s'hagin detectat delecions i no insercions, en la variable nb
    # es guarda el total de reads dels haplotips on s'ha donat la deleció (i per tant hi ha gap)
    nb[!ins.pos & del.pos] <- apply(ntmat[,!ins.pos & del.pos,drop=FALSE],2,
                                    function(vnt) sum(nr[vnt=="-"]))
  }
  # Genera un gràfic amb la funció 'plot()' per representar els InDels de les seqs dels
  # haplotips (color vermell insercions i blau delecions), en forma de línies verticals,
  # indicant la posició de la seq en l'eix X i el nº de reads en l'eix Y
  plot(nb,type="h",col=c("black","red","blue")[tp],lwd=2,yaxt="n",
       xlab="MA position",ylab="reads",xlim=c(0,xmx),ylim=c(0,max(nb[1:xmx])))
  # Genera els intervals de l'eix Y de manera que siguin de 20 en 20
  ytk <- axTicks(2)
  # Genera l'eix Y de l'esquerra (nº reads)
  axis(side=2,at=ytk,las=2)
  # Per generar l'eix Y de la dreta (%), es divideixen els valors de l'eix Y
  # de l'esquerra i es divideix entre el total de reads de la mostra avaluada
  pct <- round(ytk/sum(nr)*100,2)
  axis(side=4,at=ytk,labels=pct,las=2,cex.axis=0.8)
  mtext("Percentage",side=4,line=3,cex=0.6)
  # El títol de cada gràfic correspon a l'ID del pacient i la regió amplificada
  title(main=pnms[i],line=1)
}

dev.off()

return(cat("The generated files are in the reports folder.\n"))
}
