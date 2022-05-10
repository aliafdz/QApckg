#' @title Compute global yield by step
#'
#' @description Generates global yield reports for each evaluated pool from previous results.
#'
#' @param samples Data frame with relevant information about the samples of the sequencing experiment, including
#'   \code{Patient.ID, MID, Primer.ID, Region, RefSeq.ID}, and \code{Pool.Nm} columns.
#' @param filtres The data frame returned by \code{\link{FiltbyQ30}} function.
#' @param pm.res The list returned by \code{\link{demultiplexPrimer}}, including \code{fileTable}
#'   and \code{poolTable} data frames.
#' @param int.res The data frame returned by \code{\link{ConsHaplotypes}} function.
#'
#' @note This function is designed to be applied at the end of the quality assessment analysis and requires
#'   the previous execution of \code{\link{FiltbyQ30}},  \code{\link{demultiplexPrimer}} and \code{\link{ConsHaplotypes}}
#'   and functions from the same package.
#' @return After execution, two report files will be saved in the reports folder:
#'   \enumerate{
#'   \item{\code{GlobalYieldBarplots.pdf}: Includes some barplots representing the yield (in nº of
#'     reads and percentage) by each step of the quality assessment pipeline. These
#'     representation is done for all pools included in the analysis and also for global results.}
#'   \item{\code{GlobalYield-SumRprt.txt}: Summary report including global yield by analysis step in
#'     number of reads, in percentage by step and percentage referred to raw reads. }}
#' @seealso \code{\link{FiltbyQ30}},  \code{\link{demultiplexPrimer}}, \code{\link{ConsHaplotypes}}
#' @author Alicia Aranda
#' @export
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' ## Execute FLASH extension
#' runDir <- "./run"
#' runfiles <- list.files(runDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' flash <- "./FLASH/flash.exe"
#' flashres <- R1R2toFLASH(runfiles,flash)
#'
#' ## Execute Q30 filtering
#' flashDir <- "./flash"
#' flashfiles <- list.files(flashDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' filtres <- FiltbyQ30(max.pct=0.05,flashfiles,flashres)
#'
#' ## Execute demultiplexing by MID with default parameters
#' flashFiltDir <- "./flashFilt"
#' flashffiles <- list.files(flashFiltDir,recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
#' # Get data
#' samples <- read.table("./data/samples.csv", sep="\t", header=T,
#'                       colClasses="character",stringsAsFactors=F)
#' mids <- read.table("./data/mids.csv", sep="\t", header=T,
#'                    stringsAsFactors=F)
#' dem.res<-demultiplexMID(flashffiles,samples,mids)
#'
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
#' GblYield(samples, filtres, pm.res, int.res)
#'


GblYield <- function(samples,filtres,pm.res,int.res){

  # Si la taula samples no s'ha inclòs o no existeix, retorna un missatge d'error
  if(missing(samples)|length(samples)==0) {
    stop("Samples file not found.\n")
  }

  # Si la taula filtres no s'ha inclòs o no existeix, retorna un missatge d'error
  if(missing(filtres)|!exists('filtres')|length(filtres)==0) {
    stop("The data frame obtained by R1R2toFLASH function is needed.\n")
  }
  if(!is.data.frame(filtres)){
    stop("filtres argument must be a data frame.\n")
  }

  # Si la taula pm.res no s'ha inclòs o no existeix, retorna un missatge d'error
  if(missing(pm.res)|!exists('pm.res')|length(pm.res)==0) {
    stop("The list obtained by demultiplexPrimer function is needed.\n")
  }
  if(!is.list(pm.res)){
    stop("pm.res argument must be a list.\n")
  }

  # Si la taula int.res no s'ha inclòs o no existeix, retorna un missatge d'error
  if(missing(int.res)|!exists('int.res')|length(int.res)==0) {
    stop("The data frame obtained by ConsHaplotypes function is needed.\n")
  }
  if(!is.data.frame(int.res)){
    stop("int.res argument must be a data frame.\n")
  }

  # Si la carpeta reports no s'ha creat a l'entorn de treball, la genera per emmagatzemar els
  # gràfics de resultats.
  if(!dir.exists("./reports")) {
    dir.create("./reports")}
  repDir <- "./reports"


# A partir de la llista pm.res només es requereix la taula de resultats per pool, poolTable
PoolTbl<- pm.res$poolTable

# La funció 'str_extract()' permet treure, dels noms de les files de la taula, el terme
# final _S2 o _S1 (s'entén que correspon al terme "^[A-Za-z0-9\\.-]+")
rownames(filtres) <- str_extract(rownames(filtres),"^[A-Za-z0-9\\.-]+")

###  Rendiment per pas de cada pool
# Guarda els noms de les files de la taula PoolTbl (noms dels pools o regions avaluades)
p.nms <- rownames(PoolTbl)
# Guarda la taula filtres, en concret les files corresponents als pools anteriors
# (en aquest cas, tota la taula)
filtres <- filtres[p.nms,]

## Genera una taula on s'incloguin tots els resultats guardats anteriorment en columnes:
# raw reads per pool, reads estesos per FLASH, reads restants després del filtrat per Q30,
# reads restants després del demultiplexat per MIDS (total per cada pool),
# reads restants després del demultiplexat per primers
gbly <- data.frame(raw=filtres[,"Raw"],flash=filtres[,"Extended"],fQ30=filtres[,"FiltQ30"],
              MID=PoolTbl[p.nms,"MIDReads"],primer= PoolTbl[p.nms,"PrimerReads"])

# Guarda el noms dels pools corresponents a cada mostra (?
# which(samples$Patient.ID==int.res$Pat.ID[i] & samples$Primer.ID==int.res$Ampl.Nm[i])[1]
# indica quin identificador de la taula samples correspon al pacient corresponent de la taula int.res
# i també a la regió avaluada segons el valor i
# Guarda el pool corresponent després d'aplicar el condicional per a totes les files de la taula int.res (mostres)
int.res.pool <- sapply(1:nrow(int.res), function(i)
  samples$Pool.Nm[ which(samples$Patient.ID==int.res$Pat.ID[i] &
                           samples$Primer.ID==int.res$Ampl.Nm[i])[1] ])

# Guarda el total de reads després del filtre per longitud dels haplotips
gbly$filt <- tapply(int.res$all,int.res.pool,sum)[p.nms]
# Guarda el total de reads després de la intersecció entre haplotips
gbly$ints <- tapply(int.res$Fn.rd,int.res.pool,sum)[p.nms]

# Guarda el sumatori dels resultats dels dos pools avaluats
Tgbl <- colSums(gbly)
# Calcula el rendiment (en %) de tots els passos respecte l'anterior (pools sumats)
pct.gbl <- round(Tgbl[-1]/Tgbl[-length(Tgbl)]*100,2)


## Genera el fitxer pdf on es guardaran els gràfics de resultats a la carpeta reports
pdf.flnm3 <- "GlobalYieldBarplots.pdf"
pdf(file.path(repDir,pdf.flnm3),paper="a4",width=6.5,height=10)
par(mfrow=c(2,1))
# Genera les paletes de colors
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
# Total de columnes de la taula gbly, és a dir els passos totals de l'anàlisis
m <- ncol(gbly)
# Defineix el límit de l'eix Y a partir del màxim de la taula de resultats
ymx <- max(gbly)*1.2
## Gràfic de barres representant el nº de reads de cada pas de l'anàlisis per pool
bp <- barplot(t(gbly),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,ymx),xaxt="n",cex.axis=0.8,ylab="# reads")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)
abline(h=0)
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)
title(main="Yield on pools by analysis step")

## Gràfic de barres representant el rendiment (en %) de cada pas de l'anàlisis per pool
bp <- barplot(t(gbly/gbly[,1]*100),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,117),xaxt="n",cex.axis=0.8,ylab="Percentage")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)
grid(nx=NA,ny=NULL)
abline(h=0)
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)
title(main="Yield on pools by analysis step")

## Gràfic de barres representant el rendiment global de cada pas de l'anàlisis (pools sumats)
par(mar=c(5,8,4,6))
bp <- barplot(pct.gbl,col="lavender",border="navy",ylim=c(0,max(pct.gbl)),
              ylab="yield (%)")
text(bp,50,pct.gbl,col="navy",font=2,cex=0.8)
title(main="Global yield by step")

## Diagrama de caixa representant la cobertura final en funció dels reads de l'últim pas
# 'ifelse()' permet obtenir un resultat en funció del test del primer argument (TRUE o FALSE)
# En aquest cas, si la columna final reads de la taula int.res és major a 0, retorna el logaritme en base 10,
# i si no retorna com a resultat 0.
logReads <- ifelse(int.res$Fn.rd>0,log10(int.res$Fn.rd),0)
boxplot(logReads,border="gray",ylab="log10(# of reads)",outline=FALSE,
        ylim=range(logReads))
points(jitter(rep(1,nrow(int.res)),a=0.10),logReads,pch="+",cex=0.8)
title(main="Final coverage")

dev.off()


## Genera el fitxer .txt on es guardaran els resultats finals a la carpeta reports
txt.flnm2 <- "GlobalYield-SumRprt.txt"
sink(file=file.path(repDir,txt.flnm2))

cat("\n   Global yield by analysis step")
cat("\n===================================\n")
cat("\nIn number of reads:\n\n")
# Guarda la taula de nº de reads segons el pas de l'anàlisi per pool, i afegeix el sumatori
# dels dos pools a la tercera fila
# Nota: es podria substituir per la variable Tgbl
gbly <- rbind(gbly,TOTAL=Tgbl)
print(gbly)

cat("\nIn percentage by step:\n\n")
# Guarda els resultats de % dels passos en funció de l'anterior, per pool
print(round(gbly[,-1]/gbly[,-ncol(gbly)]*100,2))

cat("\nIn percentage referred to raw reads:\n\n")
# Guarda els resultats de % dels passos en funció dels raw reads, per pool
print(round(gbly/gbly[,1]*100,2))
sink()

return(cat("The generated files are in the reports folder.\n"))
}
