#' QCscores
#'
#' Applies \code{\link{FastqStreamer}} over a fastq file and returns quality
#'   control mesures (Phred scores), either by position or by read.
#'
#' This function is only defined for correct execution of \code{\link{PoolQCbyPos}} and \code{\link{PoolQCbyRead}}
#'  functions from the same package.
#'
#'
#' @param flnm Fastq file to be evaluated.
#' @param ln Amplicon length (by default 301 nt).
#' @param byPos Logical indicating if QC by position should be done.
#' @param byRead Logical indicating if QC by read should be done.
#' @note Arguments \code{byPos} and \code{byRead} are mutually exclusive.
#' @return If argument \code{byPos=TRUE}, the function returns a list including
#'   the following parameters:
#'   \enumerate {
#'   \item \code{fvnq}: A matrix with Phred quality scores across each nucleotide base in the reads.
#'     Columns indicate base position and rows indicate 0.05, 0.25, 0.5, 0.75 and 0.95 Phred quantiles.
#'   \item \code{fvnl}: A vector with normalized read lengths for each Phred quantile.
#'   \item \code{all.ln}: A vector with all read lengths.}
#' @return If argument \code{byRead=TRUE}, the function returns a list including
#'   the following parameters:
#'   \enumerate {
#'   \item \code{all.ln}: A vector with all read lengths.
#'   \item \code{all.ln30}: A vector with the number of bases below Q30 for each read.}
#' @seealso \code{\link{PoolQCbyPos}}, \code{\link{PoolQCbyRead}}, \code{\link{QCplot}}
#' @export

QCscores <- function(flnm,ln=301,byPos=FALSE,byRead=FALSE) # ln és 301 per defecte, però en realitat és la longitud de l'amplicó
{ # Defineix una matriu de dimensions 5xln de 0s
  fnm.q <- matrix(0,nrow=5,ncol=ln)
  # Vector de 5 0s
  fnm.l <- numeric(5)
  # Variable numèrica
  nrds <- numeric()
  # Variable de nombre enter
  all.ln <- integer()

  ### Aplica streamer (iteració) en el fitxer fastq. Cada iteració serà de la mida de la variable chunck.sz
  strm <- FastqStreamer(flnm,n=chunck.sz) # chunck.sz definit en el fitxer principal

  # Variable per a la iteració sobre els fastq
  nchk <- 0

# Si volem obtenir el QC per posició s'executa aquest codi:
if(byPos){

  ### Carrega el fitxer fastq per chuncks
  while(length(sqq <- yield(strm)))
  { nchk <- nchk+1 # Nº de cada iteració (chunck del fastq avaluat)
  nrds[nchk] <- length(sqq) # Guarda el nº de reads del chunck avaluat

  ###  Phred scores. Codi ASCII-33
  # Funció 'quality()' retorna el valor de qualitat per posició dels reads
  # Funció 'as()' permet fer coerció del resultat a matriu
  phrsc <- as(quality(sqq),"matrix")

  # Guarda el valor mínim entre la variable ln i les columnes de la matriu amb les qualitats dels strings
  nc <- min(ln,ncol(phrsc))
  # De les columnes 1 al mínim assignat abans, aplica els quantils del vector a les columnes
  # de la matriu amb les qualitats phrsc, i ho multiplica pel total de reads (normalització per chuncks)
  # Cada fila es un quantil i cada columna es un cicle de seqüenciació (posició del read)
  fnm.q[,1:nc] <- fnm.q[,1:nc] +
    apply(phrsc,2,quantile,p=c(0.05,0.25,0.5,0.75,0.95),
          na.rm=TRUE)[,1:nc] * nrds[nchk]

  ###  Longituds de seqüència
  sqln <- width(sqq) # Cicles de seqüenciació (longitud)
  all.ln <- c(all.ln,sqln) # Vector amb el nº de cicles (longitud) com a nombres enters. S'afegiran els nous valors amb cada iteració
  # Aplica els quantils per guardar-los en la llista de 5 columnes de 0 (1 columna per quantil)
  # i guarda els resultats multiplicats per la longitud de la seq (normalització)
  fnm.l <- fnm.l + quantile(sqln,p=c(0.05,0.25,0.5,0.75,0.95),
                            na.rm=TRUE) * nrds[nchk]
  }
  # Retorna una llista amb 3 matrius:
  # 1- Quantils de phred score per posició entre total de reads
  # 2- Quantils aplicats a la longitud dels reads entre el total
  # 3- Nº de cicles de seqüenciació (longitud dels reads)
  result <- list(fvnq=fnm.q/sum(nrds),fvnl=fnm.l/sum(nrds),all.ln=all.ln)}

# Si volem obtenir el QC per read i no per posició, s'executa aquest codi:
else if(byRead){
  # Variable de nombre enter
  all.nl30 <- integer()
  # Itera sobre els chuncks del fastq
  while(length(sqq <- yield(strm)))  {
  ### Phred scores. Codi ASCII-33
  # Funció 'quality()' retorna el valor de qualitat dels strings
  # Funció 'as()' permet fer coerció del resultat a matriu
  phrsc <- as(quality(sqq),"matrix")

  ### A la matriu amb els Phred scores aplica el sumatori de les bases amb puntuació menor a 30 (per sota de Q30)
  nl30 <- apply(phrsc,1,function(x) sum(x<30,na.rm=TRUE)) # na.rm per no tenir en compte missing values
  # La variable de nombre enter ara inclourà les bases per sota de Q30
  all.nl30 <- c(all.nl30,nl30)
  ###  Longituds de seqüència
  all.ln <- c(all.ln, width(sqq)) # Vector amb el nº de cicles com a nombres enters
  }
  # Retorna una llista formada per les 2 matrius: amb el nº de cicles (longitud dels reads)
  # i les bases per sota de Q30
  result <- list(all.ln=all.ln,all.nl30=all.nl30)
}
  close(strm)
  return(result)
}
