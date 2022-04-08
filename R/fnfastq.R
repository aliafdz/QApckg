#' fn.fastq
#'
#' none
#'
#'
# Funció que es cridarà després sobre les dades abans i després de flash
fn.fastq <- function(flnm,ln=301) # ln és 301 per defecte, però en realitat és la longitud de l'amplicó
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
  ### Carrega el fitxer fastq per chuncks
  nchk <- 0
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
  all.ln <- c(all.ln,sqln) # Vector amb el nº de cicles com a nombres enters. S'afegiran els nous valors amb cada iteració
  # Aplica els quantils per guardar-los en la llista de 5 columnes de 0 (1 columna per quantil)
  # i guarda els resultats multiplicats per la longitud de la seq (normalització)
  fnm.l <- fnm.l + quantile(sqln,p=c(0.05,0.25,0.5,0.75,0.95),
                            na.rm=TRUE) * nrds[nchk]
  }
  close(strm)
  result <- list(fvnq=fnm.q/sum(nrds),fvnl=fnm.l/sum(nrds),all.ln=all.ln)
  # Retorna una llista amb 3 matrius:
  # 1- Quantils de phred score per posició entre total de reads
  # 2- Quantils aplicats a la longitud dels reads entre el total
  # 3- Nº de cicles de seqüenciació (longitud dels reads)
  return(result)#list(fvnq=fnm.q/sum(nrds),fvnl=fnm.l/sum(nrds),all.ln=all.ln))
}
