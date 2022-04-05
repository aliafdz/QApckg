
#p1<- executeFLASH(R1.flnms[1],R2.flnms[1],out.flnms[1])
#p2<- executeFLASH(R1.flnms[2],R2.flnms[2],out.flnms{2})

#res[1,] <- p1
#res[2,] <- p2


#' @param R1 File path with R1 reads
#' @param R2 File path with R2 reads
#' @param outfile File path for FLASH output

executeFLASH <- function(R1,R2,outfile) {

  # Concatena la ruta de l'executable flash, els paràmetres definits en el fitxer QA i la ruta dels fitxers
  # R1 i R2 del pool avaluat
  command <- paste(flash,flash.opts,R1,R2,collapse=" ")

  # La funció system() invoca una comanda especificada en l'argument
  # Executa el programa FLASH -> es guarden els fitxers següents a la carpeta global: out.extendedFrags.fastq,
  # out.hist, out.histogram, out.notCombined_1.fastq i out.notCombined_2.fastq
  es <- system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
               ignore.stdout=FALSE,invisible=TRUE)
  if(es!=0) next
  # Copia el fitxer de 'from' en 'to'. Es guarda el fitxer resultant (to) en la carpeta flash
  file.copy(from="out.extendedFrags.fastq",to=outfile,overwrite=TRUE)

  ## Llegeix el fitxer out.hist, una taula que recull els valors de 2 variables V1 i V2
  # V1 correspon a la longitud del read i V2 al nº de reads
  hstln <- read.table("out.hist",header=FALSE)

  # Recull en la columna el sumatori de la V2 de la variable anterior (extended reads)
  ext <- sum(hstln[,2])

  #  Aplica la funció FastqStreamer per iterar sobre el fitxer fastq dels reads no units per FLASH
  strm <- FastqStreamer("out.notCombined_1.fastq",n=chunck.sz)
  #  Carrega el fitxer fastq per chuncks i actualitza en nº total de reads
  noext <- 0
  while(length(sqq <- yield(strm)))
    noext <- noext + length(sqq)
  close(strm)

  # Guarda en una matriu el nº de reads units i no units per FLASH
  result <- matrix(c(ext,noext),byrow=T,ncol=2,dimnames=list('',c('Extended','NoExtd')))

return(result)

}
