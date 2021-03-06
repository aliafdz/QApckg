#' @title Run FLASH to extend R1 and R2 reads.
#'
#' @description  With R1 and R2 files generated by paired-end sequencing, the function executes the FLASH program
#'   to obtain the number of extended and not-extended reads.
#' @param R1 Path for R1 reads file.
#' @param R2 Path for R2 reads file.
#' @param flash Folder path containing the FLASH executable.
#' @param flash.opts Character indicating FLASH options that will be part of the execution command.
#' @param outfile File path for FLASH output. If it is not specified, the fastq file
#'   generated will be saved in the current working directory.
#' @return This function returns a matrix containing the number of reads extended
#'   and not extended by FLASH. Additionaly, a fastq file with extended reads
#'   will be saved to outfile path. Further FLASH output files will be saved in a new folder named "tmp".
#'
#' @note This function is defined for correct execution of \code{\link{R1R2toFLASH}} function from
#'   the same package, where all arguments are defined automatically.
#' @import ShortRead
#' @seealso \code{\link{R1R2toFLASH}}
#' @export
#' @author Alicia Aranda


executeFLASH <- function(R1,R2,flash,flash.opts,outfile='./flash.fastq') {

  if(!dir.exists("./tmp")) {
    dir.create("./tmp")}

  # Concatena la ruta de l'executable flash, els paràmetres definits en el fitxer QA i la ruta dels fitxers
  # R1 i R2 del pool avaluat. També afegeix una opció per guardar els arxius generats en una carpeta nova.
  command <- paste(flash,flash.opts,R1,R2,"--output-directory=./tmp",collapse=" ")

  # La funció system() invoca una comanda especificada en l'argument
  # Executa el programa FLASH -> es guarden els fitxers següents a la carpeta global: out.extendedFrags.fastq,
  # out.hist, out.histogram, out.notCombined_1.fastq i out.notCombined_2.fastq
  es <- system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
               ignore.stdout=FALSE,invisible=TRUE)
  if(es!=0) stop()
  # Copia el fitxer de 'from' en 'to'. Es guarda el fitxer resultant (to) en la carpeta flash
  file.copy(from="./tmp/out.extendedFrags.fastq",to=outfile,overwrite=TRUE)

  ## Llegeix el fitxer out.hist, una taula que recull els valors de 2 variables V1 i V2
  # V1 correspon a la longitud del read i V2 al nº de reads
  hstln <- read.table("./tmp/out.hist",header=FALSE)

  # Recull en la columna el sumatori de la V2 de la variable anterior (extended reads)
  ext <- sum(hstln[,2])

  #  Aplica la funció FastqStreamer per iterar sobre el fitxer fastq dels reads no units per FLASH
  strm <- FastqStreamer("./tmp/out.notCombined_1.fastq")
  #  Carrega el fitxer fastq per chuncks i actualitza en nº total de reads
  noext <- 0
  while(length(sqq <- yield(strm)))
    noext <- noext + length(sqq)
  close(strm)

  # Guarda en una matriu el nº de reads units i no units per FLASH
  result <- matrix(c(ext,noext),byrow=T,ncol=2,dimnames=list('',c('Extended','NoExtd')))

return(result)

}
