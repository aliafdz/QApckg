#' @title Save consensus haplotypes after sorting by mutation and abundance
#'
#' @description Sorts and renames haplotypes by the number of mutations with respect to the
#'   dominant haplotype, and by abundance, and saves their sequences in a FASTA file.
#'
#' @details This function is similar to \code{\link{SortByMutations}} function from \code{QSutils}
#'   package but has new features. For example, in this case haplotypes with a huge number
#'   of mutations (defined by \code{max.difs}) with respect to the dominant one are discarted, and columns with all gaps
#'   are eliminated. Also, the final sequences are saved in a FASTA file.
#'
#' @param flnm File path of the FASTA file that will be generated with haplotype sequences.
#' @param bseqs Character object with the haplotype alignment.
#' @param nr Vector with the haplotype counts.
#' @param max.difs Maximum number of mismatches allowed with respect to the dominant one
#'
#' @return A list containing the following:
#'   \item{bseqs}{DNAStringSet or AAStringSet with the haplotype sequences.}
#'   \item{nr}{Vector of the haplotype counts.}
#'   \item{nm}{Vector of the number of differences of each haplotype with respect to the dominant haplotype.}
#'
#'   After execution, a FASTA file named as \code{flnm} with the \code{bseqs} element will be generated.
#'
#' @import Biostrings
#' @seealso \code{\link{SortByMutations}}, \code{\link{ConsHaplotypes}}, \code{\link{IntersectStrandHpls}}
#' @export
#' @author Alicia Aranda


### Dóna nom als haplotips en funció del nombre de mutacions respecte la màster
### i de la seva frequencia poblacional, i salva les seqüències en un fitxer FASTA.
# Els arguments corresponen al nom del fitxer on es desaran els resultats, les seqs
# dels haplotips coincidents en ambdues cadenes, el sumatori de reads per haplotip coincident
# i el nº màxim de diferències permeses
SaveHaplotypes <- function(flnm="./SavedHaplotypes.fna",bseqs,nr,max.difs=250)
{
  if(!is(bseqs, "character"))
    stop("The input object must be character.\n")
  if(length(bseqs) != length(nr))
    stop("The input objects must have the same length.\n")

  ## Si només hi ha un haplotip coincident:
  if(length(bseqs)==1)
  { # Separa la seqüència per nucleòtids
    bnts <- strsplit(bseqs[1],split="")[[1]]
    # Guarda únicament els nucleòtids de la seqüència, elimina els gaps (-)
    bnts <- bnts[bnts!="-"]
    # Torna a generar la seqüència de DNA tornant a unir els nucleòtids
    bseqs <- paste(bnts,collapse="")
    # Assigna a la seq el nom corresponent: haplotip nº1 amb 0 mutacions, indicant el nº de reads,
    # i la freq serà del 100% (perquè només queda 1 haplotip)
    names(bseqs) <- paste("Hpl.0.0001",nr,100,sep="|")
    # Guarda la seq en un fitxer FASTA en la ruta indicada en l'argument de la funció
    writeXStringSet(DNAStringSet(bseqs),flnm)
    # Retorna una llista amb la seqüència i el nº de reads de l'haplotip
    return( list(bseqs=bseqs,nr=nr,nm=0) )
  }

  ## Determinar diferències respecte la seqüència (haplotip) màster
  # bseqs[which.max(nr)] permet assignar la màster com l'haplotip amb major nº de reads
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[which.max(nr)])
  # Guarda el nº de diferències obtingudes dels aliniaments
  nm <- nmismatch(psa)

  ## Elimina les seqs amb massa diferències respecte la màster
  # Subset de seqüències d'haplotip que tenen un nombre de diferències menor al màxim permès
  bseqs <- bseqs[nm<=max.difs]
  # També actualitza el nº de reads només dels haplotips amb menys diferències del màxim permès
  nr <- nr[nm<=max.difs]
  # Elimina les dades dels aliniaments amb major nº de diferències del permès
  nm <- nm[nm<=max.difs]

  # Si després d'eliminar els haplotips amb múltiples diferències només en queda un, realitza el procés
  # del principi per guardar aquella seqüència en el fitxer FASTA
  if(length(bseqs)==1)
  { bnts <- strsplit(bseqs[1],split="")[[1]]
  bnts <- bnts[bnts!="-"]
  bseqs <- paste(bnts,collapse="")
  names(bseqs) <- paste("Hpl.0.0001",nr,100,sep="|")
  writeXStringSet(DNAStringSet(bseqs),flnm)
  return( list(bseqs=bseqs,nr=nr,nm=0) )
  }

  ## Ordenar per nombre de diferències
  # Ordena el nombre de mutacions respecte la màster en ordre ascendent
  o <- order(nm)
  # Ordena les seqs segons el seu nº de mutacions respecte la màster
  bseqs <- bseqs[o]
  # Ordena també les freqüències (reads) segons el nº de mutacions de la seqüència
  nr <- nr[o]
  # Ordena la variable amb el nº de mismatches segons les mutacions (ordre ascendent)
  nm <- nm[o]

  ## Numero d'ordre dins de cada nombre de mutacions:
  # Guarda les vegades que apareix cada nº de mismatches
  tnm <- table(nm)
  # length(tnm) calcula el nº d'agrupacions de la taula tnm, és a dir el nº màxim de mutacions que s'han trobat
  # 1:tnm[i] s'aplica sobre els valors de 1 fins al total de mutacions trobades. Per cada nº de mutacions, retorna un
  # conjunt de nombres que van de l'1 al total de vegades que s'ha donat aquell nº de mutacions
  # 'unlist()' concatena tots els valors per tots el nº de mutacions
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))

  ## Ordena per freqüència descendent dins de cada nombre de mutacions:
  # as.integer(names(tnm)) retorna els valors de 1 fins al total de mutacions trobades
  for(i in as.integer(names(tnm)))
  { # Indexs de les seqüències que presenten aquell nº de mutacions (i) respecte la màster
    idx <- which(nm==i)
    # Pels valors de freqüència pertanyents a les seqs amb 'i' mutacions, retorna la seva posició
    # assignada en ordenar-les en ordre descendent
    o <- order(nr[idx],decreasing=TRUE)
    # Ordena les seqüències amb 'i' mutacions segons freq en ordre descendent
    bseqs[idx] <- bseqs[idx[o]]
    # Ordena les freqüències de les seqs amb 'i' mutacions en ordre descendent
    nr[idx] <- nr[idx[o]]
  }

  ## Calcula la freqüència relativa dels haplotips coincidents en FW i RV
  frq <- round(nr/sum(nr)*100,2)

  ## Nom complet per cada haplotip
  # Defineix el nom que rebrà cadascun dels haplotips coincidents en les cadenes
  # nm= mutacions respecte seq màster
  nms <- paste("Hpl", nm, sprintf("%04d",isq), sep=".")

  ## Genera la capçalera FASTA amb el nom de l'haplotip, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")


  ### Eliminar columnes de tot gaps
  # Per cadascuna de les seqs dels haplotips coincidents, aplica la funció 'strsplit()' per
  # separar-les per nucleòtids
  nuc.mat <- t(sapply(bseqs,function(s) strsplit(s,split="")[[1]]))
  # Guarda els gaps de totes les seqüències
  fl <- apply(nuc.mat,2,function(st) all(st=="-"))
  # Si es troba algun gap,
  if(sum(fl))
  { # Guarda només els nucleòtids sense gaps
    nuc.mat <- nuc.mat[,!fl]
    # Torna a ensamblar els nucleòtids per generar les seqs de DNA
    bseqs <- apply(nuc.mat,1,paste,collapse="")
  }

  ## Genera el fitxer FASTA amb els haplotips que han coincidit i sense gaps
  writeXStringSet(DNAStringSet(bseqs),flnm)
  # Retorna una llista amb les seqs dels haplotips coincidents, el nº de reads i el nº de mutacions
  result <- list(bseqs=bseqs,nr=nr,nm=nm)
return(result)
}
