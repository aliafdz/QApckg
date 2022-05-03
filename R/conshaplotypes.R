############################################################
###    RAW CONSENSUS HAPLOTYPES by MULTIPLE ALIGNMENT    ###
############################################################


### Funció per separar els noms dels haplotips amb el recompte de reads i la freqüència
split.fasta.names <- function(seqs,mnr=2)
{ # Guarda els noms assignats a les seqüències de l'objecte DNAStringSet de l'argument
  IDstr <- names(seqs)
  # Guarda el total de noms assignats del pas anterior, és a dir el total d'haplotips
  n <- length(IDstr)
  # Guarda un vector buit amb tantes entrades com haplotips
  nms <- character(length=n)
  # Guarda una matriu buida amb n files (1 per haplotip) i dues columnes, definides amb
  # 'colnames()', on es guardaran en nº de reads de l'haplotip i la seva freqüència
  sts <- matrix(0,nrow=n,ncol=2)
  colnames(sts) <- c("nseqs","pct1")
  # Iteració sobre el total d'haplotips del fitxer avaluat
  for(j in 1:n)
  { # Separa el nom de l'identificador de l'haplotip, que inclou 3 elements:
    # 1) Nº de mutacions i ordre de l'haplotip dins del grup amb aquestes mutacions
    # 2) Nº de reads de l'haplotip
    # 3) Freq relativa dins de la mostra
    strs <- strsplit(IDstr[j],split="\\|")[[1]]
    # Guarda a la taula nms el primer element de l'identificador de l'haplotip corresponent
    nms[j] <- strs[1]
    # Guarda a la 1a i 2a columna de la taula sts el 2n i 3r elements de l'identificador
    sts[j,] <- as.numeric(strs[2:3])
  }
  # Guarda en un data frame els 3 elements de l'identificador de tots els haplotips separats per
  # columnes
  # NOTA: es podria simplificar per tal de guardar les variables directament a partir del bucle?
  IDs <- data.frame(ID=nms,sts,stringsAsFactors=FALSE)
  # Guarda el nº de files del data frame (total d'haplotips)
  nall <- nrow(IDs)
  # Realitza el sumatori de la columna amb el nº de reads per haplotip
  tnr <- sum(IDs$nseqs)
  # Guarda el nº de reads majors al valor definit de minim nombre de reads
  flags <- IDs$nseqs >= mnr
  # Retorna una llista amb un subset dels haplotips que presentin un nº de reads major al mínim definit,
  # les seqüències dels haplotips que compleixen la condició, el nº de total d'haplotips i el sumatori de reads
  return(list(IDs=IDs[flags,],seqs=seqs[flags],nall=nall,tnr=tnr))
}


### Funció per llegir els .fna amb els amplicons ja classificats i sense primers
read.ampl.seqs <- function(flnm,mnr=2)
{ # Guarda les seqüències del fitxer fasta de l'argument obtingudes amb la funció 'readDNAStringSet()'
  seqs <- as.character(readDNAStringSet(flnm))
  # Aplica la funció definida anteriorment per obtenir tots els noms dels haplotips amb un nº de reads
  # major al mínim definit i les seves dades
  return(split.fasta.names(seqs,mnr)) # NOTA: es podrien fusionar les 2 funcions?
}


###  Funció per a l'aliniament múltiple amb muscle de les seqs
###  Les seqüècies s'entren i es tornen com a DNAStringSet
# Guarda la ruta de l'executable muscle (aquest es posarà al fitxer global)
muscle <- "C:\\Muscle\\muscle3.8.31_i86win32.exe"
# Guarda el nom del fitxer d'opcions de muscle que es guardarà a l'entorn global del projecte
muscle.cl.opts <- c("-log muscle.log")

doMuscle <- function(seqs)
{ # Genera un fitxer al directori tmp per introduir les seqüències a alinear
  tmp.file <- file.path(tempDir,"muscleInFile.fna")
  # Genera el fitxer al directori tmp que resultarà de l'alineament amb muscle
  res.file <- file.path(tempDir,"muscleOutFile.fna")
  # Si ja existeix l'arxiu de resultats, s'elimina per generar-lo de nou
  if(file.exists(res.file)) file.remove(res.file)
  # Aplica la funció d'abans per escriure el fitxer fasta d'entrada amb les seqüències de l'argument
  writeXStringSet( seqs, tmp.file )
  # Genera els fitxers in i out per a executar muscle
  in.file <- paste("-in ",tmp.file,sep="")
  out.file <- paste("-out ",res.file,sep="")
  # Genera la comanda per executar muscle amb el fitxer .exe, els arxius in i out i el fitxer d'opcions
  command <- paste(muscle,in.file,out.file,
                   paste(muscle.cl.opts,collapse=" "),sep=" ")
  # Empra la funció 'system()' per executar la comanda indicada abans i fer l'aliniament múltiple
  system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
         ignore.stdout=FALSE,invisible=TRUE)
  # Retorna el fitxer de resultats com a objecte DNAStringSet gràcies a la funció de llegir arxius fasta
  if( file.exists(res.file) )
    return( readDNAStringSet(res.file) )
  return(NULL)
}


### Alinear les distribucions d'haplotips de dues poblacions
# Els arguments de la funció corresponen al nº de reads dels haplotips FW o RV (nA o nB)
# i les seqüències en si (seqsA i seqsB)
PopsAlgnHist <- function(nA,seqsA,nB,seqsB)
{ # Relaciona les seqüències dels haplotips FW (A) i RV (B) amb el seu nº de reads
  names(nA) <- seqsA
  names(nB) <- seqsB
  # Combina les dades dels dos conjunts de seqüències, eliminant aquelles que estiguin duplicades entre les cadenes
  nms <- union(seqsA,seqsB) # tots: A + B
  # Guarda la longitud (nº de seqs) de la nova combinació -> No correspon a la suma de seqsA i seqsB perquè elimina els duplicats
  nb <- length(nms)
  # Guarda un vector buit de longitud igual al nº de seqüències totals, per guardar les coincidències amb FW
  pA <- integer(nb)
  # Indica els indexs d'haplotips FW dins del conjunt global
  idx <- which(nms %in% seqsA)
  # Tenint en compte els indexs obtinguts, assigna el nº de reads dels haplotips FW en les seves corresponents posicions
  pA[idx] <- nA[nms[idx]]
  # Guarda un vector buit de longitud igual al nº de seqüències totals, per guardar les coincidències amb RV
  pB <- integer(nb)
  # Indica els indexs d'haplotips RV dins del conjunt global
  idx <- which(nms %in% seqsB)
  # Tenint en compte els indexs obtinguts, assigna el nº de reads dels haplotips RV en les seves corresponents posicions
  # Si ara es compara pA i pB, les posicions on hi ha més de 0 reads en ambdòs vectors corresponen als haplotips coincidents!
  pB[idx] <- nB[nms[idx]]
  # Genera una llista amb el nº de seqs dels haplotips FW i RV, i també el conjunt global amb totes les seqs
  # IMPORTANT! En aquest pas encara no s'ha fet la intersecció per aliniament
  list(pA=pA,pB=pB,Hpl=nms)
}


### Confrontació d'haplotips FW i RV
# Els arguments de la funció corresponen al nº de seqüències dels haplotips FW o RV (nA o nB)
# i les seqüències en si (seqsA i seqsB)
Intersect.FWRV <- function(nA,seqsA,nB,seqsB)
{ ### Aliniament de distribucions
  # Aplica la funció d'abans per obtenir el conjunt de totes les seqs d'haplotips FW i RV i el seu nº de reads
  lst <- PopsAlgnHist(nA,seqsA,nB,seqsB)

  ### Haplotips comuns i solapament global
  # Indica en quins casos coincideixen les seqs dels haplotips FW i RV, ja que seran els indexs on el nº de reads
  # per A i per B sigui major a 0
  fl <- lst$pA>0 & lst$pB>0
  # sum(lst$pA+lst$pB) correspon al total de reads d'ambdues cadenes que han entrat a la intersecció després del filtre per abundància
  # sum(lst$pA[fl]+lst$pB[fl]) correspon al total de reads que coincideixen en FW i RV
  # Per tant, el resultat de ov.a és la fracció de reads comuns entres les cadenes respecte el total
  # Es considera que el sumatori de reads en ambdues cadenes correspon a la cobertura de l'amplicó
  ov.a <- sum(lst$pA[fl]+lst$pB[fl])/sum(lst$pA+lst$pB)

  ### Freqüències i solapament per intersecció
  # Calcula la freq relativa (nº de reads d'un haplotip entre el total dels que entraven a la intersecció d'aquella cadena)
  # per aquells haplotips coincidents en ambdues cadenes
  pFW <- (lst$pA/sum(lst$pA))[fl]
  pRV <- (lst$pB/sum(lst$pB))[fl]
  # Indica el valor mínim de freq al comparar cada parella d'haplotips coincidents -> freq mínima d'intersecció
  p <- pmin(pFW,pRV)
  # Calcula el sumatori de tots els valors mínims de freq calculats entre les parelles coincidents
  # Aquest valor correspondrà al percentatge de superfície que solapa entre ambdues cadenes
  ov.i <- sum(p)
  # Ara renormalitza les dades dividint cada valor mínim de freq entre el sumatori de totes
  p <- p/sum(p)
  ## Retorna una llista que inclou, en ordre:
  # 1) Freq relatives dels haplotips normalitzades
  # 2) Seqüències dels haplotips coincidents en FW i RV
  # 3) Total de superfície que solapa entre ambdues cadenes
  # 4) Fracció de reads comuns entre les cadenes respecte el total
  # 5) Vector amb el nº de seqs dels haplotips FW que entraven a la intersecció (després de filtrar per abundància)
  # 6) Vector amb el nº de seqs dels haplotips RV que entraven a la intersecció
  list(p=p,seqs=lst$Hpl[fl],ov.i=ov.i,ov.a=ov.a,pA=lst$pA,pB=lst$pB)
}


### Funció per representar les distribucions aliniades dels haplotips
PlotHplHistos <- function(tt,pA,pB,p)
{ # Divideix el nº de seqs de cada haplotip de la cadena FW entre el total dels que entraven a la intersecció
  pA <- pA/sum(pA)
  # Divideix el nº de seqs de cada haplotip de la cadena RV entre el total
  pB <- pB/sum(pB)
  # Defineix el límit de l'eix Y en funció de les freqüències calculades
  ymx <- max(c(pA,pB))
  # Gràfic de barres per representar els haplotips de la cadena FW i la seva freq relativa
  barplot(pA,ylim=c(0,ymx)); abline(h=0)
  title(main="FW strand haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  # Gràfic de barres per representar els haplotips de la cadena RV i la seva freq relativa
  barplot(pB,ylim=c(0,ymx)); abline(h=0)
  title(main="RV strand haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  # Gràfic de barres per representar els haplotips coincidents en ambdues cadenes i la seva freq
  barplot(p); abline(h=0)
  title(main="Intersected haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
}



### Dóna nom als haplotips en funció del nombre de mutacions respecte la màster
### i de la seva frequencia poblacional, i salva les seqüències en un fitxer fasta.
# Els arguments corresponen al nom del fitxer on es desaran els resultats, les seqs
# dels haplotips coincidents en ambdues cadenes, el sumatori de reads per haplotip coincident
# i el nº màxim de diferències permeses
SaveHaplotypes <- function(flnm,bseqs,nr,max.difs=250)  #,seq0)
{
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
    # Guarda la seq en un fitxer fasta en la ruta indicada en l'argument de la funció
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
  # del principi per guardar aquella seqüència en el fitxer fasta
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

  ## Genera la capçalera fasta amb el nom de l'haplotip, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ###  Afegir-hr la RefSeq -> no s'executa
  #bseqs <- c(seq0,bseqs)

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

  ## Genera el fitxer fasta amb els haplotips que han coincidit i sense gaps
  writeXStringSet(DNAStringSet(bseqs),flnm)
  # Retorna una llista amb les seqs dels haplotips coincidents, el nº de reads i el nº de mutacions
  list(bseqs=bseqs,nr=nr,nm=nm)
}


#-------------------------------------------------------------------#

#' ConsHaplotypes
#'
#' Consensus haplotypes by multiple alignment
#'
#' @param thr Minimum threshold of abundance to enter in multiple alignment.
#' @param min.seq.len Minimum length of sequences to enter in intersection.
#' @param min.rd Minimum nº of reads required in fasta files.
#' @param pm.res List... resultat de la funció anterior
#'
#' @import Biostrings
#' @import ape
#' @examples
#' min.seq.len <- 150
#' min.rd <-   1
#' thr <- 0.2
#'
#' @export

conshaplotypes <- function(samples, primers, thr, min.seq.len,min.rd,pm.res){
## No cal introduir els arxius fasta de la carpeta trim perquè es treuen de la taula fileTable de la funció anterior.
# En lloc de carregar el RData es carreguen les taules fileTable i poolTable de la variable de la funció anterior.

  FlTbl<- pm.res$fileTable
  PoolTbl<- pm.res$poolTable


# Retorna els indexs de la taula FlTbl derivats d'ordenar les mostres en funció de l'ID del
# pacient, de l'amplicó (regió) avaluat i la cadena forward o reverse
o <- order(FlTbl$Pat.ID,FlTbl$Ampl.Nm,FlTbl$Str)
# Reordena les entrades de la taula en funció dels indexs anteriors, de manera que s'agrupen
# les entrades corresponents al mateix pacient i en segon ordre a la mateixa regió avaluada
FlTbl <- FlTbl[o,]

### Guarda els noms dels fitxers fasta inclosos al directori trim i els separa en funció
# de la cadena asignada als reads
in.files <- file.path(trimDir,FlTbl$File.Name)
idx.fw <- which(FlTbl$Str=="fw")
idx.rv <- which(FlTbl$Str=="rv")

# Guarda els noms dels fitxers resultants d'aquest script, que correspondran a la concatenació
# del terme MACHpl02, l'ID del pacient i el nom de la regió avaluada, indicats de manera que només
# s'obtingui un fitxer .fna per a cada regió del mateix pacient (per això el condicional fw)
out.flnms <- paste("MACHpl02",FlTbl$Pat.ID[idx.fw],
                   FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
# Guarda la ruta on es desaran els fitxers .fna, a la carpeta MACH
out.flnms <- file.path(mach.Dir,out.flnms)

# Guarda els noms d'uns altres fitxers resultants, que correspondran a la concatenació
# del terme MAfwrv, l'ID del pacient i el nom de la regió avaluada, indicats de manera que només
# s'obtingui un fitxer .fna per a cada regió del mateix pacient (per això el condicional fw)
ma.flnms <- paste("MAfwrv",FlTbl$Pat.ID[idx.fw],
                  FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
# Guarda la ruta on es desaran els fitxers .fna, a la carpeta MACH
ma.flnms <- file.path(mach.Dir,ma.flnms)

# Guarda la meitat de la longitud de fitxers presents a la carpeta trim del fitxer previ
n <- length(in.files)/2

# Genera un data frame inicialment buit, amb n files (2 entrades per pacient) i 5 columnes
# que inclouran els resultats obtinguts en pasos posteriors per a les cadenes fw
rdf.fw <- data.frame(fw.all=integer(n),fw.lowf=integer(n),
                     fw.in=integer(n),fw.unq=integer(n),fw.com=integer(n))
# Genera un altre data frame també buit, amb n files i 5 columnes per incloure resultats
# de les cadenes rv
rdf.rv <- data.frame(rv.all=integer(n),rv.lowf=integer(n),
                     rv.in=integer(n),rv.unq=integer(n),rv.com=integer(n))
# Genera un altre data frame buit, amb n files i 6 columnes, per incloure els resultats globals
rdf.gbl <- data.frame(all=integer(n),lowf=integer(n),unq=integer(n),
                      ovrlp=numeric(n),common=numeric(n),Fn.rd=integer(n))

# Genera el primer fitxer pdf on s'inclouran diverses representacions de resultats generats en el bucle for
# Per cada regió amplificada de cada pacient, es representen els haplotips de cada cadena i els que han
# coincidit entre les dues, amb les corresponents freqüències
pdf(file.path(repDir,"MA.Intersects.plots.pdf"),paper="a4",
    width=7,height=11)
par(mfrow=c(3,1))

# Bucle for per iterar sobre totes les mostres (2 per pacient)
for(i in 1:n)
{ # Si no existeix el fitxer de la mostra avaluada a la carpeta trim, torna a començar la iteració següent
  if(!file.exists(in.files[idx.fw[i]])) next
  if(!file.exists(in.files[idx.rv[i]])) next

  # Aplica la funció 'ReadAmplSeqs()' del paquet QSutils que permet obtenir els haplotips i les
  # seves freqüències a partir d'un fitxer fasta indicat.
  lstf <- ReadAmplSeqs(in.files[idx.fw[i]])
  # Guarda les seqüències de longitud major al mínim permès i filtra el nº de reads
  # eliminant els de les seqüències que presenten longitud menor al mínim permès
  lstf$nr <- lstf$nr[nchar(lstf$hseqs) >= min.seq.len]

  # Aplica la funció 'GetQSData()' del paquet QSutils per llegir els amplicons del
  # fitxer amb els valores d'abundància, filtrar per mínima abundància i ordenar
  # segons les mutacions.
  s1 <- GetQSData(in.files[idx.fw[i]],a.cut)
  # Guarda les seqüències de longitud major al mínim permès
  ffilt <- nchar(s1$seqs) >= min.seq.len
  # Filtra les seqüències per eliminar les que presenten longitud menor al mínim permès
  seqsf <- s1$seqs[ffilt]
  nrf <- s1$nr[ffilt]

  # Calcula la freqüència relativa de cada amplicó filtrat per abundància respecte el total
  # de reads del fitxer inicial.
  frq_lst <- round(nrf/sum(lstf$nr)*100,2)
  # Donat que la funció del paquet QSutils retorna els noms dels haplotips diferent al desitjat, es canvien els
  # noms substituint els caràcters '_' per '.', i generant la capçalera del fitxer fasta amb nom de l'haplotip,
  # nombre de reads i freq relativa. També es substitueix el terme Hpl per HplFw amb la funció 'gsub()'.
  names(seqsf) <- paste(gsub("Hpl","HplFw",gsub('_','.',names(seqsf))),nrf,frq_lst,sep="|")

  ## Guarda a la taula de resultats per a cadenes fw:
  # El total de reads de la mostra avaluada
  rdf.fw$fw.all[i]  <- sum(lstf$nr)
  # El nº de reads amb baixa freqüència (menor al mínim permès)
  rdf.fw$fw.lowf[i] <- sum(lstf$nr)-sum(nrf)

  ## Aplica el mateix procés per a les cadenes classificades reverse de la mostra avaluada

  # Aplica la funció 'ReadAmplSeqs()' del paquet QSutils que permet obtenir els haplotips i les
  # seves freqüències a partir d'un fitxer fasta indicat.
  lstr <- ReadAmplSeqs(in.files[idx.rv[i]])
  # Guarda les seqüències de longitud major al mínim permès i filtra el nº de reads
  # eliminant els de les seqüències que presenten longitud menor al mínim permès
  lstr$nr <- lstr$nr[nchar(lstr$hseqs) >= min.seq.len]

  # Aplica la funció 'GetQSData()' del paquet QSutils per llegir els amplicons del
  # fitxer amb els valores d'abundància, filtrar per mínima abundància i ordenar
  # segons les mutacions.
  s2 <- GetQSData(in.files[idx.rv[i]],a.cut)
  # Guarda les seqüències de longitud major al mínim permès
  rfilt <- nchar(s2$seqs) >= min.seq.len
  # Filtra les seqüències per eliminar les que presenten longitud menor al mínim permès
  seqsr <- s2$seqs[rfilt]
  nrr <- s2$nr[rfilt]

  # Calcula la freqüència relativa de cada amplicó filtrat per abundància respecte el total
  # de reads del fitxer inicial.
  frq_lst <- round(nrr/sum(lstr$nr)*100,2)
  # Donat que la funció del paquet QSutils retorna els noms dels haplotips diferent al desitjat, es canvien els
  # noms substituint els caràcters '_' per '.', i generant la capçalera del fitxer fasta amb nom de l'haplotip,
  # nombre de reads i freq relativa. També es substitueix el terme Hpl per HplFw amb la funció 'gsub()'
  names(seqsr) <- paste(gsub("Hpl","HplRv",gsub('_','.',names(seqsr))),nrr,frq_lst,sep="|")

  ## Guarda a la taula de resultats per a cadenes rv:
  # El total de reads de la mostra avaluada
  rdf.rv$rv.all[i]  <- sum(lstr$nr)
  # El nº de reads amb baixa freqüència (menor al mínim permès)
  rdf.rv$rv.lowf[i] <- sum(lstr$nr)-sum(nrr)

  # Guarda en un vector les seqs de les dues cadenes
  aseqs <- c(seqsf,seqsr)


  ###  Aliniament múltiple per muscle dels haplotips fw i rv
  # Guarda el resultat de l'aliniament múltiple realitzat amb la funció definida al principi
  seqs <- doMuscle(DNAStringSet(aseqs))
  # Copia el fitxer resultant de l'aliniament múltiple al directori MACH per guardar l'aliniament
  # dels haplotips de la mostra avaluada
  file.copy(file.path(tempDir,"muscleOutFile.fna"),ma.flnms[i],
            overwrite=TRUE)

  # Aplica la funció per obtenir els noms dels haplotips i les seves dades
  lst <- split.fasta.names(as.character(seqs))
  # Guarda el vector que inclou el nº de seqüències dels haplotips
  nr <- lst$IDs$nseqs
  # Guarda les seqüències dels haplotips amb més reads del mínim permès (per defecte a la funció mnr=2)
  seqs <- lst$seqs

  ## Separar seq FW i RV ja aliniades
  # Després de l'aliniament múltiple de totes les seqüències de la mostra, separa en 2 variables les que
  # corresponen a cadena forward o reverse
  ifw <- grep("^HplFw",names(seqs))
  irv <- grep("^HplRv",names(seqs))

  # Aplica una altra funció definida al principi per calcular la intersecció entre els haplotips
  lst <- Intersect.FWRV(nr[ifw],seqs[ifw],nr[irv],seqs[irv])
  # Indica en quins casos coincideixen les seqs dels haplotips FW i RV, ja que seran els indexs on el nº de reads
  # per A i per B sigui major a 0
  fl <- lst$pA>0 & lst$pB>0  # Nota: Aquesta variable es podria obtenir directament de la funció

  ## Guarda a les taules de resultats de cadenes FW i RV:
  # El nº total de reads dels haplotips que entraven a la intersecció
  rdf.fw$fw.in[i]  <- sum(lst$pA)
  rdf.rv$rv.in[i]  <- sum(lst$pB)
  # El nº de reads dels haplotips coincidents entre ambdues cadenes
  rdf.fw$fw.com[i] <- sum(lst$pA[fl])
  rdf.rv$rv.com[i] <- sum(lst$pB[fl])
  # El nº de reads dels haplotips no coincidents (únics d'una cadena)
  rdf.fw$fw.unq[i] <- sum(lst$pA[!fl])
  rdf.rv$rv.unq[i] <- sum(lst$pB[!fl])

  ## Guarda a la taula de resultats globals (que en aquest punt encara està buida):
  # El total de reads FW+RV de la mostra avaluada
  rdf.gbl$all[i]    <- rdf.fw$fw.all[i]+rdf.rv$rv.all[i]
  # Sumatori de tots els valors mínims de freq calculats entre les parelles coincidents: intersecció entre reads
  # Percentatge de superfície que solapa entre ambdues cadenes
  rdf.gbl$ovrlp[i]  <- round(lst$ov.i*100,2)
  # Percentatge de reads dels haplotips coincidents en FW i RV: % de reads comuns respecte el total
  # Coincideix amb el resultat de dividir els reads que han coincidit (suma de fw.com i rv.com) entre els reads
  # que han entrat en la intersecció (suma de fw.in i rv.in)
  rdf.gbl$common[i] <- round(lst$ov.a*100,2)
  # Sumatori del nº de reads totals dels haplotips coincidents d'ambdues cadenes: reads finals
  rdf.gbl$Fn.rd[i]  <- rdf.fw$fw.com[i]+rdf.rv$rv.com[i]
  # Sumatori d reads dels haplotips amb baixa freqüència d'ambues cadenes
  rdf.gbl$lowf[i]   <- rdf.fw$fw.lowf[i] + rdf.rv$rv.lowf[i]
  # Sumatori dels reads únics d'una cadena de DNA (no coincidents)
  rdf.gbl$unq[i]    <- rdf.fw$fw.unq[i] + rdf.rv$rv.unq[i]

  ### Si no es detecta superposició entre els haplotips d'ambdues cadenes, salta a la seqüent iteració
  if(sum(fl)==0)
  { cat("\n--------------------------------------------------\n")
    next
  }

  ### Gràfic de barres dels haplotips alineats
  # Concatenació de l'ID del pacient i la regió avaluada en la iteració
  tt <- paste(FlTbl$Pat.ID[idx.fw[i]]," - ",FlTbl$Ampl.Nm[idx.fw[i]],sep="")
  # Suma el nº de seqs dels haplotips FW i RV independentment de si coincideixen o no
  p <- lst$pA+lst$pB
  # Substitueix el nº de reads d'aquells haplotips no coincidents per 0: Només queda el sumatori de reads
  # per cada haplotip que ha coincidit en ambdues cadenes
  p[ lst$pA==0 | lst$pB==0 ] <- 0
  # Calcula la freq relativa del nº de reads de cada haplotip entre el total dels que han coincidit
  p <- p/sum(p)
  # Aplica la funció local per representar els haplotips de cada cadena i els aliniats amb les seves freq
  PlotHplHistos(tt,lst$pA,lst$pB,p)

  ### Suma de reads de cada haplotip coincident en FW i RV
  rds <- lst$pA[fl]+lst$pB[fl]

  ### Aplica la funció local del principi per guardar els haplotips coincidents en ambdues cadenes
  # amb les seves freqüències en els fitxers MACHpl02
  SaveHaplotypes(out.flnms[i],lst$seqs,rds)  #,seq0)
} # Fi del loop sobre totes les mostres (2 per pacient)

dev.off()

# Genera el fitxer .txt de resultats de les interseccions
sink(file=file.path(repDir,"MA.Intersects-SummRprt.txt"))

cat("\n   FW + RV HAPLOTYPES INTERSECTIONS")
cat("\n======================================\n")
cat("\nCutting FW and RV at ",a.cut,"% followed by haplotypes intersection.\n",
    sep="")
cat("\nFrequencies as sum of reads FW+RV.\n\n")

## Resultats per total de reads
cat("\nSUMMARY RESULTS BY READS NUMBER")
cat("\n===============================\n\n")
# A partir de la taula FlTbl (provinent de l'arxiu RData del pas anterior del pipeline), guarda
# les entrades corresponents a la cadena forward
fl.fw <- FlTbl$Str=="fw"
# Genera un data frame amb els ID dels pacients, la regió amplificada, i els resultats obtinguts
# al llarg de la intersecció d'haplotips per la cadena FW
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                   rdf.fw,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
# Genera el mateix data frame d'abans pels haplotips de la cadena RV per incloure-la al fitxer
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                   rdf.rv,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
# Genera un altre data frame per incloure els resultats globals d'ambdues cadenes (FW+RV)
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                   rdf.gbl,stringsAsFactors=FALSE)
print(frdf)

## Mostra els resultats d'aplicar el sumatori sobre les columnes de resultats
# de cadascun dels 3 data frames anteriors
cat("\n\nTotal counts:\n\n")
tots.fw <- apply(rdf.fw,2,sum)
print(tots.fw)
cat("\n")
tots.rv <- apply(rdf.rv,2,sum)
print(tots.rv)
cat("\n")
tots.gbl <- apply(rdf.gbl[,c(1:3,6)],2,sum)
print(tots.gbl)
cat("\n")

## Aplica els percentatges sobre els resultats dels sumatoris anterior per FW, RV y global
cat("\nPercentage:\n\n")
ptts <- round(tots.fw/tots.fw[3]*100,2) ## Nota: revisar aquests resultats (% sobre la columna 3)
print(ptts)
cat("\n")
ptts <- round(tots.rv/tots.rv[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.gbl/tots.gbl[1]*100,2)
print(ptts)
cat("\n")

sink()

## Guarda l'últim data frame amb els resultats globals en un arxiu .RData
save(frdf,file=file.path(repDir,"IntersectedReads.RData"))

## Genera el pdf que es guardarà a la carpeta reports on s'inclouran els gràfics
# de resultats globals de la intersecció d'haplotips
pdf.flnm2 <- "IntersectBarplots.pdf"
pdf(file.path(repDir,pdf.flnm2),paper="a4r",width=10,height=6)
par(mar=c(7,4,4,2)+.1)

## Genera la paleta de colors
pal <- cls[1:2] # cls <- brewer.pal(8,"Dark2") definit a l'arxiu 'seqsanalfns.v4.5.R'
# Guarda els noms de les mostres a partir de la concatenació de l'ID del pacient i
# la regió avaluada (amb el condicional FW per no repetir els noms dos cops)
bnms <- paste(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw])
# Genera una taula amb les columnes corresponents als valors de reads comuns en FW i RV
vals <- cbind(rdf.fw$fw.com,rdf.rv$rv.com)
# Defineix el límit màxim de l'eix Y a partir del sumatori de reads comuns (aplicat sobre files)
ymx <- max(rowSums(vals))*1.2
# Gràfic de barres on es representen els reads comuns en ambdues cadenes en funció de la mostra
# names.arg correspon a un vector on s'indiquen els noms a representar sota les barres del gràfic
barplot(t(vals),col=pal,ylim=c(0,ymx),names.arg=bnms,
        las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("fw","rv"),cex=0.8)
title(main="Intersected reads")

## Genera una altra taula on s'inclouen les columnes de reads de baixa freqüència i únics d'una cadena
# per als haplotips FW i RV, fent el sumatori per files (per mostra)
vals <- cbind(rowSums(rdf.fw[,c("fw.lowf","fw.unq")]),
              rowSums(rdf.rv[,c("rv.lowf","rv.unq")]))
# Defineix el límit màxim de l'eix Y a partir del sumatori aplicat sobre files de la taula anterior
ymx <- max(rowSums(vals))*1.2
# Gràfic de barres on es representen els reads filtrats totals (per abundància o per no trobar intersecció)
# en funció de la mostra
barplot(t(vals),col=pal,ylim=c(0,ymx),names.arg=bnms,
        las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("fw","rv"),cex=0.8)
title(main="Filtered out reads")

## Defineix la nova paleta de colors
pal <- cls[c(1,3,2,4)]
# Guarda en una matriu les columnes amb els valors de reads finals, únics de cadena i de baixa freq
# (total per ambdues cadenes)
mbp <- data.matrix(frdf[,c(8,5,4)])
# Guarda els noms de les mostres a partir de la concatenació de l'ID del pacient i
# la regió avaluada separats per punt (aquest cop els agafa de la taula de resultats)
mbp.nms <- paste(frdf$Pat.ID,frdf$Ampl.Nm,sep=".")
omar <- par(mar=c(6,3.5,3,2)) # Defineix els marges de pàgina
# Gràfic que representa, per cada mostra, una barra indicant en diferents colors els reads
# comuns (els finals després de la intersecció), els únics de cadena i els de baixa freq
barplot(t(mbp),col=pal,ylim=c(0,max(rowSums(mbp))*1.2),names.arg=mbp.nms,
        las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("Common","Unique","Low freq."),
       cex=0.8)

## Genera un altre gràfic de barres per representar els percentatges de reads
# comuns per cada mostra (rendiment de la intersecció)
par(mar=omar)
barplot(frdf$common,col="lavender",ylim=c(0,100),names.arg=mbp.nms,
        las=2,cex.names=0.6,cex.axis=0.8)
title(main="FW + RV intersection yield",line=1)

### Carregar PoolTbl i Fltbl
# Carrega el fitxer RData generat en el pas anterior del pipeline, on s'han eliminat els primers
# de totes les seqüències les quals es classifiquen en forward o reverse
# Nota: això ja es fa al principi del script
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))

## Guarda els noms concatenants l'ID dels pacients i la regió amplificada separats per punt
all.nms <- paste(FlTbl$Pat.ID,FlTbl$Ampl.Nm,sep=".")
# Realitza el sumatori dels reads FW i RV de cadascuna de les mostres
# Nota: es podria obtenir de la taula rdf.gbl (columna all)
trds <- tapply(FlTbl$Reads,all.nms,sum)
# Actualitza el vector amb els noms per eliminar els duplicats, ja que al fer el sumatori
# pasa a tenir una única entrada per mostra (i no 2 quan tenia FW i RV)
all.nms <- names(trds)
# Guarda la columna de reads finals per mostra (després de la intersecció)
frds <- frdf$Fn.rd
# Els noms de la columna de reads finals seran de nou la concatenació de ID de pacient i regió
# Nota: el paste() seria equivalent a posar la variable all.nms
names(frds) <- paste(frdf$Pat.ID,frdf$Ampl.Nm,sep=".")
# Genera una matriu buida amb tantes files com noms de mostres s'han generat, i 2 columnes
ytbl <- matrix(0,nrow=length(all.nms),ncol=2)
# Els noms de les files de la matriu correspondran als noms de les mostres
rownames(ytbl) <- all.nms
# Els noms de les columnes corresponen als reads totals i els que han quedat
# després de la intersecció
colnames(ytbl) <- c("All","Passed")
# Guarda els valors de la columna de reads totals
ytbl[,1] <- trds
# Assigna a cada mostra els valors de la columna de reads finals (que han passat)
ytbl[names(frds),2] <- frds
omar <- par(mar=c(6,3.5,3,2)) # Marges de pàgina
# Gràfic de barres que representa, per cada mostra, una barra indicant els reads totals
# i una altra els reads que han passat
barplot(t(ytbl),beside=TRUE,col=cls[1:2],ylim=c(0,max(ytbl)*1.2),las=2,
        names.arg=rownames(ytbl),cex.names=0.6,cex.axis=0.8)
abline(h=0)
legend("top",horiz=TRUE,fill=cls,legend=colnames(ytbl),cex=0.8)

## Últim gràfic per representar el rendiment de la intersecció, dividint els
# reads que han passat respecte els totals (els inicials del pas)
par(mar=omar)
barplot(as.vector(ytbl[,2]/ytbl[,1]*100),col="lavender",ylim=c(0,100),
        names.arg=rownames(ytbl),las=2,cex.names=0.6,cex.axis=0.8)
abline(h=0)
title(main="Global yield",line=1)

dev.off()

}
