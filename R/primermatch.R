
primermatch <- function(idx,flnms,pool,PoolTbl){

# foreach(i=1:length(pools)) %do%
#   { # Identifica les mostres que corresponen al pool avaluat
#     idx <- which(samples$Pool.Nm==pools[i])
#
#     # Guarda els fitxers de la carpeta split que contenen reads assignats a algun MID (exclou els fitxers MID0)
#     indexmids <- which(splitfiles!=splitfiles[grep('MID0',splitfiles)])
#     # Actualitza els fitxers a avaluar excloent els que no s'han assignat a MID
#     splitfiles <- splitfiles[indexmids]
#     # Guarda el noms dels fitxers de les mostres (MID) corresponents al pool avaluat
#     flnms <- splitfiles[grep(pools[i],splitfiles)]

    ### Loop sobre les mostres del pool avaluat -> cadascuna amb un MID concret
#    for(j in 1:length(idx))
#    { # Guarda l'index de la mostra avaluada
      j <- idx #[j]

      ### Carrega el fitxer fastq de la carpeta split corresponent al MID de la mostra avaluada
      # Si el fitxer no es troba al directori indicat, executa la següent iteració
      if(! file.exists(file.path(splitDir,flnms[j]))) stop(paste('File',flnms[j],'not found.'))

      # Llegeix el fitxer .fna del MID de la mostra avaluada
      seqs <- readDNAStringSet(file.path(splitDir,flnms[j]))

      ### Per guardar el nº total de reads abans de demultiplexar per primers:
      totalseqs <- length(seqs)

      # Afegeix el total de reads del fitxer del MID concret al total del pool al qual correspon
      PoolTbl[pool,1] <- PoolTbl[pool,1] + length(seqs)

      ### Retorna la posició en la columna Ampl.Nm del primer ID que s'està avaluant
      # És a dir, de la mostra dins del pool que estem avaluant, indica quin és el seu primer ID i el busca
      # al fitxer de primers
      ipr <- grep(samples$Primer.ID[j],primers$Ampl.Nm)

      # Concatena els caràcters d'associació entre el pool que s'està avaluant, el MID i la regió de l'amplicó
      tt <- paste("Pool",pools[i],"  MID",samples$MID[j],
                  "  Ampl",primers$Ampl.Nm[ipr])
      cat("\n\n",tt,"\n",sep="")

      ### Elimina les seqüències del fitxer del MID avaluat que tinguin longitud menor a 180
      seqs <- seqs[width(seqs)>min.len]

      ### Coincidències cadena forward
      # Guarda la seq del primer FW específic de la regió avaluada
      pr.up <- primers$Primer.FW[ ipr ]
      # Busca la seq del primer FW en la regió 5' (posicions 1-100 definides al fitxer de paràmetres)
      # de les seqüències del fitxer .fna del MID avaluat. prmm defineix el màxim de mismatches permesos
      # Guarda totes les coincidències amb les posicions inicial, final, i longitud
      up.matches <- vmatchPattern(pattern=pr.up,
                                  subject=subseq(seqs,start=target.io,end=target.in),
                                  max.mismatch=prmm,fixed=FALSE)
      # Resta 1 a la posició inicial de cerca del primer (per defecte, 1)
      delta <- target.io-1
      # Aplica la funció que està definida al fitxer principal
      # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer FW
      flags <- elementNROWS(up.matches)>=1
      # Variable buida de nombre enter
      nr <- integer()
      # Calcula quantes seqs s'han eliminat per ser menor a la longitud establerta
      shorts <- totalseqs - length(seqs) ####
      trim.len <- 0
      # Genera un fitxer .fna, el nom del qual està format per l'ID del pacient, la regió (amplicó) avaluada,
      # i la consecució PrFW.fna (primer forward), que es guarda a la carpeta trim
      up.flnm <- paste(samples$Patient.ID[j],primers$Ampl.Nm[ipr],"PrFW.fna",
                       sep=".")

      seqs.up <- "" # Coincidències cadena FW
      if(sum(flags)) # Sumatori de seqüències amb 1 o més coincidències amb la seq del primer
      { # 'startIndex()' retorna una llista amb les posicions inicials de les coincidències del patró cercat
        # Cal tenir cura, ja que si no hi ha coincidència la llista guarda un valor NULL
        # Aplicar el condicional flags sobre la llista de posicions inicials permet eliminar els valors NULL
        # Guarda en una variable tots els valors de posició inicial del primer FW sobre les seqüències
        pos <- sapply(startIndex(up.matches)[flags],function(x) x[[1]])+delta

        ### Guarda les seqs amb coincidències del primer FW
        seqs.up <- seqs[flags]
        ### Actualitza el nº de seqs amb les que no presentaven coincidències
        seqs <- seqs[!flags]

      ### Guarda el nº de seqs totals a partir de les quals buscar el primer RV més endavant
      treads <- length(seqs)

        ### Retalla el primer 5'
        # Suma a la posició inicial del primer FW en la seq la seva longitud
        st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])
        # Retalla les seqüències des de la posició on acaba el primer FW fins al final
        seqs.up <- subseq(seqs.up,start=st,end=width(seqs.up))

        ### Localitza el primer per l'altre extrem --> a la mateixa cadena!
        # Guarda la reversa complementària del primer RV de la regió avaluada
        pr.p3 <- as.character(
          reverseComplement(DNAString(primers$Primer.RV[ipr])))
        # Calcula la longitud de les seqüències que tenien coincidència amb el primer FW després de retallar-lo
        # Seria equivalent a dir que es calcula la posició final del primer RV (parlant en sentit 5'-3')
        endp <- width(seqs.up)
        # Resta a la longitud de les seqs la posició final on hem de buscar el primer
        # Com estem parlant del primer RV, que es troba a 3', l'haurem de buscar en les 100 últimes posicions
        fstp <- endp-target.in
        # Calcula la diferència entre la posició on comença el primer RV i la posició inicial de la seqüència
        delta <- fstp-1
        # Busca la seq del primer RV en la regió 3' (100 últimes posicions) de les seqüències on
        # ja s'ha retallat el primer FW. prmm defineix el màxim de mismatches permesos
        # Guarda totes les coincidències amb les posicions inicial, final, i longitud
        p3.matches <- vmatchPattern(pattern=pr.p3,
                                    subject=subseq(seqs.up,start=fstp,end=endp),
                                    max.mismatch=prmm,fixed=FALSE)
        # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer RV
        flags <- elementNROWS(p3.matches)>=1

        ### Els reads que passen aquí tenen l'amplicó sencer, ja que han coincidit els dos primers
        if(sum(flags))
        { ### Retalla el primer RV per deixar l'amplicó net
          seqs.up <- seqs.up[flags]
          # Torna a fer el procés d'abans i guarda tots els valors de posició inicial del primer RV
          # sobre les seqüències
          pos <- sapply(startIndex(p3.matches)[flags],function(x) x[[1]])+
            delta[flags]
          # Retalla les seqs des de la posició 1 fins on es detecta el primer RV
          seqs.up <- subseq(seqs.up,start=1,end=pos-1)

          ### Colapsa les seqüències 'up' en haplotips + freqüències
          # aplicant la funció del paquet QSutils
          col <- Collapse(seqs.up)
          # Guarda les seqs dels haplotips
          bseqs <- col$hseqs
          # Guarda les freqüències dels haplotips
          nr <- col$nr

          # Aplica la funció del paquet QSutils per ordenar els haplotips en funció del nº de mutacions amb la màster
          # sobre les seqüències i les seves freqüències
          # Cal aplicar DNAStringSet per convertir les seqüències, ja que la funció no les accepta si no són d'aquest tipus
          lst <- SortByMutations(bseqs,nr)
          # Guarda les freqüències dels haplotips obtinguts
          nr <- lst$nr

          # Calcula la freqüència relativa de cada amplicó pertanyent al MID avaluat
          frq_lst <- round(nr/sum(nr)*100,2)
          # Donat que la funció del paquet QSutils retorna els noms dels haplotips diferent al desitjat, es canvien els
          # noms substituint els caràcters '_' per '.', i generant la capçalera del fitxer fasta amb nom de l'haplotip,
          # nombre de reads i freq relativa
          names(lst$bseqs) <- paste(str_replace_all(names(lst$bseqs),'_','.'),nr,frq_lst,sep="|")
          # Generació del fitxer fasta amb tots els haplotips assignats a la carpeta trim
          writeXStringSet(lst$bseqs,file.path(trimDir,up.flnm))

          cat("\nForward seqs, table of read lengths (over 10 rd)\n")
          # 'tapply()' en aquest cas calcula el sumatori (sum) de les freqüències (nº reads) segons les diferents
          # longituds de seqüència
          tbl.len <- tapply(nr,nchar(bseqs),sum)
          # Filtra les longituds amb més de 10 reads (per això el títol amb la funció 'cat()')
          print(tbl.len[tbl.len>=10])
          # Gràfic representant les diferents longituds en X i les freqüències (nº reads) en Y
          plot(as.integer(names(tbl.len)),tbl.len,type="h", # "h" per representar histogrames en forma de línies verticals
               xlab="Read length",ylab="# reads")
          title(main=paste(tt," Str FW"))
        }
      }
      # Actualitza k per indicar l'index a la taula de reports: per cada mostra (2 per pacient) hi haurà 2 resultats
      # corresponents a les cadenes forward i reverse
      k <- k+1

    # Genera un data frame amb el doble de files que el total de mostres i 9 columnes amb els noms dels arguments
    # 'character()' aplicat sobre un nombre enter retorna tants caràcters buits com el nº indiqui
    # 'integer()' sobre un nombre enter retorna tants 0s com el nº indiqui
    FlTbl <- data.frame(File.Name=up.flnm,Pat.ID=samples$Patient.ID[j],
                        Ampl.Nm=primers$Ampl.Nm[ipr],Pr.ID=ipr,
                        Str="fw",Pos=primers$FW.tpos[ipr],Len=mean(width(seqs.up)),
                        Reads=sum(nr),Hpls=length(nr),
                        stringsAsFactors=FALSE)

      # # Nom del fitxer (pacient, regió, PrFw)
      # FlTbl$File.Name[k] <- up.flnm
      # # Identificació del pacient
      # FlTbl$Pat.ID[k] <- samples$Patient.ID[j]
      # # Coordenades de la egió amplificada (5'X o preS1)
      # FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
      # # Identificador del primer (1 o 2) segons la regió que amplifica
      # FlTbl$Pr.ID[k] <- ipr
      # # Cadena (en aquest cas forward)
      # FlTbl$Str[k] <- "fw"
      # # Posició del genoma de HBV en la que comença l'amplicó (després d'eliminar els primers)
      # FlTbl$Pos[k] <- primers$FW.tpos[ipr]
      # # Mitjana de longitud de les seqüències després de retallar-les
      # # Comprovant els valors de width(seqs.up) es pot veure que no totes les seqüències tenen
      # # la mateixa longitud, per tant es fa la mitjana de totes
      # FlTbl$Len[k] <- mean(width(seqs.up)) ####
      # # Sumatori del total de reads en el MID avaluat (mostra) que s'han assignat
      # FlTbl$Reads[k] <- sum(nr)
      # # Total d'haplotips detectats
      # FlTbl$Hpls[k] <- length(nr)

      # Afegeix el nº de reads (del MID concret) idenficats amb la cadena FW al total del pool al qual correspon
      PoolTbl[pool,2] <- PoolTbl[pool,2] + sum(nr)

    pr.res <- data.frame(Tot.reads=totalseqs,matches=sum(flags),shorts=sum(shorts),final.reads=sum(nr))
      # # A l'altra taula de reports afegeix, per la iteració i cadena avaluada:
      # ### En aquest punt aquesta variable correspon als reads que no s'han assignat a la cadena FW!!
      # pr.res[k,1] <- totalseqs
      # # Reads que s'han pogut associar a la cadena FW
      # pr.res[k,2] <- sum(flags)
      # # La variable shorts és 0 osigui que el sumatori també ho és
      # # Hauria de correspondre als reads que no cobreixen l'amplicó sencer
      # pr.res[k,3] <- sum(shorts) ####
      # # Sumatori del total de reads en el MID avaluat (mostra) que s'han assignat
      # pr.res[k,4] <- sum(nr)

      ### Coincidències cadena reverse
      # Aquestes seqüències corresponen a la cadena reverse, és a dir, presenten el primer RV en la regió 5' i la reversa
      # complementària del primer FW a la regió 3'
      # Cal tenir en compte que la variable seqs ara està actualitzada amb les no assignades a cadena FW!
      # Guarda la seq del primer RV específic de la regió avaluada
      pr.dn <- primers$Primer.RV[ ipr ]
      # Busca la seq del primer RV en la regió 5' (posicions 1-100 definides al fitxer de paràmetres)
      # de les seqüències del fitxer .fna del MID avaluat. prmm defineix el màxim de mismatches permesos
      # Guarda totes les coincidències amb les posicions inicial, final, i longitud
      dn.matches <- vmatchPattern(pattern=pr.dn,
                                  subject=subseq(seqs,start=target.io,end=target.in),
                                  max.mismatch=prmm,fixed=FALSE)
      # Resta 1 a la posició inicial de cerca del primer (per defecte, 1)
      delta <- target.io-1
      # Aplica la funció que està definida al fitxer principal
      # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer RV
      flags <- elementNROWS(dn.matches)>=1
      # Torna a buidar aquestes variables per guardar les noves dades
      nr <- integer()
      shorts <- integer()
      # Genera un fitxer .fna, el nom del qual està format per l'ID del pacient, la regió (amplicó) avaluada,
      # i la consecució PrRV.fna (primer reverse)
      dn.flnm <- paste(samples$Patient.ID[j],primers$Ampl.Nm[ipr],"PrRV.fna",
                       sep=".")

      seqs.dn <- "" # Coincidències cadena RV
      if(sum(flags)) # Sumatori de seqüències amb 1 o més coincidències amb la seq del primer
      { # 'startIndex()' retorna una llista amb les posicions inicials de les coincidències del patró cercat
        # Cal tenir cura, ja que si no hi ha coincidència la llista guarda un valor NULL
        # Aplicar el condicional flags sobre la llista de posicions inicials permet eliminar els valors NULL
        # Guarda en una variable tots els valors de posició inicial del primer RV sobre les seqüències
        pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta

        ### Guarda les seqs amb coincidències del primer RV
        seqs.dn <- seqs[flags]


        ### Actualitza el nº de seqs amb les que no presentaven coincidències
        seqs <- seqs[!flags]

        ### Retalla el primer RV de la regió 5'
        # Suma a la posició inicial del primer RV en la seq la seva longitud
        st <- pos + (primers$RV.pos[ipr]-primers$RV.tpos[ipr])
        # Retalla les seqüències des de la posició on acaba el primer RV fins al final
        seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))
        # Guarda només les seqüències amb longitud major a 105 (la posició final on busquem el primer +5)
        seqs.dn <- seqs.dn[width(seqs.dn)>target.in+5]

        ### Actualitza les seqüències de cadena RV fent la seva reversa complementària
        # Important perquè ara les tenim en sentit 3'-5' per poder associar-les a les seves cadenes up (forward)!
        seqs.dn <- reverseComplement(seqs.dn)

        ### Busca ara el primer FW original a la regió inicial (5'). En lloc de fer la reversa complementària del primer,
        # ho fa de totes les seqüències on s'ha trobat el primer RV. Per tant, cal buscar la seq original del primer FW
        # Resta 5 unitats a la posició inicial de cerca del primer (per defecte 1) i calcula el màxim entre aquest valor i 1
        io <- max(target.io-5,1)
        # Calcula la diferència entre aquest màxim i la posició inicial de la seqüència
        delta <- io-1
        # Busca la seq del primer FW en la regió 5' (100 primeres posicions) de les seqüències on
        # ja s'ha retallat el primer RV. prmm defineix el màxim de mismatches permesos
        # Guarda totes les coincidències amb les posicions inicial, final, i longitud
        dn.matches <- vmatchPattern(pattern=pr.up,
                                    subject=subseq(seqs.dn,start=io,end=target.in+5),
                                    max.mismatch=prmm,fixed=FALSE)

        # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer RV
        flags <- elementNROWS(dn.matches)>=1

        ### Els reads que passen aquí tenen l'amplicó sencer, ja que han coincidit els dos primers
        if(sum(flags))
        { ### Retalla el primer FW (ara a 5') per deixar l'amplicó net
          seqs.dn <- seqs.dn[flags]
          # Torna a fer el procés d'abans i guarda tots els valors de posició inicial del primer FW
          # sobre les seqüències
          pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta
          # Suma a la posició inicial del primer FW en la seq la seva longitud
          st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])
          # Retalla les seqüències des de la posició on acaba el primer FW fins al final
          seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))

          ### Colapsa les seqüències 'down' en haplotips + freqüències
          # aplicant la funció del paquet QSutils
          col <- Collapse(seqs.dn)
          # Guarda les seqs dels haplotips
          bseqs <- col$hseqs
          # Guarda les freqüències dels haplotips
          nr <- col$nr

          # Aplica la funció del paquet QSutils per ordenar els haplotips en funció del nº de mutacions amb la màster
          # sobre les seqüències i les seves freqüències
          # Cal aplicar DNAStringSet per convertir les seqüències, ja que la funció no les accepta si no són d'aquest tipus
          lst <- SortByMutations(bseqs,nr)
          # Guarda les freqüències dels haplotips obtinguts
          nr <- lst$nr

          # Calcula la freqüència relativa de cada amplicó pertanyent al MID avaluat
          frq_lst <- round(nr/sum(nr)*100,2)
          # Donat que la funció del paquet QSutils retorna els noms dels haplotips diferent al desitjat, es canvien els
          # noms substituint els caràcters '_' per '.', i generant la capçalera del fitxer fasta amb nom de l'haplotip,
          # nombre de reads i freq relativa
          names(lst$bseqs) <- paste(str_replace_all(names(lst$bseqs),'_','.'),nr,frq_lst,sep="|")
          # Generació del fitxer fasta amb tots els haplotips assignats a la carpeta trim
          writeXStringSet(lst$bseqs,file.path(trimDir,dn.flnm))


          cat("\nReverse seqs, table of read lengths (over 10 rd)\n")
          # 'tapply()' en aquest cas calcula el sumatori (sum) de les freqüències (nº reads) segons les diferents
          # longituds de seqüència
          tbl.len <- tapply(nr,nchar(bseqs),sum)
          # Filtra les longituds amb més de 10 reads (per això el títol amb la funció 'cat()')
          print(tbl.len[tbl.len>=10])
          # Gràfic histograma en forma de línies verticals representant les diferents longituds en X i les freqüències
          # (nº reads) en Y, per les seqs reverse
          plot(as.integer(names(tbl.len)),tbl.len,type="h",
               xlab="Read length",ylab="# reads")
          title(main=paste(tt," Str RV"))
        }
      }

    FlTbl[nrow(FlTbl) + 1,] <- c(dn.flnm,samples$Patient.ID[j],primers$Ampl.Nm[ipr],
                                 ipr,"rv",primers$FW.tpos[ipr],mean(width(seqs.dn)),
                                 sum(nr),length(nr))
      # # Actualitza de nou el valor k per a la propera iteració (mostra) que tornarà a començar per les cadenes up
      # k <- k+1
      # # Guarda tots els resultats igual que en el cas de les cadenes up (forward) excepte columna $Len
      # FlTbl$File.Name[k] <- dn.flnm
      # FlTbl$Pat.ID[k] <- samples$Patient.ID[j]
      # FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
      # FlTbl$Pr.ID[k] <- ipr
      # FlTbl$Str[k] <- "rv"
      # # Després de fer la reversa complementària de les cadenes dn (reverse), ambdues cadenes
      # # es troben el mateix sentit, per tant la posició inicial és la mateixa
      # FlTbl$Pos[k] <- primers$FW.tpos[ipr]
      #
      # FlTbl$Len[k] <- mean(width(seqs.dn))
      # FlTbl$Reads[k] <- sum(nr)
      # FlTbl$Hpls[k] <- length(nr)

      # Afegeix el nº de reads (del MID concret) idenficats amb la cadena FW al total del pool al qual correspon
      PoolTbl[pool,2] <- PoolTbl[pool,2] + sum(nr)

      # Guarda tots els resultats igual que en el cas de les cadenes up (forward) a l'altra taula (nº de reads per pas)
    pr.res[nrow(pr.res) + 1,] <- c(treads,sum(flags),sum(shorts),sum(nr))
      # pr.res[k,1] <- length(seqs) ####
      # pr.res[k,2] <- sum(flags)
      # pr.res[k,3] <- sum(shorts) ####
      # pr.res[k,4] <- sum(nr)

result <- list(pr.res=pr.res,PoolTbl=PoolTbl,FlTbl=FlTbl)

}
