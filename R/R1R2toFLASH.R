
min.len <- 200   # Longitud mínima per considerar una seqüència
min.ov <- 20     # Mínim solapament (en nt) entre R1 and R2
max.ov <- 300    # Màxim solapament (en nt) entre R1 and R2
err.lv <- 0.10   # Fracció de diferències acceptades en el solapament
flash <- "C:/FLASH/flash.exe" # Carpeta on es troba l'executable FLASH.

# Indicador de les variables min, max i error level
flash.opts <- paste("-m",min.ov,"-M",max.ov,"-x",err.lv)

### Chunck size to be used by FastqStreamer() --> Buscar info
# Nombre de registres successius a retornar a cada rendiment (yield)
chunck.sz <- 1.e6

runfiles <- list.files(runDir)

#--------- Documentació de la funció----------#

#' @param runfiles Fastq files from Illumina (run directory folder)
#' @param min.len Minimum length to consider a sequence
#' @param min.ov Minimum overlap (in nt) between R1 and R2
#' @param max.ov Maximum overlap (in nt) between R1 and R2
#' @param err.lv Mismatch fraction accepted in overlapping
#' @param flash Folder path containing FLASH executable
#' @param chunck.sz Chunck size to be used by FastqStreamer() function

R1R2toFLASH <- function(runfiles,min.len=200,min.ov=20,max.ov=300,err.lv=0.10,flash,chunck.sz=1.e6)
{

# La funció sub() permet substituir un patró pel que indiquem com 2n argument
# En aquest cas, les variables runfiles i snms son idèntiques (de moment)
snms <- sub("\\.fastq$","",runfiles)

# Taula on es disposa els nom de cada fitxer fastq en la primera columna i després es
# separa cadascun del termes del nom del fitxer en successives columnes
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
# Elimina les columnes 3 i 5 de la taula (que indiquen L001 i l'extensió)
parts <- parts[,-c(3,5)]
# Assigna els noms a les columnes: ID pacient, ID mostra/pool, read
colnames(parts) <- c("PatID","SmplID","Read")

# Guarda a la variable la ruta dels fitxers que es troben a la carpeta run
flnms <- file.path(runDir,runfiles)

# Dels 4 fitxers que hi havia, guarda en 2 variables els que corresponen a R1 i R2
R1.flnms <- flnms[parts[,3]=="R1"] # parts[,3] = columna 3 (read) de la taula parts
R2.flnms <- flnms[parts[,3]=="R2"]

# Guarda dos noms de fitxer corresponents a R1: pool, ID de mostra i extensió flash.fastq
# Com hi ha 2 fitxers amb R1, vol dir que avaluem 2 regions o pools
out.flnms <- paste(parts[parts[,3]=="R1",1],parts[parts[,3]=="R1",2],
                   "flash.fastq",sep="_")

# Genera la ruta dels fitxers definits abans a la carpeta flash
out.flnms <- file.path(flashDir,out.flnms)

# Guarda de la taula només els R1
parts <- parts[parts[,3]=="R1",,drop=FALSE]
# Construeix una matriu 2x2 amb tot 0
res <- matrix(0,nrow=length(out.flnms),ncol=2) # length(out.flnms)= 2 en aquest cas, que son els pools
# Assigna com a nom de fila el pool (regió de VHB) i com a columnes els reads segons si s'ha donat o no extensió
rownames(res) <- paste(parts[,1],parts[,2],sep="_")
colnames(res) <- c("Extended","NoExtd")

# Defineix les opcions necessàries per a l'execució de FLASH
flash.opts <- paste("-m",min.ov,"-M",max.ov,"-x",err.lv)

#----- Caldria revisar aquesta part --------#
## Aplica la funció 'executeFLASH()' del paquet, que permet realitzar l'extensió
# dels reads R1 i R2 i guardar el nº de reads units (extended) i no units (no extended)
# Parteix dels fitxers R1 i R2 de cadascun dels pools, així com el nom del fitxer fastq
# resultant que es guardarà a la carpeta flash
p1 <- executeFLASH(R1.flnms[1],R2.flnms[1],out.flnms[1])
p2 <- executeFLASH(R1.flnms[2],R2.flnms[2],out.flnms{2})

# Guarda a la taula de resultats (res) el nº de reads units i no units per FLASH
# per a cadascun dels pools avaluats
res[1,] <- p1
res[2,] <- p2
#------------------------------------#

# Guarda un dataframe amb les dades resultants del FLASH (taula res), afegint la
# columna Yield calculada dividint els reads extended entre el total *100
df.res <- data.frame(res,Yield=round(res[,1]/(res[,1]+res[,2])*100,1))
# Guarda el fitxer de report de FLASH en format .txt
txt.flnm <- file.path(repDir,"FLASH_report.txt")
# Comandes per omplir el fitxer txt generat. Recull els paràmetres indicats al FLASH
# i la taula df.res que conté els resultats obtinguts
sink(txt.flnm)
cat("\nExtending Illumina reads by FLASH\n")
cat("\nFLASH parameters:")
cat("\n    Minimum overlap:",min.ov)
cat("\n    Maximum overlap:",max.ov)
cat("\n        Error level:",err.lv,"\n\n")
print(df.res)
sink() # Tanca el fitxer

flash.res <- df.res # Necessari??
# Guarda a la carpeta reports la taula en format RData
save(flash.res,file=file.path(repDir,"FLASH_table.RData"))

# Genera el pdf que contindrà el gràfic barplot dels resultats FLASH
pdf.flnm <- file.path(repDir,"FLASH_barplot.pdf")
pdf(pdf.flnm,paper="a4",width=6,height=10)
par(mfrow=c(2,1))

library(RColorBrewer)

pal=brewer.pal(8,"Dark2") # Crea la paleta de colors
M <- data.matrix(df.res[,1:2]) # Guarda les columnes 1 i 2 de la taula de resultats
ymx <- max(M)*1.15 # Defineix límit superior del gràfic
barplot(t(M),beside=TRUE,las=2,col=pal[1:2],ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal[1:2],cex=0.8,
       legend=c("Extended","Not extended"))
title(main="FLASH results on paired-ends",line=2.5)
title(main="Yield in number of reads",cex.main=1,line=1)

# Genera un altre barplot on s'inclouen els resultats del Yield de FLASH
bp <- barplot(df.res$Yield,col="Lavender",ylim=c(0,100),ylab="Yield",
              names.arg=rownames(df.res),las=2)
text(bp,10,paste(df.res$Yield,"%",sep=""),srt=90,col="navy")
title(main="Yield in percentage",cex.main=1,line=1)

dev.off()

}
