#1
library(viridis)
library(Biostrings)
library(DECIPHER)
library(ade4)
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(viridis)
library(ggplot2)

#2
#Toma el ID del virus y a partir de read.GenBank encuentra su genoma en el sitio web NCBI virus
virus<-c("NC_045512","OX646281.1","LR991699","MT324062","PP519172","MT126808","OQ430688","MT012098","OM918219")
virus_read<-read.GenBank(virus)

#3
str(virus_read)

#4
#Escribe un archivo fasta con todos los genomas juntos
write.dna(virus_read,file="Coronavirus_seq.fasta",format="fasta")

#5
#Lee el archivo y guarda la información en virus_seq
virus_seq<-readDNAStringSet("Coronavirus_seq.fasta",format="fasta")
class(virus_seq)
virus_seq

#6
#Orienta los nucleotidos en el mismo sentido
virus_seq_notalign<-OrientNucleotides(virus_seq)

#7
#Alinea los nucleotidos
virus_seq_align<-AlignSeqs(virus_seq_notalign)
#Visualiza resultado de alineamiento en navegador
BrowseSeqs(virus_seq_align)

#8
#Genera archivo de nucleotidos alineados
writeXStringSet(virus_seq_align,file="coronavirus_seq_align.fasta")

#9
#Leer archivo
virus_aligned<-read.alignment("coronavirus_seq_align.fasta",format="fasta")

#10
#Matriz de distancia entre sec alineadas
matriz_dist<-dist.alignment(virus_aligned,matrix="similarity")
temp<-as.data.frame(as.matrix(matriz_dist))
table.paint(temp,cleg=0,clabel.col=.5)+scale_color_viridis()

#11
virus_tree<-nj(matriz_dist)
class(virus_tree)

#12
virus_tree<-ladderize(virus_tree)
plot(virus_tree)