library(seqinr)
library(RColorBrewer)
options(max.print = 1000000)

#Leer files
wuhan<-read.fasta(file="wuhan.fasta")
uk<-read.fasta(file="uk.fasta")
alfa<-read.fasta(file="alfa.fasta")
southa<-read.fasta(file="sudafrica.fasta")
beta<-read.fasta(file="beta.fasta")
brazil<-read.fasta(file="brazil.fasta")
gamma<-read.fasta(file="gamma.fasta")
india<-read.fasta(file="india.fasta")
delta<-read.fasta(file="delta.fasta")

#Obtener longitud de genoma
wuhanlen<-getLength(wuhan)
uklen<-getLength(uk)
alfalen<-getLength(alfa)
southalen<-getLength(southa)
betalen<-getLength(beta)
brazillen<-getLength(brazil)
gammalen<-getLength(gamma)
indialen<-getLength(india)
deltalen<-getLength(delta)

#Obtener dos vectores, uno con longitud y otro con nombres de virus
len<-c(wuhanlen,uklen,alfalen,southalen,betalen,brazillen,gammalen,indialen,deltalen)
names<-c("Wuhan","Reino Unido","Alfa","Sudafrica", "Beta", "Brazil","Gamma","India","Delta")

#utilizar función para generar grafico de barras
color1<-brewer.pal(9,"Spectral")

bp <- barplot(len, beside=TRUE, col=color1, ylab="Longitud", ylim=c(29800, 29940), main= "Longitud de los genomas")

legend("topright", legend=names, fill=color1, cex=0.65)

text(x = bp, y = len, label = len, pos = 1, cex = 0.8)

#Obtener vector ADN de cada genoma
function_vector<-function(nchar,genoma){
  genomavector<-c("")
  for (n in 1:nchar){
    genomavector<-append(genomavector,genoma[[1]][n])
  }
  genoma_total<-(paste(genomavector,collapse=""))
  genoma_total<-(toupper(genoma_total))
  return(genoma_total)
}
wuhan_gen<-function_vector(wuhanlen,wuhan)
uk_gen<-function_vector(uklen,uk)
alfa_gen<-function_vector(alfalen,alfa)
southa_gen<-function_vector(southalen,southa)
beta_gen<-function_vector(betalen,beta)
brazil_gen<-function_vector(brazillen,brazil)
gamma_gen<-function_vector(gammalen,gamma)
india_gen<-function_vector(indialen,india)
delta_gen<-function_vector(deltalen,delta)

#Composición de nucleotidos
function_porcentaje<-function(vector_adn, nchar){
  
  conteo_a <- 0
  conteo_c <- 0
  conteo_g <- 0
  conteo_t <- 0
  
  for (i in 1:nchar) {
    nuc <- substr(vector_adn, i, i)
    if (nuc == "A") {
      conteo_a <- conteo_a + 1
    } else if (nuc == "C") {
      conteo_c <- conteo_c + 1
    } else if (nuc == "G") {
      conteo_g <- conteo_g + 1
    } else if (nuc == "T") {
      conteo_t <- conteo_t + 1
    }
  }
  
  nucleotides_num<-c(conteo_a,conteo_c,conteo_g,conteo_t)
  cat("Existen", conteo_a, "A en el genoma\n")
  cat("Existen", conteo_c, "C en el genoma\n")
  cat("Existen", conteo_g, "G en el genoma\n")
  cat("Existen", conteo_t, "T en el genoma\n")
  
  porcentaje_A <- (conteo_a/nchar)*100
  porcentaje_C <- (conteo_c/nchar)*100
  porcentaje_G <- (conteo_g/nchar)*100
  porcentaje_T <- (conteo_t/nchar)*100
  nucleotides_porc<-c(porcentaje_A,porcentaje_C,porcentaje_G,porcentaje_T)
  
  cat("Porcentaje de A:", porcentaje_A, "%\n")
  cat("Porcentaje de C:", porcentaje_C, "%\n")
  cat("Porcentaje de G:", porcentaje_G, "%\n")
  cat("Porcentaje de T:", porcentaje_T, "%\n")
  return(nucleotides_num)
}

cat("Composición del genoma SARS-CoV-2 proveniente de Wuhan:\n")
wuhan_comp<-function_porcentaje(wuhan_gen, wuhanlen)

cat("Composición del genoma SARS-CoV-2 proveniente del Reino Unido:\n")
uk_comp<-function_porcentaje(uk_gen, uklen)

cat("Composición del genoma SARS-CoV-2 variante Alfa:\n")
alfa_comp<-function_porcentaje(alfa_gen, alfalen)

cat("Composición del genoma SARS-CoV-2 proveniente de Sudafrica:\n")
southa_comp<-function_porcentaje(southa_gen, southalen)

cat("Composición del genoma SARS-CoV-2 variante Beta:\n")
beta_comp<-function_porcentaje(beta_gen, betalen)

cat("Composición del genoma SARS-CoV-2 proveniente de Brazil:\n")
brazil_comp<-function_porcentaje(brazil_gen, brazillen)

cat("Composición del genoma SARS-CoV-2 variante Gamma:\n")
gamma_comp<-function_porcentaje(gamma_gen, gammalen)

cat("Composición del genoma SARS-CoV-2 proveniente de India:\n")
india_comp<-function_porcentaje(india_gen, indialen)

cat("Composición del genoma SARS-CoV-2 variante Delta:\n")
delta_comp<-function_porcentaje(delta_gen, deltalen)

#cree un vector en donde guarde los nombres de los nucleotidos
Nucleotides<-c("Adenina","Citosina","Guanina","Timina")

#cree un dataframe con los compuestos de nucleotidos de cada virus
df<-data.frame(Wuhan=wuhan_comp,UK=uk_comp,A=alfa_comp,Sudafrica=southa_comp,B=beta_comp,Brazil=brazil_comp,G=gamma_comp,India=india_comp,D=delta_comp)
head(df)
#converti el dataframe en una matriz para que la pudiera utilizar con barplot
df_matrix<-data.matrix(df)

color<-brewer.pal(4,"PuRd")
#guarde mi barplot en una variable
comp_n<-barplot(df_matrix,col=color,main="Composición de Nucleotidos",ylim=c(0,39000))
#la leyenda del barplot
legend("topright",legend=Nucleotides,fill=color, cex=0.5)

#Esta función toma como argumento un vector con la composición de nucleotidos y el numero de orden del virus en la grafica
functiontext<-function(vector_nucleotidos,num){
  #guarda el primer dato del vector en un vector nuevo
  vector_nuevo<-c(vector_nucleotidos[1])
  #este loop le suma el segundo dato (cant de nucleotido) al anterior para poder tener la ubicación en y para posicionar el texto
  #Lo realiza de 2 a 4 porque ya tenemos nuestra primera posición en "y" y queremos tener las posiciones de los 3 compuestos restantes
  for (i in (2:4)){
    value<-vector_nucleotidos[i]+vector_nuevo[length(vector_nuevo)]
    #toma el valor del nucleotido y se lo suma al primer valor que esta actualmente en el vector nuevo
    vector_nuevo<-append(vector_nuevo,value)
  }
  print(vector_nuevo)
  #num viene siendo el numero de orden en el que esta en la grafica, comp_n son las posiciones en x, vector_nuevo ahora tiene las posiciones en y, label toma el valor de cada dato
  for (i in (1:4)){
    text(x = comp_n[num], y = vector_nuevo[i], label = vector_nucleotidos[i], pos = 1,cex = 0.5)
  }
}


