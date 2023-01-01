#########/////////Ejercicio 08/////////#########
install.packages("combinat")
library(combinat)
#Ejercico 8. Resuelva el problema del viajero.
#Obtenga de Google Maps las coordenadas (latitud y longitud) de las siguientes ciudades, Pátzcuaro, 
#CDMX, Cuernavaca, Acapulco, Oaxaca, Guanajuato, Celaya, Puerto Escondido, Puerto Vallarta y 
#Guadalajara investigue como obtener la distancia en “línea recta” de cada una de estas ciudades a 
#cualquier otra mediante una transformación adecuada (verifique sus cálculos en Maps) y resuelva el 
#problema del viajero.

#Calculo de la matriz de distancias.

inRadian <- function(k){  #función para convertir grados a radianes
  k <- k*3.141593/180
  return (k)
}
distancia <- function(lat1,long1,lat2,long2){ #función para calcular la distancia entre dos puntos en la tierra.
  R = 6371                                    #radio terrestre medio en kilómetros
  d = 2*R*asin(sqrt(sin((inRadian(lat2)-inRadian(lat1))/2)^2 + cos(inRadian(lat1))*cos(inRadian(lat2))*sin((inRadian(long2)-inRadian(long1))/2)^2))
  return (d)
}
#coordenadas de las 10 ciudades, se guardan en un data frame.
#1 Patzcuaro, 2 CDMX, 3 Cuernavaca, 4 Acapulco, 5 Oaxaca, 6 Guanajuato, 7 Celaya, 8 Puerto Escondido, 9 Puerto Vallarta, 10 Guadalajara
lat <- c( 19.516389 , 19.419444, 18.918611, 16.86287 , 16.898056, 21.018889 , 20.528889, 15.861944, 20.616667 , 20.666111 )  
lon <- c(-101.609722,-99.145556,-99.234167,-99.887009,-96.414167,-101.262778,-100.815  ,-97.071667,-105.233333,-103.351944)
data <- data.frame("Latitud" = lat, "Longitud" = lon)
rownames(data)  <- c("Patzcuaro","CDMX", "Cuernavaca","Acapulco","Oaxaca","Guanajuato","Celaya","Puerto Escondido","Puerto Vallarta","Guadalajara")
matriz <- matrix( rep(0,10), nrow=10, ncol=10) 
for (i in 1:10){                             #Calculo de las distancias entre cada una de las 10 ciudades, los valores se almacenan en una matriz.
  for (j in 1:10){
    matriz[i,j] = distancia(data[i,1],data[i,2],data[j,1],data[j,2]) 
  }
}

#Obtener los tours.

n   <- 10               #Número de ciudades
Sn  <- permn(seq(1,n))  #Grupo de permutaciones en n elementos. Su cardinalidad es: n!
val <- rep(0,fact(n))   #Para cada permutación este vector guarda la cardinalidad de cada órbita.
length(Sn)              #La longitud de Sn es n!

#Sn[[2]]                 #Esta es la lista con la permutación 2 por ejemplo: 1 2 3 4 5 6 7 8 10 9
#Sn[[2]][i]              #Podemos acceder al décimo elemento de la permutación.

for (i in 1:fact(n)){
  x      <- Sn[i]      #A x se le asigna la i-ésima permutación. 
  tmp    <- rep(0,n)   #Vector temporal donde se guarda la órbita de "1" #Es decir, tmp guardará el "recorrido" definido por la permutación i
  tmp[1] <- x[[1]][1]  #La primer entrada de tmp será igual a la primer entrada de la permutación i.
  for(j in 2:n){tmp[j] <-x[[1]][tmp[j-1]]} #Órbita de la i-esima permutación.
  val[i] <- length(unique(tmp))                      #Las entradas que valen n corresponden a los tours.
}
#table(val)               #Validación
#Tabla de cardinalidades de órbitas de "1". Aquellas con cardinalidad 10 corresponden a los tours.

#Calcular la distancia recorrida en cada tour y tomar el tour con la mínima distancia.

data <- as.data.frame(cbind(Sn,val)) #Juntar en un data frame las permutaciones y el vector de validación 
data1 <- subset(data,val==10)        #Extraer los renglones que corresponden a los tours.
#table(data1$val)
#data1$Sn[[1]][10]

tours <- data1$hSn                    #El vector tours contiene todos los posibles tours 
#str(tours)
#tours[[2]][1]
distancia_tours <- rep(0,fact(n-1))  #El vector distancia_tours guardará la distancia recorrida por cada tour.
length(distancia_tours)
#Calculo de las distancias recorridas
for (i in 1:fact(n-1)){ #i es el número de tours. Esta función calcula la distancia de cada tour
  d <- 0                #d es una variable auxiliar para almacenar las distancias
  for (j in 1:n){       #j es el número de ciudades.
    d <- d + matriz[j,tours[[i]][j]] 
  }
  distancia_tours[[i]] <- d    #Almacena la distancia del i-esimo tour.
}
min(distancia_tours)         #La mínima distancia es 2331.350 km
m <- which(abs(distancia_tours - min(2331.350))<=0.1) #Buscamos los tours que corresponde a esa distancia

tours[[ m[[1]] ]]            #El primer tour esta definido por  1 2 3 4 5 6 7 8 9  10 
d <- 0                       #                                  9 3 5 1 8 7 2 4 10 6
#El tour sería Patzcuaro-PuertoVallarta-Guadalajara-Guanajuato-Celaya-CDMX-Cuernavaca-Oaxaca-PuertoEscondido-Acapulco-Patzcuaro
for (j in 1:n){              
  d <- d + matriz[j,tours[[ m[[1]] ]][j]] 
}
d                            #La distancia recorrida por este tour es justo d=2331.350 km

tours[[ m[[2]] ]]            #El segundo tour esta definido por 1 2 3 4 5 6  7 8 9 10
#                                  4 7 2 8 3 10 6 5 1 9
d <- 0                       #El tour sería Patzcuaro-Acapulco-PuertoEscondido-Oaxaca-Cuernavaca-CDMX-Celaya-Guanajuato-Guadalajara-PuertoVAllarta-Patzcuaro
for (j in 1:n){              
  d <- d + matriz[j,tours[[ m[[2]] ]][j]] 
}
d                            #La distancia recorrida por este tour es justo d=2331.350 km
