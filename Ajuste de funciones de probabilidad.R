#Debe de hacer esto para cada vector de montos que le pasemos
#Familia            parametro1         parametro2   parmetro3 EstadisticoKS P-valueKS EstadisticoAD P-ValueAD  
#Exponencial        rate               ----         ----      ----          ----      ----          ----                       
#Gamma              shape              rate         ----      ----          ----      ----          ----                        
#LogNormal          meanlog            sdlog        ----      ----          ----      ----          ----                       
#Weibull            shape              scale        ----      ----          ----      ----          ----                       
#Pareto             shape              scale        ----      ----          ----      ----          ----                       
#Pareto Tipo 1      shape              min          ----      ----          ----      ----          ----                       
#Generalized Pareto shape1             shape2       rate      ----          ----      ----          ----                       
#Burr               theta              alfa         beta      ----          ----      ----          ----                       

#Las 3 parametrizaciones de la Pareto que use
#rpareto(n, shape, scale) #de actuar
#rpareto1(n, shape, min)  #de actuar
#rgenpareto(n, shape1, shape2, scale) #de actuar

library(actuar)
library(MASS)         #De aquí usamos fitdistr
library(survival)
library(fitdistrplus) #De aquí usamos fitdist

#Función que estima parametros y hace las pruebas de bondad de ajuste.
#Recibe un vector como parametro
#Devuelve una matriz con los parametros estimados, estadisticas  y p-value

myfunction <- function(x){
  
  #Quitamos los NA'S
  x <- x[!is.na(x)]
  
  #EXPONENCIAL
  ##################################################################################################
  #set.seed(100)
  #x <- rexp(500,rate = 5) #Vector de datos
  #Estimador por momentos
  #mme <- fitdist(x,"exp",method ="mme") 
  #mme
  #Estimador por maxima verosimilitud                                     
  #mle <- fitdist(x,"exp",start=c(rate=mme$estimate[[1]]),method ="mle") 
  #mle
  #Estimador por cuantiles
  qme <- fitdist(x,"exp", probs=c(1/3),method ="qme") 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"exp", gof = "AD",start=c(rate=qme$estimate[[1]]),method ="mge") 
  mge
  estim <- mge

  #Pruebas de bondad de ajuste
  #ks.test(x, "pexp",rate=estim$estimate[[1]]) 
  KS <- ks.test(x, "pexp",rate=estim$estimate[[1]])$statistic[[1]]
  #KS   #Estadistico Kolomogorov Smirnov
  KSPV <- ks.test(x, "pexp",rate=estim$estimate[[1]])$p.value
  #KSPV #P-value Kolomgorov Smirnov
  #ad.test(x, pexp, rate=estim$estimate[[1]])
  AD <- ad.test(x, pexp, rate=estim$estimate[[1]])$statistic[[1]]
  #AD   #Estadistico Anderson Darling
  ADPV <- ad.test(x, pexp, rate=estim$estimate[[1]])$p.value
  #ADPV #P-value Anderson Darling

  tabla <- matrix("",8,7) #Matriz con datos
  tabla
  colnames(tabla) <- c("Parametro 1", "Parametro 2", "Parametro 3", "Estadistico KS", "P-value KS", "Estadistico AD", "P-value AD")
  rownames(tabla) <- c("Exponencial", "Gamma", "LogNormal", "Weibull", "Pareto tipo 1", "Pareto tipo 2", "Pareto tipo 3", "Burr")
  tabla[1,1]<-estim$estimate[[1]];tabla[1,2]<-"--------";tabla[1,3]<-"--------";tabla[1,4]<-KS;tabla[1,5]<-KSPV;tabla[1,6]<-AD;tabla[1,7]<-ADPV;tabla
  rm(qme,mge,estim)
  
  #GAMMA
  ##################################################################################################
  #set.seed(100)
  #x <- rgamma(200, shape = 1, rate = 1/2)
  #Estimador por momentos
  #mme <- fitdist(x, "gamma",method ="mme") 
  #mme
  #Estimdor por maxima verosimilitud
  #mle <- fitdist(x,"gamma",start=c(shape = mme$estimate[[1]],rate=mme$estimate[[2]]),method ="mle") 
  #mle
  #Estimador por cuantiles
  qme <- fitdist(x,"gamma", probs=c(1/3,2/3),method ="qme") 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"gamma", gof = "AD",start=c(shape=qme$estimate[[1]],rate=qme$estimate[[2]]),method ="mge") 
  mge
  estim <- mge

  #Pruebas de bondad de ajuste
  #ks.test(x, "pgamma",shape=estim$estimate[[1]],rate=estim$estimate[[2]]) 
  KS <- ks.test(x, "pgamma",shape=estim$estimate[[1]],rate=estim$estimate[[2]])$statistic[[1]]
  #KS #Estadistico KS
  KSPV <- ks.test(x, "pgamma",shape=estim$estimate[[1]],rate=estim$estimate[[2]])$p.value
  #KSPV #P-value KS
  #ad.test(x, pgamma, shape=estim$estimate[[1]], rate=estim$estimate[[2]])
  AD <- ad.test(x, pgamma, shape=estim$estimate[[1]], rate=estim$estimate[[2]])$statistic[[1]]
  #AD #Estadistico AD
  ADPV <- ad.test(x, pgamma, shape=estim$estimate[[1]], rate=estim$estimate[[2]])$p.value
  #ADPV #P-value AD

  tabla[2,1]<-estim$estimate[[1]];tabla[2,2]<-estim$estimate[[2]];tabla[2,3]<-"--------";tabla[2,4]<-KS;tabla[2,5]<-KSPV;tabla[2,6]<-AD;tabla[2,7]<-ADPV;tabla
  tabla
  rm(qme,mge,estim)
  
  
  #LOGNORMAL
  ##################################################################################################
  #set.seed(100)
  #x <- rlnorm(100, meanlog = 2, sdlog = 4)
  #Estimador por momentos
  #mme <- fitdist(x, "lnorm",method ="mme") #Estimadores
  #mme
  #Estimador por maxima verosimilitud
  #mle <- fitdist(x,"lnorm",start=c(meanlog = mme$estimate[[1]],sdlog=mme$estimate[[2]]),method ="mle") 
  #mle
  #Estimador por cuantiles
  qme <- fitdist(x,"lnorm", probs=c(1/3,2/3),method ="qme") 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"lnorm", gof = "AD",start=c(meanlog=qme$estimate[[1]],sdlog=qme$estimate[[2]]),method ="mge") 
  mge
  estim <- mge

  #ks.test(x, "plnorm",meanlog=estim$estimate[[1]],sdlog=estim$estimate[[2]]) #Pruebas de bondad de ajuste
  KS <- ks.test(x, "plnorm",meanlog=estim$estimate[[1]],sdlog=estim$estimate[[2]])$statistic[[1]]
  #KS #Estadistico KS
  KSPV <- ks.test(x, "plnorm",meanlog=estim$estimate[[1]],sdlog=estim$estimate[[2]])$p.value
  #KSPV #P-value KS
  #ad.test(x, plnorm, meanlog=estim$estimate[[1]],sdlog=estim$estimate[[2]])
  AD <- ad.test(x, plnorm, meanlog=estim$estimate[[1]],sdlog=estim$estimate[[2]])$statistic[[1]]
  #AD #Estadistico AD
  ADPV <- ad.test(x, plnorm, meanlog=estim$estimate[[1]],sdlog=estim$estimate[[2]])$p.value
  #ADPV #P-value AD

  tabla[3,1]<-estim$estimate[[1]];tabla[3,2]<-estim$estimate[[2]];tabla[3,3]<-"--------";tabla[3,4]<-KS;tabla[3,5]<-KSPV;tabla[3,6]<-AD;tabla[3,7]<-ADPV;tabla
  rm(qme,mge,estim)
  
  #WEIBULL
  #################################################################################################
  #set.seed(100)
  #x <- rweibull(200,shape=1, scale = 1/2)
  #Estimador por momentos
  #memp  <-  function(x, order) mean(x^order) #Defino el momento empirico
  #mme <- fitdist(x, "weibull", order=c(1, 2), memp=memp,method ="mme") #Estimadores
  #mme
  #Estimador por maxima verosimilitud
  #mle <- fitdist(x,"weibull",start=c(shape = mme$estimate[[1]],scale=mme$estimate[[2]]),method ="mle") 
  #mle
  #Estimador por cuantiles
  qme <- fitdist(x,"weibull", probs=c(1/3,2/3),method ="qme") 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"weibull", gof = "AD",start=c(shape=qme$estimate[[1]],scale=qme$estimate[[2]]),method ="mge") 
  mge
  estim <- mge 
  
  #Pruebas de bondad de ajuste
  #ks.test(x, "pweibull",shape=estim$estimate[[1]],scale=estim$estimate[[2]]) 
  KS <- ks.test(x, "pweibull",shape=estim$estimate[[1]],scale=estim$estimate[[2]])$statistic[[1]]
  #KS #Estadistico KS
  KSPV <- ks.test(x, "pweibull",shape=estim$estimate[[1]],scale=estim$estimate[[2]])$p.value
  #KSPV #P-value KS
  #ad.test(x, pweibull,shape=estim$estimate[[1]],scale=estim$estimate[[2]])
  AD <- ad.test(x, pweibull,shape=estim$estimate[[1]],scale=estim$estimate[[2]])$statistic[[1]]
  #AD #Estadistico AD
  ADPV <- ad.test(x, pweibull,shape=estim$estimate[[1]],scale=estim$estimate[[2]])$p.value
  #ADPV #P-value AD 

  tabla[4,1]<-estim$estimate[[1]];tabla[4,2]<-estim$estimate[[2]];tabla[4,3]<-"--------";tabla[4,4]<-KS;tabla[4,5]<-KSPV;tabla[4,6]<-AD;tabla[4,7]<-ADPV;tabla
  rm(qme,mge,estim)
  
  #PARETO
  ###################################################################################################
  #set.seed(100)
  #x <- rpareto(500, shape=50, scale=50)
  #x
  #Estimador por maxima verosimilitud
  #mle <- fitdist(x,"pareto",method ="mle") 
  #mle
  #Estimador por cuantiles
  qme <- fitdist(x,"pareto", probs=c(1/3,2/3),start=c(shape = 1,scale=1),method ="qme") 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"pareto", gof = "AD",start=c(shape=qme$estimate[[1]],scale=qme$estimate[[2]]),method ="mge") 
  mge

  if( (exists("qme") && !is.na(qme$estimate[[1]]) && !is.na(qme$estimate[[2]])) || (exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]]) )){
    if( exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]])){
      estim <- mge
    }else{
      estim <- qme
    }
    #Pruebas de bondad de ajuste
    #ks.test(x, "ppareto",shape=estim$estimate[[1]],scale=estim$estimate[[2]]) 
    KS <- ks.test(x, "ppareto",shape=estim$estimate[[1]],scale=estim$estimate[[2]])$statistic[[1]]
    #KS #Estadistico KS
    KSPV <- ks.test(x, "ppareto",shape=estim$estimate[[1]],scale=estim$estimate[[2]])$p.value
    #KSPV #P-value KS
    #ad.test(x, ppareto,shape=estim$estimate[[1]],scale=estim$estimate[[2]])
    AD <- ad.test(x, ppareto,shape=estim$estimate[[1]],scale=estim$estimate[[2]])$statistic[[1]]
    #AD #Estadistico AD
    ADPV <- ad.test(x, ppareto,shape=estim$estimate[[1]],scale=estim$estimate[[2]])$p.value
    #ADPV #P-value AD
    
    tabla[5,1] <- estim$estimate[[1]]; tabla[5,2] <- estim$estimate[[2]];tabla[5,3]<-"--------"; tabla[5,4] <- KS; tabla[5,5] <- KSPV; tabla[5,6] <- AD; tabla[5,7] <- ADPV;tabla
  }else { #Si no se pueden calcular los estimadores
    tabla[5,1]<-"--------";tabla[5,2]<-"--------";tabla[5,3]<-"--------";tabla[5,4]<-"--------";tabla[5,5]<-"--------";tabla[5,6]<-"--------";tabla[5,7]<-"--------";tabla
  }
  rm(qme,mge,estim)
  
  #PARETO1
  ###################################################################################################
  #set.seed(200)
  #x <- rpareto1(500, shape=500, min=0.5)
  #Estimador por momentos
  #memp  <-  function(x, order) mean(x^order) #Defino el momento empirico
  #mme <- fitdist(x,"pareto1", order=c(1, 2), memp=memp,start=c(shape = 10,min=min(x)),method ="mme") 
  #mme
  #Estimador por cuantiles
  qme <- fitdist(x,"pareto1", probs=c(1/3,2/3),start=c(shape = 1,min=min(x)),method ="qme") 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"pareto1", gof = "AD",start=c(shape=qme$estimate[[1]],min=qme$estimate[[2]]),method ="mge") 
  mge
  
  
  if( (exists("qme") && !is.na(qme$estimate[[1]]) && !is.na(qme$estimate[[2]])) || (exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]]) )){
    if( exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]])){
      estim <- mge
    }else{
      estim <- qme
    }
    #Pruebas de bondad de ajuste
    #ks.test(x, "ppareto1",shape=estim$estimate[[1]],min=estim$estimate[[2]]) 
    KS <- ks.test(x, "ppareto1",shape=estim$estimate[[1]],min=estim$estimate[[2]])$statistic[[1]]
    #KS #Estadistico KS
    KSPV <- ks.test(x, "ppareto1",shape=estim$estimate[[1]],min=estim$estimate[[2]])$p.value
    #KSPV #P-value KS
    #ad.test(x, ppareto1,shape=estim$estimate[[1]],min=estim$estimate[[2]])
    AD <- ad.test(x, ppareto1,shape=estim$estimate[[1]],min=estim$estimate[[2]])$statistic[[1]]
    #AD #Estadistico AD
    ADPV <- ad.test(x, ppareto1,shape=estim$estimate[[1]],min=estim$estimate[[2]])$p.value
    #ADPV #P-value AD
    
    tabla[6,1]<-estim$estimate[[1]];tabla[6,2]<-estim$estimate[[2]];tabla[6,3]<-"--------";tabla[6,4]<-KS;tabla[6,5]<-KSPV;tabla[6,6]<-AD;tabla[6,7]<-ADPV;tabla
  }else { 
    #Si no se pueden calcular los estimadores
    tabla[6,1]<-"--------";tabla[6,2]<-"--------";tabla[6,3]<-"--------";tabla[6,4]<-"--------";tabla[6,5]<-"--------";tabla[6,6]<-"--------";tabla[6,7]<-"--------";tabla
  }
  rm(qme,mge,estim)
  
  
  

  #GENPARETO
  ###################################################################################################
  #set.seed(200)
  #x <- rgenpareto(500, shape1=2, shape2=4, scale=8)
  #Estimador por momentos
  #memp  <-  function(x, order) mean(x^order) #Definir el momento empirico
  #mme <- fitdist(x,"pareto1", order=c(1, 2), memp=memp,start=c(shape = 10,min=min(x)),method ="mme") 
  #mme
  #Estimador por momentos
  #memp  <-  function(x, order) mean(x^order) #Defino el momento empirico
  #mme <- fitdist(x, "genpareto", order=c(1,2,3), memp=memp,start=c(shape1=qme$estimate[[1]],shape1=qme$estimate[[2]],scale=qme$estimate[[3]]),method ="mme") #Estimadores
  #mme
  #Estimador por cuantiles
  qme <- fitdist(x,"genpareto", probs=c(1/4,2/4,3/4),start=c(shape1=mean(x),shape2=mean(x),scale=min(x)),method ="qme") 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"genpareto", gof = "AD",start=c(shape1=qme$estimate[[1]],shape1=qme$estimate[[2]],scale=qme$estimate[[3]]),method ="mge") 
  mge

  
  if( (exists("qme") && !is.na(qme$estimate[[1]]) && !is.na(qme$estimate[[2]]) && !is.na(qme$estimate[[3]])) || (exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]]) && !is.na(mge$estimate[[3]]))){
    if( exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]]) && !is.na(mge$estimate[[3]]) ){
      estim <- mge
    }else{
      estim <- qme
    }   
    #Pruebas de bondad de ajuste
    #ks.test(x, "pgenpareto",shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]]) 
    KS <- ks.test(x, "pgenpareto",shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$statistic[[1]]
    #KS #Estadistico KS
    KSPV <- ks.test(x, "pgenpareto",shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$p.value
    #KSPV #P-value KS
    #ad.test(x, pgenpareto,shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])
    AD <- ad.test(x, pgenpareto,shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$statistic[[1]]
    #AD #Estadistico AD
    ADPV <- ad.test(x, pgenpareto,shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$p.value
    #ADPV #P-value AD
    
    
    tabla[7,1]<-estim$estimate[[1]];tabla[7,2]<-estim$estimate[[2]];tabla[7,3]<-estim$estimate[[3]];tabla[7,4]<-KS;tabla[7,5]<-KSPV;tabla[7,6]<-AD;tabla[7,7]<-ADPV;tabla
  }else { #Si no se pueden calcular los estimadores
    tabla[7,1]<-"--------";tabla[7,2]<-"--------";tabla[7,3]<-"--------";tabla[7,4]<-"--------";tabla[7,5]<-"--------";tabla[7,6]<-"--------";tabla[7,7]<-"--------";tabla
  }
  rm(qme,mge,estim)
  
  #BURR
  ################################################################################################################################
  #set.seed(100)
  #x <- rburr(500, shape1=1, shape2=1, scale = 1)
  #mean(x)
  #Estimador por maxima verosimilitud
  #mle <- fitdist(x, "burr",method ="mle",start=c(shape1=30,shape2=30,scale=5)) 
  #mle$estimate[[3]]
  #Estimador por cuantiles
  qme <- fitdist(x,"burr", probs=c(1/4,2/4,3/4),method ="qme",start=c(shape1=30,shape2=30,scale=mean(x))) 
  qme
  #Estimador por maxima bondad de ajuste
  mge <- fitdist(x,"burr", gof = "AD",start=c(shape1=qme$estimate[[1]],shape2=qme$estimate[[2]],scale=qme$estimate[[3]]),method ="mge") 
  mge
  
  if( (exists("qme") && !is.na(qme$estimate[[1]]) && !is.na(qme$estimate[[2]]) && !is.na(qme$estimate[[3]])) || (exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]]) && !is.na(mge$estimate[[3]]))){
    if( exists("mge") && !is.na(mge$estimate[[1]]) && !is.na(mge$estimate[[2]]) && !is.na(mge$estimate[[3]]) ){
      estim <- mge
    }else{
      estim <- qme
    }   
    #Pruebas de bondad de ajuste
    #ks.test(x, "pburr",shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]]) 
    KS <- ks.test(x, "pburr",shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$statistic[[1]]
    #KS #Estadistico KS
    KSPV <- ks.test(x, "pburr",shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$p.value
    #KSPV #P-value KS
    #ad.test(x,pburr,shape1=estim$estimate[[1]], shape2=estim$estimate[[2]],scale=estim$estimate[[3]])
    AD <- ad.test(x, pburr,shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$statistic[[1]]
    #AD #Estadistico AD
    ADPV <- ad.test(x, pburr,shape1=estim$estimate[[1]],shape2=estim$estimate[[2]],scale=estim$estimate[[3]])$p.value
    #ADPV #P-value AD
    
    tabla[8,1]<-estim$estimate[[1]];tabla[8,2]<-estim$estimate[[2]];tabla[8,3]<-estim$estimate[[3]];tabla[8,4]<-KS;tabla[8,5]<-KSPV;tabla[8,6]<-AD;tabla[8,7]<-ADPV;tabla
  }else { #Si no se pueden calcular los estimadores
    tabla[8,1]<-"--------";tabla[8,2]<-"--------";tabla[8,3]<-"--------";tabla[8,4]<-"--------";tabla[8,5]<-"--------";tabla[8,6]<-"--------";tabla[8,7]<-"--------";tabla
  }
  rm(qme,mge,estim)
  
  return(tabla)
}






#Función final, hace lo que queriamos.
#Función que recibe n vectores.
#Imprime n matrices con los estimados, estadisticas y p-values de cada uno de los n vectores.
ajuste <- function(...){   #Recibe n vectores
  u <- list(...) #Guarda los n vectores en una lista, esta lista se le asigna a la variable u. 
  for(i in 1:length(u)) {  #Para cada vector en la lista u, imprime la matriz de resultados.
    v <- myfunction(u[[i]])
    writeLines("\n\n")
    print(paste0("Ajuste del vector: ", i))
    print(v)
    writeLines("\n\n")
  }
}

#Prueba para 2 vectores, uno es de datos exp, el otro es de datos burr.
set.seed(100)
#Vectores de datos
x <- rexp(500,rate = 5) 
z <- rburr(500, shape1=10, shape2=10, scale = 15)
ajuste(x,z)






