

datos<-read.csv("C:\\Users\\Administrator\\Desktop\\casosCOVID.csv",sep = ";",dec='.',header=TRUE)
h <- hist(datos$casos, probability = T, main = "Normalidad Covid19 gomez",xlab = "", ylab = "",col = "red")
lines(density(datos$casos,na.rm = T), lwd = 2, col = "green")
mu <- mean(datos$casos, na.rm = T)
sigma <- sd(datos$casos,na.rm = T)
x <- seq(min(h$mids,na.rm = T), max(h$mids,na.rm = T), length = length(datos$casos))
y <- dnorm(x, mu, sigma)
lines(x,y,lwd =2, col = "blue")


qqnorm(datos$casos)
qqline(datos$casos,col = "red")


library(nortest)
nortest::pearson.test(datos$casos)
nortest::lillie.test(datos$casos)
nortest::cvm.test(datos$casos)


airpass<-ts(datos$casos,start=c(2020,66),freq=366)
plot(airpass,type="o",col="red")
plot(diff(airpass),type="o") # Muestra posibles atipicidades

casos_covid<-ts(datos$casos,start=c(2020,66),freq=366)
res<-spec.pgram(casos_covid,log="no")


order(res$spec,res$freq,decreasing = TRUE)

max<-res$freq[2]
max2<-res$freq[4]
max3<-res$freq[6]

periodo1<-366/max
periodo2<-366/max2
periodo3<-366/max3

abline(v=max,lty="dotted",col="red")
abline(v=max2,lty="dotted",col="blue")
abline(v=max3,lty="dotted",col="green")

cleandata<-read.csv("C:\\Users\\Administrator\\Desktop\\Cleaned-Data.csv",sep = ",",dec='.',header=TRUE)

library(corrplot)
cleandata$Country<-NULL
correlaciones <- cor(cleandata)
corrplot(correlaciones, method='circle', shade.col=NA, tl.col='black',
         tl.srt=40, addCoef.col='black', order='AOE', type = 'lower')

corrplot(correlaciones)

corrplot(cor(correlaciones), diag = FALSE, order = "FPC",
         tl.pos = "td", tl.cex = 0.5, method = "color", type = "upper")

corrplot(correlaciones, type="lower")



w<-datos$casos
x<-w
for(t in 2:length(datos$casos)){
  x[t]<-x[t-1]+w[t]
} 

par(mfrow=c(2,2))
ts.plot(x, main="Camino aleatorio X_t")
acf(x, main="Autocorrelaci?n Simple de Xt",ylim=c(-1,1),col="black")
d<-diff(x)
ts.plot(d,main="Primera diferencia de X_t")
acf(d, main="Autocorrelaci?n Simple de d",ylim=c(-1,1),ci.col="black")


library(nortest)
library(forecast)

casosARMA<-ts(datos$casos,start=c(2020,66),freq=366)


nortest::pearson.test(datos$casos)
nortest::lillie.test(datos$casos)
nortest::cvm.test(datos$casos)

acf(casosARMA, main="Autocorrelaci?n Simple",col="black",ylim=c(-1,1))
pacf(casosARMA,main="Autocorrelaci?n Parcial",col="black",ylim=c(-1,1))

fit<-arima(casosARMA,order=c(0,2,1))

auto.arima(casosARMA)

res<-spec.pgram(casosARMA,log="no")
order(res$spec,res$freq,decreasing = TRUE) 
max1<-res$freq[2]
max2<-res$freq[4]
max3<-res$freq[6]
periodo1 <- 365/max1
periodo2 <- 365/max2
periodo3 <- 365/max3

fit<-arima(casosARMA,order=c(0,2,1),seasonal=list(order=c(1,1,1),period=3))
cA.pred<-predict(fit,n.ahead=30)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(casosARMA,xlim=c(2020.1,2020.365),ylim=c(0,600),type="o")
lines(cA.pred$pred,col="red",type="o")
lines(cA.pred$pred+2*cA.pred$se,col="red",lty=3,type="o")
lines(cA.pred$pred-2*cA.pred$se,col="red",lty=3,type="o")
acf(LakeHuron, main="Autocorrelaci?n Simple",col="black",ylim=c(-1,1))
pacf(LakeHuron,main="Autocorrelaci?n Parcial",col="black",ylim=c(-1,1))




suppressMessages(library(xts))
suppressMessages(library(dygraphs))
suppressMessages(library(itsmr))
suppressMessages(library(forecast))  
#fit<-arima(casosARMA,order=c(0,2,1),seasonal=list(order=c(0,0,1),period=2))
#fit<-arima(x = casosARMA, order = c(1, 2, 2), seasonal = list(order = c(0, 0, 0), period = 2))
fit<-arima(x = casosARMA, order = c(1, 0, 2))
#fit<-arima(casosARMA,order=c(1,1,1))
cA.pred<-predict(fit,n.ahead=15)
cA.preds<-cA.pred$pred
LimInf<-cA.preds-cA.pred$se
LimSup<-cA.preds+cA.pred$se
per_1<-seq(as.Date("2020-03-06"),as.Date("2020-04-13"),"day")
per_2<-seq(as.Date("2020-04-14"),as.Date("2020-04-28"),"day")
todas.series<-cbind(casosARMA=xts(casosARMA,order.by = per_1),LimInf=xts(LimInf,order.by = per_2),Pronostico=xts(cA.preds,order.by = per_2),LimSup=xts(LimSup,order.by = per_2))
dygraph(todas.series,main="Casos COVID19",ylab="Cantidad de CASOS ")%>%
  dySeries(c("LimInf", "Pronostico", "LimSup"), label = "Casos")%>%
  dyRangeSelector(height = 20, strokeColor = "")%>%  
  dyOptions(axisLineColor = "navy", 
            gridLineColor = "lightblue")


library(lubridate)
library(tidyverse)
library(forecast)
library(dygraphs)
library(xts)


calibrar.arima <- function(entrenamiento = NULL, prueba = NULL, periodo = NA_integer_, rango = 0:2) {
  # se calculan todas las combinaciones para los parametros
  params <- cross(list(a = rango, b = rango, c = rango,
                       d = 0:1, e = 0:1, f = 0:1))
  # se calcula un modelos para cada combinacion de parametros
  arima_secure <- possibly(stats::arima, otherwise = NULL)
  models <- map(params, ~ suppressWarnings(arima_secure(entrenamiento, order = c(.$a,.$b,.$c),
                                                        seasonal = list(order = c(.$d,.$e,.$f),
                                                                        period = periodo))))
  # se eliminan los modelos fallidos
  models <- keep(models, negate(is.null))
  # se hace una prediccion para cada modelos
  predictions <-map(models, ~predict(., n.ahead = length(prueba)))
  # se calcula el error para cada prediccion
  error <- map_dbl(predictions, function(pred, real) {
    error <- sum((as.numeric(real) - as.numeric(pred$pred))^2)
    return(error)
  }, real = prueba)
  
  # se retorna el modelo con el menor error
  best_model <- models[[which.min(error)]]
  p <- params[[which.min(error)]]
  best_model$call <- call("arima",
                          x = quote(datos), order = as.numeric(c(p$a, p$b, p$c)),
                          seasonal = list(order = as.numeric(c(p$d, p$e, p$f)),
                                          period = periodo))
  return(best_model)
}


train.covid<-datos$casos[1:22]
test.covid<-datos$casos[23:nrow(datos)]

train.serie.covid <- ts(train.covid, start = c(2020, 66), frequency = 366)

calibrar.arima(train.covid, test.covid)


install.packages("FactoMineR") 
# Para generar gr?ficos alternativos
install.packages("factoextra")
library("FactoMineR") 
library("factoextra")

mi.tema <- theme_grey() + theme(panel.border = element_rect(fill = NA,color = "white"), plot.title = element_text(hjust = 0.5))


datosUSA<-read.csv("C:\\Users\\Administrator\\Desktop\\casosCOVIDUSA.csv",sep = ",",dec='.',header=TRUE)

rownames(datosUSA)<-datosUSA$State
datosUSA2<-datosUSA[-1]
datosUSA2<-datosUSA2[-25]

modelo <- prcomp(datosUSA2,scale. = TRUE,center = TRUE)

plot(modelo, axes=c(1, 2), choix="ind", col.ind="red",new.plot=TRUE)
plot(modelo, axes=c(1, 2), choix="var", col.var="blue",new.plot=TRUE)
fviz_pca_biplot(axes = c(1,2),modelo,col.var = "#2E9FDF",col.ind = "#696969",ggtheme = mi.tema)
fviz_eig(modelo)


