
#### Procedimento implementado

# 1) tomar diferencas nas covariáveis quando necessario

# 2) definir lags das variaveis de acordo com a maior correlacao

# 3) selecionar covariaveis com base no AICc

# 4) fazer analise de residuos



# rm(list = ls())
# setwd( "C:\\Users\\jafio\\OneDrive - unb.br\\Pesquisas\\Kinea" )


library(magrittr)
library(tidyverse)
library(readxl)
library(dplyr)
library(xts)
library(tseries)
library(forecast)
library(foreach)
library(norteste)


inicio <- Sys.time()

planilha <- read_excel("base.xlsx", sheet=1, col_names=TRUE)
planilha %>% head()
planilha %>% summary()

series_xts <- planilha[,-1] %>% xts(order.by = planilha$Dates)

series_xts[,1:11] %>% head()

series_xts[,1:11] %>% plot()

posicoes_NAs <- lapply(series_xts, FUN = function(x) which(is.na(x))) %>% unlist()

series_xts[posicoes_NAs, names(posicoes_NAs)]

series_xts[ , names(posicoes_NAs)] %>% plot()

series_xts[ , names(posicoes_NAs)] %>% na.aggregate() %>% plot()

# input
series_xts[ , names(posicoes_NAs)] <- series_xts[ , names(posicoes_NAs)] %>% na.aggregate()






series <- series_xts %>% ts(frequency=4)




#### estacionárias?

# raiz unitaria
raiz_unitaria <- which( apply( series, MARGIN=2, FUN=function(x) ndiffs( ts(x, frequency=4)) ) > 0 )

series[, names(raiz_unitaria)[2]] %>% plot()

series[, raiz_unitaria] <- rbind(NA, diff(series[,raiz_unitaria]))

colnames(series)[raiz_unitaria] <- paste0("Dif ", colnames(series)[raiz_unitaria] )

which( apply( series, MARGIN=2, FUN=function(x) ndiffs( ts(x, frequency=4)) ) > 0 )


# raiz unitaria sazonal
which( apply( series, MARGIN=2, FUN=function(x) nsdiffs( ts(x, frequency=4)) ) > 0 )





#### ccf
ccf_y <- foreach(i = 1:ncol(series), .combine=rbind) %do% {
  x <- series[,1]
  y <- series[,i]
  
  idx_NA <- which(is.na(y))
  if( length(idx_NA) > 0 ){
    if( max(idx_NA) > 1 ){
      x <- x[max(idx_NA):length(x)]
      y <- y[max(idx_NA):length(x)]
    }
  }  
  
  ccf_out <- ccf(x, y, lag.max=8, type = "correlation", na.action = na.omit, plot = FALSE)
  out <- ccf_out$acf
  tail(out,9) %>% round(3)
}

colnames(ccf_y) <- paste("lag", 0:(-8))
row.names(ccf_y) <- colnames(series)

 #cbind(series[,1], stats::lag(series[,2],-5)) %>% na.omit %>% cor( )

cc <- foreach(i = 1:nrow(ccf_y), .combine=rbind) %do% { c( max(ccf_y[i,]) , which.max(ccf_y[i,]));  }
colnames(cc) <- c("max","which.max")

ccf_y <- cbind( ccf_y, cc )

ccf_y %>% head()


#### series defasadas conforme maior correlacao (nao considera absoluto)
dados <- foreach(v = colnames(series), .combine=cbind) %do% { 
  stats::lag( series[,v], k= -(ccf_y[v,"which.max"]-1) )
}

colnames(dados) <- colnames(series) #c("y", paste0("x_",1:(ncol(series)-1)) )  
idx_NA <- dados[,1] %>% is.na %>% which()
dados <- dados[-idx_NA, ] %>% ts(frequency=4)

dados %>% head(9)
dados %>% tail(9)


#### Modelos #####################
dados[,1] %>% na.omit() %>% acf(main="ACF")
dados[,1] %>% na.omit() %>% pacf(main="PACF")

# rank correlacao
rank <- ccf_y[,"max"] %>% sort(decreasing = TRUE) %>% names()



# LM+ARMA --> selecao de variaveis por stepwise (forward)
max_var <- 10
var_selec <- rank[2]
fit <- Arima(dados[,1], xreg = dados[,var_selec], order=c(1,0,2), seasonal=c(0,0,0))
ic <- fit$aicc

for(i in 3:length(rank)){
  var_selec_new <- c(var_selec, rank[i])
  fit2 <- Arima(dados[,1], xreg = dados[, var_selec_new], order=c(1,0,2), seasonal=c(0,0,0))
  
  if( fit2$aicc < ic && all(fit2$coef[var_selec_new] > 0)){
    var_selec <- var_selec_new
    ic <- fit2$aicc
    fit <- fit2
  }
  
  if( length(var_selec) >= max_var )
    break
}

var_selec
ic
fit %>% summary()

fit$residuals %>% plot(main="residuos")
fit$residuals %>% na.aggregate() %>% kpss.test() # estacionariedade
fit$residuals %>% Box.test(lag=20, type="Ljung-Box")  # independencia
fit$residuals %>% na.omit() %>% acf(main="ACF residuos")
fit$residuals %>% na.omit() %>% pacf(main="PACF residuos")
fit$residuals %>% na.omit() %>% shapiro.test()

# fit$residuals %>% as.numeric() %>% na.omit() %>% qqnorm()
# fit$residuals %>% as.numeric() %>% na.omit() %>% qqline()

checkresiduals( na.omit(fit$residuals) )

ggplot(data.frame(fit$residuals), aes(sample = fit$residuals))+ stat_qq() + stat_qq_line()+theme_bw()


# modelo selecionado
lags <- (ccf_y[,"which.max"]-1)[var_selec]
EMV <- fit$coef
se <- diag(fit$var.coef)^(1/2)

nn <- length( setdiff(names(EMV), names(lags)) )

z0 <- EMV/se 
pvalor <- 2*pnorm(abs(z0), lower.tail = FALSE)

modelo <- data.frame(EMV=EMV, SE=round(se,3), Z = round(z0, 3), P.values = round(pvalor,3), nlags=c( rep("-", nn), lags))
#modelo <- data.frame(EMV=EMV, SE=round(se,3), nlags=c( rep("-", nn), lags))
modelo 



# fitted
fit_LM_SARMA <- xts(fit$fitted, order.by= tail(planilha$Dates, nrow(dados)))
plot( cbind(series_xts[,1], fit_LM_SARMA), main="observado, fit LM+SARMA" )

# fit LM e fit SARMA
fit_LM <- EMV["intercept"] + foreach(v = var_selec, .combine='+') %do% {  EMV[v] * dados[,v] }

fit_SARMA <- fit$fitted - fit_LM

fit_LM <- xts(fit_LM, order.by= tail(planilha$Dates, nrow(dados)))

fit_SARMA <- xts(fit_SARMA, order.by= tail(planilha$Dates, nrow(dados)))

plot( cbind(series_xts[,1], fit_LM), main="observado, fit LM" )

plot( cbind(series_xts[,1], fit_LM, fit_SARMA), main="observado, fit LM, fit SARMA" )



# nivel de correlacao entre as variaveis selecionadas
x <- as.matrix( dados[,var_selec] ) %>% na.aggregate() %>% na.omit()
x %>% dim()
colnames(x) <- paste0("x",1:10)
cor_x = cor(x) %>% round(3)
cor_x

# ordena correlacoes
cor_x[which(cor_x < 1)] %>% unique() %>% abs() %>% sort(decreasing = TRUE) #%>% head(10)



# tempo de processamento
fim <- Sys.time()
tempo = fim - inicio
tempo



