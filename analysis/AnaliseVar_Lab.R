## Analise de variancia Individual (Dados de laboratorio)

library("lme4")
library(viridis)
library(gridExtra)
library(sommer)
library(gt)
library(dplyr)
require(ggplot2)
library(egg)
library(sjstats)

#setwd("//Users//massainebandeiraesousa//EMBRAPA//Dados Natalia//Analises cap1")
dados <- read.table("Dados_Amilose_RVA_Natalia.txt", header=T, sep="\t", na.strings="NA")
dados <- dados[!dados$Rep==4,]
dados$Rep <- as.factor(dados$Rep)
dados$genotipo <- as.factor(dados$genotipo)
dados[,9:ncol(dados)] <- sapply(dados[,9:ncol(dados)],as.numeric)
dados$Umidade <- NULL

# Analise de variancia individual
traits <- colnames(dados[9:ncol(dados)])
parametros <- matrix(NA,ncol=6,nrow=length(traits))

## modelo (SOMMER)
#model1 <- mmer(fixed=Peak1 ~ 1,random = ~ genotipo, getPEV = T, data = dados)

model1 <- list()
resultados <- matrix(NA,length(traits),8)
r <- c(2,2,2,2,2,2,2,3)

# ## salvar os parametros
for(i in 1:8){
  model1 <- lm(dados[,8+i] ~ genotipo, data=dados)
  VarG <- (anova(model1)$`Mean Sq`[1] - anova(model1)$`Mean Sq`[2])/r[i]
  VarE <- anova(model1)$`Mean Sq`[2]
  H2 <- VarG / (VarG + (VarE/r[i])) #H2
 
  resultados[i,1] <- round(anova(model1)$`Mean Sq`[1],2)
  resultados[i,2] <- round(anova(model1)$`Mean Sq`[2],2)
  resultados[i,3] <- round(anova(model1)$`Pr(>F)`[1],4)
  resultados[i,4] <- round(VarG,2)
  resultados[i,5] <- round(H2,2)
  resultados[i,6] <- round(mean(dados[,8+i],na.rm=T),2) 
}

# Coef de Variação
resultados[1,7] <- round(cv(lm(Peak_1 ~ genotipo, data = dados)),2)
resultados[2,7] <- round(cv(lm(Trough_1 ~ genotipo, data = dados)),2)
resultados[3,7] <- round(cv(lm(Breakdown ~ genotipo, data = dados)),2)
resultados[4,7] <- round(cv(lm(Final_Visc ~ genotipo, data = dados)),2)
resultados[5,7] <- round(cv(lm(Setback ~ genotipo, data = dados)),2)
resultados[6,7] <- round(cv(lm(Peak_Time ~ genotipo, data = dados)),2)
resultados[7,7] <- round(cv(lm(Pasting_Temp ~ genotipo, data = dados)),2)
resultados[8,7] <- round(cv(lm(Amilose ~ genotipo, data = dados)),2)


## Tabela de resultados
colnames(resultados) <- c("QMGen","QMRes","Prob","VarGenetica","H2_amplo","Media","CV%","trait")
resultados[,8] <- traits
resultados <- as.data.frame(resultados)
resultados %>% gt()  

