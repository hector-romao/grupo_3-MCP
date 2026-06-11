
library(ape)
library(phytools)
library(adephylo)
library(phylobase)
library(picante)
library(geiger)
library(TreeSim)
library(phylosignal)
library(caper)
library(phylolm)
library(here)


###Importação dos dados

#Importação da árvore
phy <-read.tree(here::here("data/insecta_tree.txt"))
plot(phy)
is.ultrametric(phy)

#Importação da tabela com os tratis analisados
traits <-read.table("data/insecta_traits.tsv",h=T,row.names=1)
head(traits)

###Preparo dos objetos para as analises
c_type_count <-as.matrix(traits[,1])
head(c_type_count)
i_type_count <-as.matrix(traits[,2])
head(i_type_count)


rownames(c_type_count) <-phy$tip.label #para ape/geiger...
rownames(i_type_count) <-phy$tip.label #para ape/geiger...
head(c_type_count)
head(i_type_count)

#PHYLOGENETIC SIGNAL
#MODEL-BASED METHODS

#BLOMBERG's K
plot(phy)
traits
names(c_type_count) <-phy$tip.labels
names(i_type_count) <-phy$tip.labels
names(range) <-phy$tip.labels

log_c_type_count <- log(c_type_count[,1]+1)
names(log_c_type_count) <- rownames(c_type_count)

log_i_type_count <- log(i_type_count[,1]+1)
names(log_i_type_count) <- rownames(i_type_count)

k_obs_c_type <- Kcalc(log_c_type_count, phy)
k_obs_i_type <- Kcalc(log_i_type_count, phy)

#A note on randomization of traits for analysis with ape 
#Species names will be randomized with traits, so match will eliminate randomizations. So,...

bs <-as.matrix(traits[,1]) #Obtenção do trait do C-type como matriz 
bs <-sample(bs,279,r=F) #Amostragem aleatoria do conjunto de dados
names(bs) <-phy$tip.label

#Visualização da randomização
traits[, 1, drop = FALSE]
print(bs)

Kbm <-numeric()
Krand <-numeric()
Knorm <-numeric ()

for(i in 1:10){
  
  T.brownian <- fastBM(phy, a = 10, sig2 = 10, internal = F, nsim = 1)
  Kbm[i] <-Kcalc(T.brownian,phy)
  
  y.norm <-as.matrix(rnorm(nrow(traits),0,1))
  rownames(y.norm) <- phy$tip.label
  Knorm[i] <-Kcalc(y.norm[,1],phy)
  
  bs <-as.matrix(traits[,1])
  bs <-sample(bs,279,r=F)
  names(bs) <-phy$tip.label
  y.rand <-bs
  Krand[i] <-Kcalc(y.rand,phy)
  
}
hist(Kbm,nclass=50,col="grey",xlab="Blomberg's K under Brownian motion",main="")
mean(Kbm)
quantile(Kbm,c(0.025,0.975))

sum(ifelse(Kbm > k_obs_c_type[,1],1,0)) #Quantas simulações de BM produziram um K maior que o observado

mean(Krand)
hist(Krand,nclass=100,col="grey",xlab="Blomberg's K under randomization",main="")
quantile(Krand,0.95)

mean(Knorm)
hist(Knorm,nclass=50,col="grey",xlab="Blomberg's K under randomization",main="")
quantile(Knorm,0.95)

#################################################################################
#Alternatva com phytools
#################################################################################
#Phytools...Blomberg's K and Pagel's Lambda

phylosig(phy,log_c_type_count,method="K",test=T,nsim=1000)
quantile(Krand,c(0.05,0.9995))#check with previous simulations simulations...

phylosig(phy,log_c_type_count,method="lambda",test=T) #likelihood test for lambda=0


#Também podemos fazer a analise para o i_type
#phylosig(phy,log_i_type_count,method="K",test=T,nsim=1000)
#quantile(Krand,c(0.05,0.9995))#check with previous simulations simulations...

#phylosig(phy,log_i_type_count,method="lambda",test=T) #likelihood test for lambda=0



#################################################################################
#"STATISTICAL" EMPIRICAL METHODS
#################################################################################

#Moran
phy.cor <-vcv(phy,model="Brownian",cor=T)
diag(phy.cor) <-0

Moran.I(log_i_type_count,phy.cor) #i_type_count
Moran.I(log_c_type_count,phy.cor) #c_type_count


sim.trait <- fastBM(phy, a = 10, sig2 = 10, internal = F, nsim = 100) #Simula a evolução de um traço sob BM


#Variaveis
I.sim <-numeric() #Sob BM
I.simP <-numeric()
I.sim2 <-numeric() #Aleatorios
I.sim2P <-numeric()
I.sim3 <-numeric() #Permutados

for(i in 1:ncol(sim.trait)){ #Para cada coluna (simulação) 
  I.sim[i] <-Moran.I(sim.trait[,i],phy.cor)$observed
  I.simP[i] <-Moran.I(sim.trait[,i],phy.cor)$p.value
  I.sim2[i] <-Moran.I(rnorm(nrow(traits),0,1),phy.cor)$observed
  I.sim2P[i] <-Moran.I(rnorm(nrow(traits),0,1),phy.cor)$p.value
  rnd.i_type_count <- sample(log_c_type_count, replace=F)
  I.sim3[i] <-Moran.I(rnd.i_type_count,phy.cor)$observed
}

###Observação dos casos
hist(I.sim,nclass=100,col="grey",xlab="Moran's I under BM",ylab="simulations",main="") #Brownian expectation
hist(I.sim2,nclass=100,col="grey",xlab="Moran's I (null normal)",ylab="simulations",main="") #Null distribution
hist(I.sim3,nclass=100,col="grey",xlab="Moran's I (null randomized)",ylab="simulations",main="") #Null distribution

length(which(I.simP < 0.05)) #Quantas simulações foram significativas sob BM 
length(which(I.sim2P < 0.01)) #Quantas simulações foram significativas sob aleatoriedade


mean(I.sim2)
median(I.sim2)
sum(ifelse(I.sim > 0.919,1,0))/1000

quantile(I.sim2,c(0.025,0.975))
quantile(I.sim3,c(0.025,0.975))


#mudando o sim.trait o que ocorre...?

#A very simple correlogram...
dist <-as.matrix(cophenetic(phy))
hist(dist)
klim <- c(0,50,100,200,400,600,800,1030) #Recortes de tempo

Im <-numeric()
Ip <-numeric()

for(i in 1:(length(klim)-1)){
  W <-ifelse(dist > klim[i] & dist < klim[i+1],1,0)
  diag(W) <-0
  I <-Moran.I(log_c_type_count,W)
  Im[i] <-I$observed
  Ip[i] <-I$p.value
}
plot(klim[2:8],Im,type='b')


#Simulations of Moran's I under OU
phy.cor <-vcv(phy,cor=T)
diag(phy.cor) <-0

alpha.ou <-numeric()
I.simou <-numeric()

phy.sim <-rcoal(279)
plot(phy.sim)
phy.cor.sim <-vcv(phy.sim,cor=T)

for(i in 1:100){
  alpha.ou[i] <-runif(1,0.01,25)
  phy2 <-rescale(phy.sim,model="depth",1) 
  sim.trait <- fastBM(phy2, a = 10, alpha = alpha.ou[i], theta=5, internal = F, nsim = 1)
  I.simou[i] <-Moran.I(sim.trait,phy.cor.sim)$observed
}
plot(alpha.ou,I.simou,cex=1,pch=16,xlab="Alpha",ylab="Moran's I")


#PVR & PSR
plot(phy)
dist <-as.matrix(cophenetic(phy))
dist <-sqrt(dist) #para PSR, deixando linear com distancias (BM)
pcoa <-pcoa(dist)
vec <-pcoa$vectors
val <-pcoa$values
eigcum <-val[1:279,4]

#3 axes
phy_pvr1 <-phylo4d(phy,vec[,1:3]) #i_type_count size and c_type_count size
table.phylo4d(phy_pvr1, show.tip.label = FALSE,center=FALSE, ratio.tree = 0.7, box = FALSE, cex.symbol = 0.4, cex.label = 1)

#Plottando somente os afideos
clado_HU <- extract.clade(phy, node = 297)
phy_pvr1 <-phylo4d(clado_HU,vec[,1:3])
table.phylo4d(phy_pvr1, show.tip.label = FALSE,center=FALSE, ratio.tree = 0.7, box = FALSE, cex.symbol = 0.4, cex.label = 1)


#3 methods for eigenvector selection baseados nos valores de r2 de acordo com a seleção
#selecting eivectors with eigenvalues > 0.30
vec30<-min(which(eigcum > 0.55)) #Pegando somente os eixos que juntam explicam 55 da variação
summary(lm(log_c_type_count~vec[,1:vec30]))
summary(lm(log_c_type_count~vec[,1:3])) # somente os 3 primeiros os  significativos
S.PVR <-lm(log_c_type_count~vec[,1:vec30])$residuals #para evolu??o correlacionada...
S.PVR <-lm(log_c_type_count~vec[,1:4])$residuals #para evolu??o correlacionada...

S.PVR <-lm(log_c_type_count~vec[,1:vec30])$residuals #Observação se os residuos possuem sinal filogenetico
Moran.I(S.PVR,phy.cor)


#stepwise
summary(step(lm(log_c_type_count~vec[,1:min(which(eigcum > 0.55))-1]),direction="forward"))
summary(step(lm(log_i_type_count~vec[,1:min(which(eigcum > 0.55))-1]),direction="forward"))

#PVR with eigevector selection by significant correlations with Y
pv <-numeric()
for(i in 1:ncol(vec)){
  pv[i] <-summary(lm(log_c_type_count~vec[,i]))$coefficients[2,4]
  #se quiser usar o cor.test...
}
selpvr <-which(pv < 0.05)
vec_sel <-vec[,selpvr]
summary(lm(log_c_type_count~vec_sel))
S.PVR1 <-lm(log_c_type_count~vec_sel)$residuals #para evolu??o correlacionada...

Moran.I(S.PVR1,phy.cor)

hist(phy.cor)
phy.cor1 <-ifelse(phy.cor > 0.75,1,0)
diag(phy.cor1) <-1
Moran.I(S.PVR,phy.cor1)


# Seleção de autovetores minimizando o I de Moran residual
forward_moran_selection <- function(y, X, W, threshold,
                                    use_abs = TRUE,
                                    verbose = TRUE,
                                    na_action = na.omit) {
  
  if (!is.data.frame(X)) {
    X <- as.data.frame(X)
  }
  
  # Junta tudo para tratar NA de forma consistente
  dat <- data.frame(y = y, X)
  dat <- na_action(dat)
  
  y <- dat$y
  X <- dat[, setdiff(names(dat), "y"), drop = FALSE]
  
  selected <- character(0)
  remaining <- names(X)
  history <- list()
  step <- 0
  
  repeat {
    if (length(remaining) == 0) {
      if (verbose) message("Nenhum preditor restante.")
      break
    }
    
    results_step <- data.frame(
      predictor = remaining,
      moran_I = NA_real_,
      criterion = NA_real_,
      stringsAsFactors = FALSE
    )
    
    models_step <- vector("list", length(remaining))
    
    for (i in seq_along(remaining)) {
      pred <- remaining[i]
      vars <- c(selected, pred)
      
      form <- as.formula(
        paste("y ~", paste(vars, collapse = " + "))
      )
      
      fit <- lm(form, data = data.frame(y = y, X))
      res <- resid(fit)
      
      I_obs <- Moran.I(res, W)$observed
      crit <- if (use_abs) abs(I_obs) else I_obs
      
      results_step$moran_I[i] <- I_obs
      results_step$criterion[i] <- crit
      models_step[[i]] <- fit
    }
    
    # escolhe o menor critério
    best_idx <- which.min(results_step$criterion)
    best_pred <- results_step$predictor[best_idx]
    best_I <- results_step$moran_I[best_idx]
    best_crit <- results_step$criterion[best_idx]
    best_model <- models_step[[best_idx]]
    
    step <- step + 1
    selected <- c(selected, best_pred)
    remaining <- setdiff(remaining, best_pred)
    
    history[[step]] <- list(
      step = step,
      tested = results_step,
      chosen = best_pred,
      moran_I = best_I,
      criterion = best_crit,
      model = best_model,
      selected_so_far = selected
    )
    
    if (verbose) {
      message(
        sprintf(
          "Passo %d | Escolhido: %s | Moran.I = %.6f | Critério = %.6f",
          step, best_pred, best_I, best_crit
        )
      )
    }
    
    if (best_crit <= threshold) {
      if (verbose) message("Threshold atingido. Encerrando seleção.")
      break
    }
  }
  
  final_model <- if (length(history) > 0) history[[length(history)]]$model else NULL
  
  list(
    selected = selected,
    X_selected = X[, selected, drop = FALSE],
    data_used = dat
  )
}  


#Uso da função de seleção...
out <- forward_moran_selection(
  y = log_c_type_count,
  X = as.data.frame(vec),
  W = phy.cor1,
  threshold = 0.15,   # exemplo
  use_abs = TRUE
)

vec.opt <-as.matrix(out$X_selected) #vetores para o PVR otimizado

summary(lm(log_c_type_count~vec.opt))



#Under Brownian motion...
sim.trait <- fastBM(phy, a = 10, sig2 = 10, internal = F, nsim = 1000) #brownian...

r2.pvr <-numeric()
for(i in 1:ncol(sim.trait)){
  r2.pvr[i] <-summary(lm(sim.trait[,i]~vec[,1:2]))$r.squared
}
hist(r2.pvr,nclass=50,col="grey",main="",xlab="R2 (PVR - 2 axes)",ylab="Simulations")
abline(v=0.98,col="red",lwd=2)
median(r2.pvr)



#PSR Curve
r2psr <-numeric()
r2bw <-numeric()
difr2 <-numeric()

for(i in 1:(nrow(dist)-1)){
  r2psr[i] <-summary(lm(log_c_type_count~vec[,1:i]))$r.squared
}

eig <-append(0,eigcum[-10])
r2psr0 <-append(0,r2psr)
plot(eig,r2psr0,pch=16,cex=1.25,type="b",xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1,lty=2,lwd=2)



phy.homo <-phy
phy.homo$edge.length[7] <- (phy.homo$edge.length[7] *  6)
phy.homo$edge.length[2] <- (phy.homo$edge.length[2] *  15)
plot(phy.homo)
dist <-as.matrix(cophenetic(phy.homo))
dist <-sqrt(dist)
pcoa <-pcoa(dist)
vec <-pcoa$vectors
val <-pcoa$values
eigcum <-val[1:279,4]

#PSR Curve
r2psr <-numeric()
r2bw <-numeric()
difr2 <-numeric()

for(i in 1:(nrow(dist)-1)){
  r2psr[i] <-summary(lm(log_c_type_count~vec[,1:i]))$r.squared
}

eig <-append(0,eigcum[-10])
r2psr0 <-append(0,r2psr)
plot(eig,r2psr0,pch=16,cex=1.25,type="b",xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1,lty=2,lwd=2)



#PSR Brownian
brown_r2 <-matrix(0,ncol(sim.trait),length(r2psr))
sim.trait <- fastBM(phy, a = 10, sig2 = 10, internal = F, nsim = 10) #brownian...

for(i in 1:ncol(sim.trait)){
  
  brown <-sim.trait[,i]
  for(j in 1:(nrow(dist)-1)){
    r2bw[j] <-summary(lm(brown~vec[,1:j]))$r.squared
    r2psr[j] <-summary(lm(log_c_type_count~vec[,1:j]))$r.squared
  }
  brown_r2[i,] <-r2bw
  
}
meanBW <-apply(brown_r2,2,mean)
maxBW <-apply(brown_r2,2,max)
minBW <-apply(brown_r2,2,min)

meanBW <-append(0,meanBW)
maxBW <-append(0,maxBW)
minBW <-append(0,minBW)

plot(eig,r2psr0,pch=16,cex=1.25,type="b",xlim=c(0,1),ylim=c(0,1))
lines(eig,meanBW,pch=16,cex=1.25,type="b",xlim=c(0,1),ylim=c(0,1),col="blue")
lines(eig,maxBW,pch=16,cex=1.25,type="l",xlim=c(0,1),ylim=c(0,1),col="red")
lines(eig,minBW,pch=16,cex=1.25,type="l",xlim=c(0,1),ylim=c(0,1),col="red")
abline(a=0,b=1,lty=2)

#Non-stationarity (deviations from Brownian)
Pnst <-numeric()
for(i in 1:(ncol(brown_r2)-1)){
  v <-brown_r2[,i+1] - brown_r2[,i]
  dr2 <-r2psr0[i+1]-r2psr0[i]
  Pnst[i] <-(sum(ifelse(which(v > dr2),1,0)))/ncol(sim.trait)
}
r2psr
Pnst



