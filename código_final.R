#####################################################################################
######################### MODELO PROBITO ############################################
#Links relevantes:
#https://jbhender.github.io/Stats506/F18/GP/Group14.html Regressão probito e dados mroz
#http://galton.uchicago.edu/~eichler/stat24600/Handouts/s05add.pdf algoritmo EM

mroz = read.csv('https://vincentarelbundock.github.io/Rdatasets/csv/carData/Mroz.csv')

### Convertendo de "yes" "no" para 0,1
mroz$lfp = ifelse(mroz$lfp=="yes", 1, 0)
mroz$wc = ifelse(mroz$wc=="yes", 1, 0)
mroz$hc = ifelse(mroz$hc=="yes", 1, 0)

init = c(0,0,0,0)

### O correto seria pegar as váriaveis mais correlacionadas com lfp...
corr_mat = cor(mroz)

X = as.matrix(cbind(1,mroz[,c(3,5,8)]))
y = mroz[,2]

myprobit = glm(lfp ~ k5+age+lwg, family=binomial(link = "probit"), control=list(maxit=500, epsilon=1e-8), data = mroz)

## Só calculando pelo R para comparação com o modelo implementado abaixo
summary(myprobit)
coefsMAIN = coef(myprobit)

EM = function(params, X, y, tol=.0001, maxits = 500){
  b = params
  mu = X%*%b
  it = 0
  conv = FALSE
  z = rnorm(length(y)) #Variável latente Z
  
  while((!conv & (it<maxits))){
    z_ant = z
    
    z = ifelse(y==1, mu+dnorm(mu)/pnorm(mu), mu-dnorm(mu)/pnorm(-mu))     # Calculando Z|X
    
    b = solve(t(X)%*%X) %*% t(X)%*%z                                      # Resolvendo o Sistema que estima os parametros
    print(b)
    mu = X%*%b                                                            # Atualizando a estimativa de mu
    conv = max(abs(z_ant - z)) <= tol
    cat(paste(format(it), "\n", sep = ""))                                # Printando nº de iterações
    it = it + 1
  }
  return(b)
}

estimativa = EM(params=init, X=X, y=y, tol=1e-8, maxit=100)
#Erro entre a função do R e o EM
coefsMAIN - estimativa

pontos_pred = pnorm(estimativa[1]+estimativa[2]*mroz$age)

grafico = ggplot(data.frame(mroz$age, pontos_pred), aes(x=mroz.age)) +
  labs(title = "Regressão Probito: Lfp x Idade") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs( x = "Idade", y ="Probabilidade") +
  geom_line(aes(y = pontos_pred), col = "red", size = 0.75)

grafico

##################################################################################
############# Mistura Univariada e Bivariada de Normais ##########################
#Links relevantes:
#https://stats.stackexchange.com/questions/72774/numerical-example-to-understand-expectation-maximization
#https://www.stat.cmu.edu/~larry/all-of-statistics/=data/faithful.dat
#https://commons.wikimedia.org/wiki/File:Em_old_faithful.gif, código do R para caso bivariado.

seed = 250
flo = rnorm(50, 2, 1)
res = rnorm(75, 7, 3)

flores = sort(c(flo,res))

#Chutando valores iniciais para os parâmetros

mu_flo_guess = 1
sd_flo_guess = 2

mu_res_guess = 4
sd_res_guess = 4

pi_flo_guess = 0.5
pi_res_guess = 0.5

#####Funções

estim_mu = function(data, peso){
  return(sum(data*peso)/sum(peso))
}

estim_sd = function(data, peso, mean){
  variance = sum( (peso * (data - mean)**2) / sum(peso))
  return(sqrt(variance))
}

em_gaussian_mix = function(data_flores, mu_flo, mu_res, sd_flo, sd_res, pi_flo, pi_res, tol, maxit){
  erro = 1
  it = 0
  
  it_reservoir = data.frame( it, mu_flo, mu_res, sd_flo, sd_res, pi_flo, pi_res)
  print(it_reservoir)
  while( (erro > tol) & (it < maxit)){
    
    likelihood_of_flo = dnorm(data_flores, mu_flo, sd_flo)  #Chance de a obs x_i estar nessa distribuição
    likelihood_of_res = dnorm(data_flores, mu_res, sd_res)
    
    likelihood_total = (likelihood_of_flo) + (likelihood_of_res) #Normalizando para que ambos componentes somem 1
    
    peso_flo = likelihood_of_flo / likelihood_total #"peso" da verossimilhança para a distribuição 1, esses pesos são mto importantes na hora de estimar mu e sigma
    peso_res = likelihood_of_res / likelihood_total
    
    theta_anterior = c(mu_flo, mu_res, sd_flo, sd_res, pi_flo, pi_res)
    
    #Estimando correlação (talvez não faça muito sentido aqui)
    pi_flo = sum(peso_flo)/ length(data_flores)
    pi_res = sum(peso_res)/ length(data_flores)
    
    mu_flo = estim_mu(data_flores, peso_flo)
    mu_res = estim_mu(data_flores, peso_res)
    
    sd_flo = estim_sd(data_flores, peso_flo, mu_flo)  
    sd_res = estim_sd(data_flores, peso_res, mu_res)
    
    theta = c(mu_flo, mu_res, sd_flo, sd_res, pi_flo, pi_res)
    erro = abs(theta_anterior - theta)
    erro = max(erro)
    it = it + 1
    it_reservoir = rbind.data.frame(it_reservoir, c(it, theta))
  }
  cat("Número de Iterações até convergência:", it)  
  return(it_reservoir)
}

tol = 10^(-8)
maxit = 500

#Chamando a func

em_flores = em_gaussian_mix(flores, mu_flo_guess, mu_res_guess, sd_flo_guess, sd_res_guess, pi_flo_guess, pi_res_guess, tol, maxit)

dens = function(x){
  return(pi_flo_estimado*dnorm(x, mu_flo_estimado, sd_flo_estimado) + pi_res_estimado*dnorm(x, mu_res_estimado, sd_res_estimado))
}

d_flo = function(x){
  return(pi_flo_estimado*dnorm(x, mu_flo_estimado, sd_flo_estimado))
}

d_res = function(x){
  return(pi_res_estimado*dnorm(x, mu_res_estimado, sd_res_estimado))
}

library(ggplot2)

for ( i in 1:length(em_flores[,1])){
  
  #Criando a animação
  mu_flo_estimado = as.vector(em_flores[,2])[i]
  mu_res_estimado = as.vector(em_flores[,3])[i]
  sd_flo_estimado = as.vector(em_flores[,4])[i]
  sd_res_estimado = as.vector(em_flores[,5])[i]
  pi_flo_estimado = as.vector(em_flores[,6])[i]
  pi_res_estimado = as.vector(em_flores[,7])[i]
  
  if (i < 10) {name = paste('000',i,'plot.png',sep='')} 
  if (i < 100 && i >= 10) {name = paste('00',i,'plot.png', sep='')}
  if (i >= 100) {name = paste('0', i,'plot.png', sep='')}

  png(name)
  grafico_it = ggplot(data.frame(flores), aes(x=flores)) +
    labs(title = "Histograma Mistura de Normais") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs( x = "x", y ="Densidade") +
    geom_histogram(aes(y=..density..), binwidth=0.5, fill = "lightsteelblue", col="blue") +
    stat_function(fun = dens, col="orange", xlim = c(min(flores)-0.5, max(flores)), size = 0.75) +
    stat_function(fun = d_flo, col="darkblue", xlim = c(min(flores)-0.5, max(flores))) + 
    stat_function(fun = d_res, col = "red", xlim = c(min(flores)-0.5, max(flores)))
  print(grafico_it)  
  dev.off()
}



#Recapturando os valores
mu_flo_estimado = em_flores[2][length(em_flores[,1]),]
mu_res_estimado = em_flores[3][length(em_flores[,1]),]
sd_flo_estimado = em_flores[4][length(em_flores[,1]),]
sd_res_estimado = em_flores[5][length(em_flores[,1]),]
pi_flo_estimado = em_flores[6][length(em_flores[,1]),]
pi_res_estimado = em_flores[7][length(em_flores[,1]),]

grafico = ggplot(data.frame(flores), aes(x=flores)) +
  labs(title = "Histograma Mistura de Normais") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs( x = "x", y ="Densidade") +
  geom_histogram(aes(y=..density..), binwidth=0.5, fill = "lightsteelblue", col="blue") +
  stat_function(fun = dens, col="orange", xlim = c(min(flores)-0.5, max(flores)), size = 0.75) +
  stat_function(fun = d_flo, col="darkblue", xlim = c(min(flores)-0.5, max(flores))) + 
  stat_function(fun = d_res, col = "red", xlim = c(min(flores)-0.5, max(flores)))
grafico

library(mvtnorm)


############################### Código Bivariado ###############################



#Dados de erupções do Geyser Old Faithful
data(faithful)

ggplot(faithful, aes(x= x)) + 
  geom_point(data = faithful,
  aes(x = eruptions, y = waiting, z = NULL), color = "red") +
  labs(title = "Conjunto de Dados Old Faithful") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs( x = "Duração da Erupção(min)", y ="Tempo de Espera(min)") 

#Estimativas inicias dos parâmetros
theta = list(
  tau = c(0.5,0.5), #Correlação de ambas normais
  mu1 = c(2.8,75),  #Médias 
  mu2 = c(3.6,58),
  sigma1 = matrix(c(0.8,7,7,70),ncol=2), #Matrizes de variâncias e covariâncias
  sigma2 = matrix(c(0.8,7,7,70),ncol=2)
)


#Passo E: Calcular probabilidade condicional das variáveis latentes
passo_e <- function(theta) {
  t(apply(cbind(
    theta$tau[1] * dmvnorm(faithful,mean=theta$mu1,sigma=theta$sigma1),
    theta$tau[2] * dmvnorm(faithful,mean=theta$mu2,sigma=theta$sigma2)
  ),1,function(x) x/sum(x))) }


#Passo M = maximizar a função de verossimilhança
passo_m <- function(T) list(
  tau = apply(T,2,mean),
  mu1 = apply(faithful,2,weighted.mean,T[,1]),
  mu2 = apply(faithful,2,weighted.mean,T[,2]),
  sigma1 = cov.wt(faithful,T[,1])$cov,
  sigma2 = cov.wt(faithful,T[,2])$cov)

tol_nbv = 0.0001
erro = 1
iter = 0
maxit = 500

while ((iter < maxit) && (erro > tol_nbv)){
  theta_old = theta
  T <- passo_e(theta)
  theta <- passo_m(T)
  iter = iter + 1
  erro_corr = max(abs(theta[['tau']] - theta_old[['tau']]))
  erro_mu1 = max(abs(theta[['mu1']] - theta_old[['mu1']]))
  erro_mu2 = max(abs(theta[['mu2']] - theta_old[['mu2']]))
  erro_sd1 = max(abs(theta[['sigma1']] - theta_old[['sigma1']]))
  erro_sd2 = max(abs(theta[['sigma2']] - theta_old[['sigma2']]))
  erro = max(erro_corr,erro_mu1,erro_mu2,erro_sd1,erro_sd2)
}

###Printando o resultado final:

plot.em <- function(theta){
  mixture.contour <- outer(xpts,ypts,function(x,y) {
    theta$tau[1]*dmvnorm(cbind(x,y),mean=theta$mu1,sigma=theta$sigma1) + theta$tau[2]*dmvnorm(cbind(x,y),mean=theta$mu2,sigma=theta$sigma2)
  })
  contour(xpts,ypts,mixture.contour,nlevels=5,drawlabel=FALSE,col="red",xlab="Eruption time (mins)",ylab="Waiting time (mins)",main="Waiting time vs Eruption time of the Old Faithful geyser")
  points(faithful)
}

plot_em_ggplot = function(theta){
  z = outer(xpts,ypts,function(x,y) {
    theta$tau[1]*dmvnorm(cbind(x,y),mean=theta$mu1,sigma=theta$sigma1) + theta$tau[2]*dmvnorm(cbind(x,y),mean=theta$mu2,sigma=theta$sigma2)
  })
  
  library(reshape2) 
  z.molten = melt(z) 
  names(z.molten) = c("x", "y", "z")
  #print(z.molten)
  
  v = ggplot(z.molten, aes(x,y,z=z))+
    stat_contour(bins=6,aes(x,y,z=z), color="black", size=0.6) 
  
  v  + geom_point(data = faithful,
               aes(x = ind, y = waiting, z = NULL), color = "red") +
  labs(title = "Curvas de Nível - Old Faithful") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs( x = "Duração da Erupção(min)", y ="Tempo de Espera(min)") + 
  scale_x_continuous( breaks = c(20,40,60,80), labels = c(2,3,4,5)) 
}

xpts <- seq(from=1,to=6, length.out = 100)
ypts <- seq(from=1,to=100, length.out = 100)

ind = numeric(272)
eru = faithful$eruptions

for(i in 1:(length(xpts)-1)){
  for(j in 1:272){
    if (eru[j] >= xpts[i] && eru[j] <= xpts[i+1]){
      ind[j] = i
    }
  }
}

#data("faithful")
faithful = cbind.data.frame(faithful, ind)
plot_em_ggplot(theta)

#####################################################################################################
##################################################### Classificação EM vs Kmeans 
###### K Means LOL https://www.datanovia.com/en/blog/k-means-clustering-visualization-in-r-step-by-step-guide/
###### https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/kmeans

class = kmeans(faithful, 2)

grupo = class[["cluster"]]
df_class = data.frame(faithful, grupo)

graf_kmeans = ggplot(df_class, aes(x = eruptions,y = waiting))
graf_kmeans + geom_point(data = df_class,
                       aes(x = eruptions, y = waiting, colour = grupo)) +
              geom_point(data.frame(class$center), mapping = aes(x = eruptions, y = waiting), shape = 10, size = 10,col = c("red", "orange"))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + labs(title = "Classificação via K means") +
  labs( x = "Duração da Erupção(min)", y ="Tempo de Espera(min)") 


################################################## Classificação via MClust.
##### EM Mclust https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html

library(mclust)
plot(faithful)

faithfulBIC <- mclustBIC(faithful)
plot(faithfulBIC)

mod1 <- Mclust(faithful, G=2, x = faithfulBIC)

df_em_class = data.frame(faithful, mod1[["classification"]])

summary(mod1, parameters = TRUE)

df_em_class$mod1...classification... = ifelse(df_em_class$mod1...classification... == 1, 2, 1)

params = data.frame(c(4.291633,2.038831), c(79.990255, 54.506423))
colnames(params) = c("E", "W")

graf_em_class = ggplot(df_em_class, aes(x = eruptions,y = waiting))
graf_em_class + geom_point(data = df_em_class,
                         aes(x = eruptions, y = waiting, colour = mod1...classification...)) +
                theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + labs(title = "Classificação via E.M.") +   
  labs( x = "Duração da Erupção(s)", y ="Tempo de Espera") +
  geom_point( data = params, mapping = aes( x= E, y=W), shape = 10, size = 10, col = c("orange", "red"))


#Gráfico em 3d
mod5 <- densityMclust(faithful[,2:3], G = 2)
plot(mod5, main="Densidade da Mistura", xlab = "Duração da Erupção(min)", ylab ="Tempo de Espera(min)", what = "density", type = "persp", theta = -10)
