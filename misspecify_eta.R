library(ZIM)
library(VGAM)
library(ggpubr)

rzinb.bbinom <- function(n,pi0,mu,size,eta,rho) {
  y <- rzinb(n,omega=pi0,k=size,lambda=mu)
  z <- rbetabinom(n,size=y,prob=eta,rho=rho)
  z
}


pi0 <- 0.1
psi <- 0.5
mu <- 30 # conside z_{ij}, here mu = s_j * mu_i
eta <- 0.1
rho <- c(0.001,0.1,0.2,0.5)

panel <-list()
for(i in 1:length(rho)) {
  K   <- c(0.5,0.75,1,1.33,2)
  set.seed(1)
  
  z <- lapply(c(3,1,2,4,5), function(j) {
    a <- (1-rho[i])/rho[i]*eta ; b <- (1-rho[i])/rho[i]*(1-eta)
    b2 <- b*K[j] ; a2 <- eta*K[j]*b/(K[j]*a-eta*a)*a
    #rzinb.bbinom(1000000,pi0=pi0,mu=mu*K[j],size=2,eta=eta/K[j],rho=min(rho[i]*(1-eta)/(K[j]-eta), 0.999))
    rzinb.bbinom(1000000,pi0=pi0,mu=mu*K[j],size=1/psi,eta=eta/K[j],rho=1/(a2+b2+1))
  })
  
  panel[[i]] <- 
    ggplot(data = data.frame(count = unlist(z),
                             eta = as.factor(rep(c('0.1', '0.2', '0.13', '0.08', '0.05'), each = 1000000))), 
           aes(count, fill = eta)) +
    geom_bar(aes(y=..prop..), position = position_dodge(), width = 0.8) + 
    scale_x_continuous(breaks=seq(0,15), limits=c(-0.5,15.5)) + ylab('Probability') +
    scale_fill_jco(name=expression(eta*minute)) +
    ggtitle(substitute(paste(rho, '=', rh), 
                       list(rh=round(rho[i], 2))))+
    theme_minimal()
}
ggarrange(panel[[1]], panel[[2]], panel[[3]], panel[[4]], nrow = 2, ncol = 2)

# 8 6
# 700 550