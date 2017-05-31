###################################################################
##########       Solución numérica - Put Americana       ##########
##########          Diego Paul Huaraca Shagñay           ##########
###################################################################

# Parametros de entrada al programa
# S: subyacente
# X: strike
# r: tipo de interes libre de riesgo
# sigma: volatilidad
# Smax: valor maximo que puede tomar el subyacente
# M: numero de incrementos en el subyacente
# N: numero de incrementos en el tiempo

# Metodo Implicito
put_imp <- function(S,X,r,sigma,Smax,M,N){
      dt <- 1
      ds <- Smax/M
      # Estimacion parametros
      a <- numeric(M+1)
      b <- numeric(M+1)
      c <- numeric(M+1)
      for(j in 2:M){
            a[j] <- 0.5*r*(j-1)*dt -0.5*sigma^2*(j-1)^2*dt
            b[j] <- 1 + sigma^2*(j-1)^2*dt + r*dt
            c[j] <- -0.5*r*(j-1)*dt -0.5*sigma^2*(j-1)^2*dt
      }
      # Condiciones
      f <- matrix(0, ncol=N+1, nrow = M+1)
      # Valor intrinseco
      for(j in 1:(M+1)){
            f[j,N+1] <- max(X-(j-1)*ds,0)
      }
      # Si el subyacente es cero
      for(i in 1:(N+1)){
            f[1,i] <- X
      }
      # Si el subyacente alcanza el maximo
      for(i in 1:(N+1)){
            f[M+1,i] <- 0
      }
      # Matriz de coeficientes
      XM <- matrix(0,ncol=M-1, nrow=M-1)
      XM[1,c(1,2)] <- c(b[2],c[2])
      for(i in 2:(M-2)){
            XM[i,c(i-1,i,i+1)] <- c(a[i+1], b[i+1], c[i+1])
      }
      XM[M-1,c(M-2,M-1)] <- c(a[M],b[M])
      # Solucion del enmallado
      for(i in N:1){
            vec <- f[2:M,i+1] - c(a[2]*f[1,i], rep(0,M-3) , c[M]*f[M+1,i])
            f[2:M,i] <- apply(solve(XM)%*%diag(vec), MARGIN = 1, sum)
      }
      return(list(Malla=round(f,4), ValorPut=f[(X/ds)+1,1]))
}

put_imp(50,50,0.02,0.2,100,10,4)


# Metodo explicito
put_exp <- function(S,X,r,sigma,Smax,M,N){
      dt <- 1
      ds <- Smax/M
      # Estimacion parametros
      a <- numeric(M+1)
      b <- numeric(M+1)
      c <- numeric(M+1)
      for(j in 1:(M+1)){
            a[j] <- (-0.5*r*(j-1)*dt + 0.5*sigma^2*(j-1)^2*dt)/(1+r*dt)
            b[j] <- (1 - sigma^2*(j-1)^2*dt)/(1+r*dt)
            c[j] <- (0.5*r*(j-1)*dt + 0.5*sigma^2*(j-1)^2*dt)/(1+r*dt)
      }
      f <- matrix(0, ncol=N+1, nrow = M+1)
      # Valor intrinseco
      for(j in 1:(M+1)){
            f[j,N+1] <- max(X-(j-1)*ds,0)
      }
      # Si el subyacente es cero
      for(i in 1:(N+1)){
            f[1,i] <- X
      }
      # Si el subyacente alcanza el maximo
      for(i in 1:(N+1)){
            f[M+1,i] <- 0
      }
      for(i in N:1){
            for(j in 2:M){
                  f[j,i] <- a[j]*f[j-1,i+1] + b[j]*f[j,i+1] + c[j]*f[j+1,i+1]
            }
      }
      return(list(Malla=round(f,4), ValorPut=f[(X/ds)+1,1]))
}

put_exp(50,50,0.02,0.2,100,10,4)
