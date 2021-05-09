## =======================================================================
##  ******* load exchange module for multi-box model (anybox)  *******  ##
##  updated by YQ- on 30/10/2020                                        ##
##                                                                      ##
## =======================================================================
exclock <- function(u, w, U, W, Carr, cBgarr,h,l,
                    nbox, nbv, nbh,
                    deltaTime, NoS, Numspe)
{
  # this programme uses the 4th order Rugge-Kutta to solve the matrix ODEs
  # Carr: solution in matrix
  # output: Carr[kbox,ibox,Species]
  # e.g., Carr[1,1,1] refers to concentrations of 1st species(NO) in Box[1,1]
  # fw: vertical fluxes
  # fu: horizontal fluxes
  # residual: should be 0 for a correct solution

  k = 1; i = 1;NoS = 1
  k1 = array(0,c(nbh,nbv)); k2 = array(0,c(nbh,nbv));
  k3 = array(0,c(nbh,nbv)); k4 = array(0,c(nbh,nbv));
  
  for (NoS in 1:Numspe) {

  ## -----------------  Runge-Kutta  ------------------ ##
 
  # 1st order EQ 
 
  # Equation for corner position
  k1[1,1] = -(0.5*(sign(W[2,1])+1)*W[2,1]*Carr[1,1,NoS]+0.5*(1-sign(W[2,1]))*W[2,1]*Carr[2,1,NoS])/h[1]-(w[2,1]*(Carr[1,1,NoS]-Carr[2,1,NoS]))/h[1]-(0.5*(sign(U[1,2])+1)*U[1,2]*Carr[1,1,NoS]+0.5*(1-sign(U[1,2]))*U[1,2]*Carr[1,2,NoS])/l[1]-(u[1,2]*(Carr[1,1,NoS]-Carr[1,2,NoS]))/l[1]
  k1[nbv,1] = (0.5*(sign(W[nbv,1])+1)*W[nbv,1]*Carr[nbv-1,1,NoS]+0.5*(1-sign(W[nbv,1]))*W[nbv,1]*Carr[nbv,1,NoS])/h[nbv]+(w[nbv,1]*(Carr[nbv-1,1,NoS]-Carr[nbv,1,NoS]))/h[nbv]-(w[nbv+1,1]*(Carr[nbv,1,NoS]-cBgarr[NoS]))/h[nbv]-(0.5*(sign(U[nbv,2])+1)*U[nbv,2]*Carr[nbv,1,NoS]+0.5*(1-sign(U[nbv,2]))*U[nbv,2]*Carr[nbv,2,NoS])/l[1]-(u[nbv,2]*(Carr[nbv,1,NoS]-Carr[nbv,2,NoS]))/l[1]
  
  k1[1,nbh] = -(0.5*(sign(W[2,nbh])+1)*W[2,nbh]*Carr[1,nbh,NoS]+0.5*(1-sign(W[2,nbh]))*W[2,nbh]*Carr[2,nbh,NoS])/h[1]-(w[2,nbh]*(Carr[1,nbh,NoS]-Carr[2,nbh,NoS]))/h[1]+(0.5*(sign(U[1,nbh])+1)*U[1,nbh]*Carr[1,nbh-1,NoS]+0.5*(1-sign(U[1,nbh]))*U[1,nbh]*Carr[1,nbh,NoS])/l[nbh]+(u[1,nbh]*(Carr[1,nbh-1,NoS]-Carr[1,nbh,NoS]))/l[nbh]
  k1[nbv,nbh] = (0.5*(sign(W[nbv,nbh])+1)*W[nbv,nbh]*Carr[nbv-1,nbh,NoS]+0.5*(1-sign(W[nbv,nbh]))*W[nbv,nbh]*Carr[nbv,nbh,NoS])/h[nbv]+(w[nbv,nbh]*(Carr[nbv-1,nbh,NoS]-Carr[nbv,nbh,NoS]))/h[nbv]-(w[nbv+1,nbh]*(Carr[nbv,nbh,NoS]-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,nbh])+1)*U[nbv,nbh]*Carr[nbv,nbh-1,NoS]+0.5*(1-sign(U[nbv,nbh]))*U[nbv,nbh]*Carr[nbv,nbh,NoS])/l[nbh]+(u[nbv,nbh]*(Carr[nbv,nbh-1,NoS]-Carr[nbv,nbh,NoS]))/l[nbh]  
  
  # Equation for edge i.e. ground or walls.. 
  if (nbv>2)
    # for wall
    for (k in 2:(nbv-1))
    {
      k1[k,1] = (0.5*(sign(W[k,1])+1)*W[k,1]*Carr[k-1,1,NoS]+0.5*(1-sign(W[k,1]))*W[k,1]*Carr[k,1,NoS])/h[k]+(w[k,1]*(Carr[k-1,1,NoS]-Carr[k,1,NoS]))/h[k]-(0.5*(sign(W[k+1,1])+1)*W[k+1,1]*Carr[k,1,NoS]+0.5*(1-sign(W[k+1,1]))*W[k+1,1]*Carr[k+1,1,NoS])/h[k]-(w[k+1,1]*(Carr[k,1,NoS]-Carr[k+1,1,NoS]))/h[k]-(0.5*(sign(U[k,2])+1)*U[k,2]*Carr[k,1,NoS]+0.5*(1-sign(U[k,2]))*U[k,2]*Carr[k,2,NoS])/l[1]-(u[k,2]*(Carr[k,1,NoS]-Carr[k,2,NoS]))/l[1]
      k1[k,nbh] = (0.5*(sign(W[k,nbh])+1)*W[k,nbh]*Carr[k-1,nbh,NoS]+0.5*(1-sign(W[k,nbh]))*W[k,nbh]*Carr[k,nbh,NoS])/h[2]+(w[k,nbh]*(Carr[k-1,nbh,NoS]-Carr[k,nbh,NoS]))/h[k]-(0.5*(sign(W[k+1,nbh])+1)*W[k+1,nbh]*Carr[k,nbh,NoS]+0.5*(1-sign(W[k+1,nbh]))*W[k+1,nbh]*Carr[k+1,nbh,NoS])/h[k]-(w[k+1,nbh]*(Carr[k,nbh,NoS]-Carr[k+1,nbh,NoS]))/h[k]+(0.5*(sign(U[k,nbh])+1)*U[k,nbh]*Carr[k,nbh-1,NoS]+0.5*(1-sign(U[k,nbh]))*U[k,nbh]*Carr[k,nbh,NoS])/l[nbh]+(u[k,nbh]*(Carr[k,nbh-1,NoS]-Carr[k,nbh,NoS]))/l[nbh] 
    }
  if (nbh>2)  
    # for ground and rooftop
    for (i in 2:(nbh-1))
    {
      k1[1,i] = -(0.5*(sign(W[2,i])+1)*W[2,i]*Carr[1,i,NoS]+0.5*(1-sign(W[2,i]))*W[2,i]*Carr[2,i,NoS])/h[1]-(w[2,i]*(Carr[1,i,NoS]-Carr[2,i,NoS]))/h[1]+(0.5*(sign(U[1,i])+1)*U[1,i]*Carr[1,i-1,NoS]+0.5*(1-sign(U[1,i]))*U[1,i]*Carr[1,i,NoS])/l[i]+(u[1,i]*(Carr[1,i-1,NoS]-Carr[1,i,NoS]))/l[i]-(0.5*(sign(U[1,i+1])+1)*U[1,i+1]*Carr[1,i,NoS]+0.5*(1-sign(U[1,i+1]))*U[1,i+1]*Carr[1,i+1,NoS])/l[i]-(u[1,i+1]*(Carr[1,i,NoS]-Carr[1,i+1,NoS]))/l[i]
      k1[nbv,i] = (0.5*(sign(W[nbv,i])+1)*W[nbv,i]*Carr[nbv-1,i,NoS]+0.5*(1-sign(W[nbv,i]))*W[nbv,i]*Carr[nbv,i,NoS])/h[nbv]+(w[nbv,i]*(Carr[nbv-1,i,NoS]-Carr[nbv,i,NoS]))/h[nbv]-(w[nbv+1,i]*(Carr[nbv,i,NoS]-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,i])+1)*U[nbv,i]*Carr[nbv,i-1,NoS]+0.5*(1-sign(U[nbv,i]))*U[nbv,i]*Carr[nbv,i,NoS])/l[i]+(u[nbv,i]*(Carr[nbv,i-1,NoS]-Carr[nbv,i,NoS]))/l[i]-(0.5*(sign(U[nbv,i+1])+1)*U[nbv,i+1]*Carr[nbv,i,NoS]+0.5*(1-sign(U[nbv,i+1]))*U[nbv,i+1]*Carr[nbv,i+1,NoS])/l[i]-(u[nbv,i+1]*(Carr[nbv,i,NoS]-Carr[nbv,i+1,NoS]))/l[i]
    }
  #
  if (nbv>2 & nbh>2)
    # for in-canyon box
    for (k in 2:(nbv-1))
    {
      for (i in 2:(nbh-1))
      {
        k1[k,i] = (0.5*(sign(W[k,i])+1)*W[k,i]*Carr[k-1,i,NoS]+0.5*(1-sign(W[k,i]))*W[k,i]*Carr[k,i,NoS])/h[k]+(w[k,i]*(Carr[k-1,i,NoS]-Carr[k,i,NoS]))/h[k]-(0.5*(sign(W[k+1,i])+1)*W[k+1,i]*Carr[k,i,NoS]+0.5*(1-sign(W[k+1,i]))*W[k+1,i]*Carr[k+1,i,NoS])/h[k]-(w[k+1,i]*(Carr[k,i,NoS]-Carr[k+1,i,NoS]))/h[k]+(0.5*(sign(U[k,i])+1)*U[k,i]*Carr[k,i-1,NoS]+0.5*(1-sign(U[k,i]))*U[k,i]*Carr[k,i,NoS])/l[i]+(u[k,i]*(Carr[k,i-1,NoS]-Carr[k,i,NoS]))/l[i]-(0.5*(sign(U[k,i+1])+1)*U[k,i+1]*Carr[k,i,NoS]+0.5*(1-sign(U[k,i+1]))*U[k,i+1]*Carr[k,i+1,NoS])/l[i]-(u[k,i+1]*(Carr[k,i,NoS]-Carr[k,i+1,NoS]))/l[i]
      }
    }
 
  # 2nd order EQ 
 
  # Equation for corner position
  k2[1,1] = -(0.5*(sign(W[2,1])+1)*W[2,1]*(Carr[1,1,NoS]+deltaTime*k1[1,1]/2)+0.5*(1-sign(W[2,1]))*W[2,1]*Carr[2,1,NoS])/h[1]-(w[2,1]*((Carr[1,1,NoS]+deltaTime*k1[1,1]/2)-Carr[2,1,NoS]))/h[1]-(0.5*(sign(U[1,2])+1)*U[1,2]*(Carr[1,1,NoS]+deltaTime*k1[1,1]/2)+0.5*(1-sign(U[1,2]))*U[1,2]*Carr[1,2,NoS])/l[1]-(u[1,2]*((Carr[1,1,NoS]+deltaTime*k1[1,1]/2)-Carr[1,2,NoS]))/l[1]
  k2[nbv,1] = (0.5*(sign(W[nbv,1])+1)*W[nbv,1]*Carr[nbv-1,1,NoS]+0.5*(1-sign(W[nbv,1]))*W[nbv,1]*(Carr[nbv,1,NoS]+deltaTime*k1[nbv,1]/2))/h[nbv]+(w[nbv,1]*(Carr[nbv-1,1,NoS]-(Carr[nbv,1,NoS]+deltaTime*k1[nbv,1]/2)))/h[nbv]-(w[nbv+1,1]*((Carr[nbv,1,NoS]+deltaTime*k1[nbv,1]/2)-cBgarr[NoS]))/h[nbv]-(0.5*(sign(U[nbv,2])+1)*U[nbv,2]*(Carr[nbv,1,NoS]+deltaTime*k1[nbv,1]/2)+0.5*(1-sign(U[nbv,2]))*U[nbv,2]*Carr[nbv,2,NoS])/l[1]-(u[nbv,2]*((Carr[nbv,1,NoS]+deltaTime*k1[nbv,1]/2)-Carr[nbv,2,NoS]))/l[1]
  
  k2[1,nbh] = -(0.5*(sign(W[2,nbh])+1)*W[2,nbh]*(Carr[1,nbh,NoS]+deltaTime*k1[1,nbh]/2)+0.5*(1-sign(W[2,nbh]))*W[2,nbh]*Carr[2,nbh,NoS])/h[1]-(w[2,nbh]*((Carr[1,nbh,NoS]+deltaTime*k1[1,nbh]/2)-Carr[2,nbh,NoS]))/h[1]+(0.5*(sign(U[1,nbh])+1)*U[1,nbh]*Carr[1,nbh-1,NoS]+0.5*(1-sign(U[1,nbh]))*U[1,nbh]*(Carr[1,nbh,NoS]+deltaTime*k1[1,nbh]/2))/l[nbh]+(u[1,nbh]*(Carr[1,nbh-1,NoS]-(Carr[1,nbh,NoS]+deltaTime*k1[1,nbh]/2)))/l[nbh]
  k2[nbv,nbh] = (0.5*(sign(W[nbv,nbh])+1)*W[nbv,nbh]*Carr[nbv-1,nbh,NoS]+0.5*(1-sign(W[nbv,nbh]))*W[nbv,nbh]*(Carr[nbv,nbh,NoS]+deltaTime*k1[nbv,nbh]/2))/h[nbv]+(w[nbv,nbh]*(Carr[nbv-1,nbh,NoS]-(Carr[nbv,nbh,NoS]+deltaTime*k1[nbv,nbh]/2)))/h[nbv]-(w[nbv+1,nbh]*((Carr[nbv,nbh,NoS]+deltaTime*k1[nbv,nbh]/2)-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,nbh])+1)*U[nbv,nbh]*Carr[nbv,nbh-1,NoS]+0.5*(1-sign(U[nbv,nbh]))*U[nbv,nbh]*(Carr[nbv,nbh,NoS]+deltaTime*k1[nbv,nbh]/2))/l[nbh]+(u[nbv,nbh]*(Carr[nbv,nbh-1,NoS]-(Carr[nbv,nbh,NoS]+deltaTime*k1[nbv,nbh]/2)))/l[nbh]  
  
  # Equation for edge i.e. ground or walls...
  if (nbv>2)
    # for wall
    for (k in 2:(nbv-1))
    {
      k2[k,1] = (0.5*(sign(W[k,1])+1)*W[k,1]*Carr[k-1,1,NoS]+0.5*(1-sign(W[k,1]))*W[k,1]*(Carr[k,1,NoS]+deltaTime*k1[k,1]/2))/h[k]+(w[k,1]*(Carr[k-1,1,NoS]-(Carr[k,1,NoS]+deltaTime*k1[k,1]/2)))/h[k]-(0.5*(sign(W[k+1,1])+1)*W[k+1,1]*(Carr[k,1,NoS]+deltaTime*k1[k,1]/2)+0.5*(1-sign(W[k+1,1]))*W[k+1,1]*Carr[k+1,1,NoS])/h[k]-(w[k+1,1]*((Carr[k,1,NoS]+deltaTime*k1[k,1]/2)-Carr[k+1,1,NoS]))/h[k]-(0.5*(sign(U[k,2])+1)*U[k,2]*(Carr[k,1,NoS]+deltaTime*k1[k,1]/2)+0.5*(1-sign(U[k,2]))*U[k,2]*Carr[k,2,NoS])/l[1]-(u[k,2]*((Carr[k,1,NoS]+deltaTime*k1[k,1]/2)-Carr[k,2,NoS]))/l[1]
      k2[k,nbh] = (0.5*(sign(W[k,nbh])+1)*W[k,nbh]*Carr[k-1,nbh,NoS]+0.5*(1-sign(W[k,nbh]))*W[k,nbh]*(Carr[k,nbh,NoS]+deltaTime*k1[k,nbh]/2))/h[2]+(w[k,nbh]*(Carr[k-1,nbh,NoS]-(Carr[k,nbh,NoS]+deltaTime*k1[k,nbh]/2)))/h[k]-(0.5*(sign(W[k+1,nbh])+1)*W[k+1,nbh]*(Carr[k,nbh,NoS]+deltaTime*k1[k,nbh]/2)+0.5*(1-sign(W[k+1,nbh]))*W[k+1,nbh]*Carr[k+1,nbh,NoS])/h[k]-(w[k+1,nbh]*((Carr[k,nbh,NoS]+deltaTime*k1[k,nbh]/2)-Carr[k+1,nbh,NoS]))/h[k]+(0.5*(sign(U[k,nbh])+1)*U[k,nbh]*Carr[k,nbh-1,NoS]+0.5*(1-sign(U[k,nbh]))*U[k,nbh]*(Carr[k,nbh,NoS]+deltaTime*k1[k,nbh]/2))/l[nbh]+(u[k,nbh]*(Carr[k,nbh-1,NoS]-(Carr[k,nbh,NoS]+deltaTime*k1[k,nbh]/2)))/l[nbh] 
    }
  if (nbh>2)  
    # for ground and rooftop
    for (i in 2:(nbh-1))
    {
      k2[1,i] = -(0.5*(sign(W[2,i])+1)*W[2,i]*(Carr[1,i,NoS]+deltaTime*k1[1,i]/2)+0.5*(1-sign(W[2,i]))*W[2,i]*Carr[2,i,NoS])/h[1]-(w[2,i]*((Carr[1,i,NoS]+deltaTime*k1[1,i]/2)-Carr[2,i,NoS]))/h[1]+(0.5*(sign(U[1,i])+1)*U[1,i]*Carr[1,i-1,NoS]+0.5*(1-sign(U[1,i]))*U[1,i]*(Carr[1,i,NoS]+deltaTime*k1[1,i]/2))/l[i]+(u[1,i]*(Carr[1,i-1,NoS]-(Carr[1,i,NoS]+deltaTime*k1[1,i]/2)))/l[i]-(0.5*(sign(U[1,i+1])+1)*U[1,i+1]*(Carr[1,i,NoS]+deltaTime*k1[1,i]/2)+0.5*(1-sign(U[1,i+1]))*U[1,i+1]*Carr[1,i+1,NoS])/l[i]-(u[1,i+1]*((Carr[1,i,NoS]+deltaTime*k1[1,i]/2)-Carr[1,i+1,NoS]))/l[i]
      k2[nbv,i] = (0.5*(sign(W[nbv,i])+1)*W[nbv,i]*Carr[nbv-1,i,NoS]+0.5*(1-sign(W[nbv,i]))*W[nbv,i]*(Carr[nbv,i,NoS]+deltaTime*k1[nbv,i]/2))/h[nbv]+(w[nbv,i]*(Carr[nbv-1,i,NoS]-(Carr[nbv,i,NoS]+deltaTime*k1[nbv,i]/2)))/h[nbv]-(w[nbv+1,i]*((Carr[nbv,i,NoS]+deltaTime*k1[nbv,i]/2)-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,i])+1)*U[nbv,i]*Carr[nbv,i-1,NoS]+0.5*(1-sign(U[nbv,i]))*U[nbv,i]*(Carr[nbv,i,NoS]+deltaTime*k1[nbv,i]/2))/l[i]+(u[nbv,i]*(Carr[nbv,i-1,NoS]-(Carr[nbv,i,NoS]+deltaTime*k1[nbv,i]/2)))/l[i]-(0.5*(sign(U[nbv,i+1])+1)*U[nbv,i+1]*(Carr[nbv,i,NoS]+deltaTime*k1[nbv,i]/2)+0.5*(1-sign(U[nbv,i+1]))*U[nbv,i+1]*Carr[nbv,i+1,NoS])/l[i]-(u[nbv,i+1]*((Carr[nbv,i,NoS]+deltaTime*k1[nbv,i]/2)-Carr[nbv,i+1,NoS]))/l[i]
    }
  #
  if (nbv>2 & nbh>2)
    # for in-canyon box
    for (k in 2:(nbv-1))
    {
      for (i in 2:(nbh-1))
      {
        k2[k,i] = (0.5*(sign(W[k,i])+1)*W[k,i]*Carr[k-1,i,NoS]+0.5*(1-sign(W[k,i]))*W[k,i]*(Carr[k,i,NoS]+deltaTime*k1[k,i]/2))/h[k]+(w[k,i]*(Carr[k-1,i,NoS]-(Carr[k,i,NoS]+deltaTime*k1[k,i]/2)))/h[k]-(0.5*(sign(W[k+1,i])+1)*W[k+1,i]*(Carr[k,i,NoS]+deltaTime*k1[k,i]/2)+0.5*(1-sign(W[k+1,i]))*W[k+1,i]*Carr[k+1,i,NoS])/h[k]-(w[k+1,i]*((Carr[k,i,NoS]+deltaTime*k1[k,i]/2)-Carr[k+1,i,NoS]))/h[k]+(0.5*(sign(U[k,i])+1)*U[k,i]*Carr[k,i-1,NoS]+0.5*(1-sign(U[k,i]))*U[k,i]*(Carr[k,i,NoS]+deltaTime*k1[k,i]/2))/l[i]+(u[k,i]*(Carr[k,i-1,NoS]-(Carr[k,i,NoS]+deltaTime*k1[k,i]/2)))/l[i]-(0.5*(sign(U[k,i+1])+1)*U[k,i+1]*(Carr[k,i,NoS]+deltaTime*k1[k,i]/2)+0.5*(1-sign(U[k,i+1]))*U[k,i+1]*Carr[k,i+1,NoS])/l[i]-(u[k,i+1]*((Carr[k,i,NoS]+deltaTime*k1[k,i]/2)-Carr[k,i+1,NoS]))/l[i]
      }
    }

  # 3rd order EQ 

  # Equation for corner position
  k3[1,1] = -(0.5*(sign(W[2,1])+1)*W[2,1]*(Carr[1,1,NoS]+deltaTime*k2[1,1]/2)+0.5*(1-sign(W[2,1]))*W[2,1]*Carr[2,1,NoS])/h[1]-(w[2,1]*((Carr[1,1,NoS]+deltaTime*k2[1,1]/2)-Carr[2,1,NoS]))/h[1]-(0.5*(sign(U[1,2])+1)*U[1,2]*(Carr[1,1,NoS]+deltaTime*k2[1,1]/2)+0.5*(1-sign(U[1,2]))*U[1,2]*Carr[1,2,NoS])/l[1]-(u[1,2]*((Carr[1,1,NoS]+deltaTime*k2[1,1]/2)-Carr[1,2,NoS]))/l[1]
  k3[nbv,1] = (0.5*(sign(W[nbv,1])+1)*W[nbv,1]*Carr[nbv-1,1,NoS]+0.5*(1-sign(W[nbv,1]))*W[nbv,1]*(Carr[nbv,1,NoS]+deltaTime*k2[nbv,1]/2))/h[nbv]+(w[nbv,1]*(Carr[nbv-1,1,NoS]-(Carr[nbv,1,NoS]+deltaTime*k2[nbv,1]/2)))/h[nbv]-(w[nbv+1,1]*((Carr[nbv,1,NoS]+deltaTime*k2[nbv,1]/2)-cBgarr[NoS]))/h[nbv]-(0.5*(sign(U[nbv,2])+1)*U[nbv,2]*(Carr[nbv,1,NoS]+deltaTime*k2[nbv,1]/2)+0.5*(1-sign(U[nbv,2]))*U[nbv,2]*Carr[nbv,2,NoS])/l[1]-(u[nbv,2]*((Carr[nbv,1,NoS]+deltaTime*k2[nbv,1]/2)-Carr[nbv,2,NoS]))/l[1]
  
  k3[1,nbh] = -(0.5*(sign(W[2,nbh])+1)*W[2,nbh]*(Carr[1,nbh,NoS]+deltaTime*k2[1,nbh]/2)+0.5*(1-sign(W[2,nbh]))*W[2,nbh]*Carr[2,nbh,NoS])/h[1]-(w[2,nbh]*((Carr[1,nbh,NoS]+deltaTime*k2[1,nbh]/2)-Carr[2,nbh,NoS]))/h[1]+(0.5*(sign(U[1,nbh])+1)*U[1,nbh]*Carr[1,nbh-1,NoS]+0.5*(1-sign(U[1,nbh]))*U[1,nbh]*(Carr[1,nbh,NoS]+deltaTime*k2[1,nbh]/2))/l[nbh]+(u[1,nbh]*(Carr[1,nbh-1,NoS]-(Carr[1,nbh,NoS]+deltaTime*k2[1,nbh]/2)))/l[nbh]
  k3[nbv,nbh] = (0.5*(sign(W[nbv,nbh])+1)*W[nbv,nbh]*Carr[nbv-1,nbh,NoS]+0.5*(1-sign(W[nbv,nbh]))*W[nbv,nbh]*(Carr[nbv,nbh,NoS]+deltaTime*k2[nbv,nbh]/2))/h[nbv]+(w[nbv,nbh]*(Carr[nbv-1,nbh,NoS]-(Carr[nbv,nbh,NoS]+deltaTime*k2[nbv,nbh]/2)))/h[nbv]-(w[nbv+1,nbh]*((Carr[nbv,nbh,NoS]+deltaTime*k2[nbv,nbh]/2)-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,nbh])+1)*U[nbv,nbh]*Carr[nbv,nbh-1,NoS]+0.5*(1-sign(U[nbv,nbh]))*U[nbv,nbh]*(Carr[nbv,nbh,NoS]+deltaTime*k2[nbv,nbh]/2))/l[nbh]+(u[nbv,nbh]*(Carr[nbv,nbh-1,NoS]-(Carr[nbv,nbh,NoS]+deltaTime*k2[nbv,nbh]/2)))/l[nbh]  
  
  # Equation for edge i.e. ground or walls...
  if (nbv>2)
    # for wall
    for (k in 2:(nbv-1))
    {
      k3[k,1] = (0.5*(sign(W[k,1])+1)*W[k,1]*Carr[k-1,1,NoS]+0.5*(1-sign(W[k,1]))*W[k,1]*(Carr[k,1,NoS]+deltaTime*k2[k,1]/2))/h[k]+(w[k,1]*(Carr[k-1,1,NoS]-(Carr[k,1,NoS]+deltaTime*k2[k,1]/2)))/h[k]-(0.5*(sign(W[k+1,1])+1)*W[k+1,1]*(Carr[k,1,NoS]+deltaTime*k2[k,1]/2)+0.5*(1-sign(W[k+1,1]))*W[k+1,1]*Carr[k+1,1,NoS])/h[k]-(w[k+1,1]*((Carr[k,1,NoS]+deltaTime*k2[k,1]/2)-Carr[k+1,1,NoS]))/h[k]-(0.5*(sign(U[k,2])+1)*U[k,2]*(Carr[k,1,NoS]+deltaTime*k2[k,1]/2)+0.5*(1-sign(U[k,2]))*U[k,2]*Carr[k,2,NoS])/l[1]-(u[k,2]*((Carr[k,1,NoS]+deltaTime*k2[k,1]/2)-Carr[k,2,NoS]))/l[1]
      k3[k,nbh] = (0.5*(sign(W[k,nbh])+1)*W[k,nbh]*Carr[k-1,nbh,NoS]+0.5*(1-sign(W[k,nbh]))*W[k,nbh]*(Carr[k,nbh,NoS]+deltaTime*k2[k,nbh]/2))/h[2]+(w[k,nbh]*(Carr[k-1,nbh,NoS]-(Carr[k,nbh,NoS]+deltaTime*k2[k,nbh]/2)))/h[k]-(0.5*(sign(W[k+1,nbh])+1)*W[k+1,nbh]*(Carr[k,nbh,NoS]+deltaTime*k2[k,nbh]/2)+0.5*(1-sign(W[k+1,nbh]))*W[k+1,nbh]*Carr[k+1,nbh,NoS])/h[k]-(w[k+1,nbh]*((Carr[k,nbh,NoS]+deltaTime*k2[k,nbh]/2)-Carr[k+1,nbh,NoS]))/h[k]+(0.5*(sign(U[k,nbh])+1)*U[k,nbh]*Carr[k,nbh-1,NoS]+0.5*(1-sign(U[k,nbh]))*U[k,nbh]*(Carr[k,nbh,NoS]+deltaTime*k2[k,nbh]/2))/l[nbh]+(u[k,nbh]*(Carr[k,nbh-1,NoS]-(Carr[k,nbh,NoS]+deltaTime*k2[k,nbh]/2)))/l[nbh] 
    }
  if (nbh>2)  
    # for ground and rooftop
    for (i in 2:(nbh-1))
    {
      k3[1,i] = -(0.5*(sign(W[2,i])+1)*W[2,i]*(Carr[1,i,NoS]+deltaTime*k2[1,i]/2)+0.5*(1-sign(W[2,i]))*W[2,i]*Carr[2,i,NoS])/h[1]-(w[2,i]*((Carr[1,i,NoS]+deltaTime*k2[1,i]/2)-Carr[2,i,NoS]))/h[1]+(0.5*(sign(U[1,i])+1)*U[1,i]*Carr[1,i-1,NoS]+0.5*(1-sign(U[1,i]))*U[1,i]*(Carr[1,i,NoS]+deltaTime*k2[1,i]/2))/l[i]+(u[1,i]*(Carr[1,i-1,NoS]-(Carr[1,i,NoS]+deltaTime*k2[1,i]/2)))/l[i]-(0.5*(sign(U[1,i+1])+1)*U[1,i+1]*(Carr[1,i,NoS]+deltaTime*k2[1,i]/2)+0.5*(1-sign(U[1,i+1]))*U[1,i+1]*Carr[1,i+1,NoS])/l[i]-(u[1,i+1]*((Carr[1,i,NoS]+deltaTime*k2[1,i]/2)-Carr[1,i+1,NoS]))/l[i]
      k3[nbv,i] = (0.5*(sign(W[nbv,i])+1)*W[nbv,i]*Carr[nbv-1,i,NoS]+0.5*(1-sign(W[nbv,i]))*W[nbv,i]*(Carr[nbv,i,NoS]+deltaTime*k2[nbv,i]/2))/h[nbv]+(w[nbv,i]*(Carr[nbv-1,i,NoS]-(Carr[nbv,i,NoS]+deltaTime*k2[nbv,i]/2)))/h[nbv]-(w[nbv+1,i]*((Carr[nbv,i,NoS]+deltaTime*k2[nbv,i]/2)-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,i])+1)*U[nbv,i]*Carr[nbv,i-1,NoS]+0.5*(1-sign(U[nbv,i]))*U[nbv,i]*(Carr[nbv,i,NoS]+deltaTime*k2[nbv,i]/2))/l[i]+(u[nbv,i]*(Carr[nbv,i-1,NoS]-(Carr[nbv,i,NoS]+deltaTime*k2[nbv,i]/2)))/l[i]-(0.5*(sign(U[nbv,i+1])+1)*U[nbv,i+1]*(Carr[nbv,i,NoS]+deltaTime*k2[nbv,i]/2)+0.5*(1-sign(U[nbv,i+1]))*U[nbv,i+1]*Carr[nbv,i+1,NoS])/l[i]-(u[nbv,i+1]*((Carr[nbv,i,NoS]+deltaTime*k2[nbv,i]/2)-Carr[nbv,i+1,NoS]))/l[i]
    }
  #
  if (nbv>2 & nbh>2)
    ## Equation for in-canyon box
    for (k in 2:(nbv-1))
    {
      for (i in 2:(nbh-1))
      {
        k3[k,i] = (0.5*(sign(W[k,i])+1)*W[k,i]*Carr[k-1,i,NoS]+0.5*(1-sign(W[k,i]))*W[k,i]*(Carr[k,i,NoS]+deltaTime*k2[k,i]/2))/h[k]+(w[k,i]*(Carr[k-1,i,NoS]-(Carr[k,i,NoS]+deltaTime*k2[k,i]/2)))/h[k]-(0.5*(sign(W[k+1,i])+1)*W[k+1,i]*(Carr[k,i,NoS]+deltaTime*k2[k,i]/2)+0.5*(1-sign(W[k+1,i]))*W[k+1,i]*Carr[k+1,i,NoS])/h[k]-(w[k+1,i]*((Carr[k,i,NoS]+deltaTime*k2[k,i]/2)-Carr[k+1,i,NoS]))/h[k]+(0.5*(sign(U[k,i])+1)*U[k,i]*Carr[k,i-1,NoS]+0.5*(1-sign(U[k,i]))*U[k,i]*(Carr[k,i,NoS]+deltaTime*k2[k,i]/2))/l[i]+(u[k,i]*(Carr[k,i-1,NoS]-(Carr[k,i,NoS]+deltaTime*k2[k,i]/2)))/l[i]-(0.5*(sign(U[k,i+1])+1)*U[k,i+1]*(Carr[k,i,NoS]+deltaTime*k2[k,i]/2)+0.5*(1-sign(U[k,i+1]))*U[k,i+1]*Carr[k,i+1,NoS])/l[i]-(u[k,i+1]*((Carr[k,i,NoS]+deltaTime*k2[k,i]/2)-Carr[k,i+1,NoS]))/l[i]
      }
    }
  
  # 4th order EQ 

  # Equation for corner position
  k4[1,1] = -(0.5*(sign(W[2,1])+1)*W[2,1]*(Carr[1,1,NoS]+deltaTime*k3[1,1])+0.5*(1-sign(W[2,1]))*W[2,1]*Carr[2,1,NoS])/h[1]-(w[2,1]*((Carr[1,1,NoS]+deltaTime*k3[1,1])-Carr[2,1,NoS]))/h[1]-(0.5*(sign(U[1,2])+1)*U[1,2]*(Carr[1,1,NoS]+deltaTime*k3[1,1])+0.5*(1-sign(U[1,2]))*U[1,2]*Carr[1,2,NoS])/l[1]-(u[1,2]*((Carr[1,1,NoS]+deltaTime*k3[1,1])-Carr[1,2,NoS]))/l[1]
  k4[nbv,1] = (0.5*(sign(W[nbv,1])+1)*W[nbv,1]*Carr[nbv-1,1,NoS]+0.5*(1-sign(W[nbv,1]))*W[nbv,1]*(Carr[nbv,1,NoS]+deltaTime*k3[nbv,1]))/h[nbv]+(w[nbv,1]*(Carr[nbv-1,1,NoS]-(Carr[nbv,1,NoS]+deltaTime*k3[nbv,1])))/h[nbv]-(w[nbv+1,1]*((Carr[nbv,1,NoS]+deltaTime*k3[nbv,1])-cBgarr[NoS]))/h[nbv]-(0.5*(sign(U[nbv,2])+1)*U[nbv,2]*(Carr[nbv,1,NoS]+deltaTime*k3[nbv,1])+0.5*(1-sign(U[nbv,2]))*U[nbv,2]*Carr[nbv,2,NoS])/l[1]-(u[nbv,2]*((Carr[nbv,1,NoS]+deltaTime*k3[nbv,1])-Carr[nbv,2,NoS]))/l[1]
  
  k4[1,nbh] = -(0.5*(sign(W[2,nbh])+1)*W[2,nbh]*(Carr[1,nbh,NoS]+deltaTime*k3[1,nbh])+0.5*(1-sign(W[2,nbh]))*W[2,nbh]*Carr[2,nbh,NoS])/h[1]-(w[2,nbh]*((Carr[1,nbh,NoS]+deltaTime*k3[1,nbh])-Carr[2,nbh,NoS]))/h[1]+(0.5*(sign(U[1,nbh])+1)*U[1,nbh]*Carr[1,nbh-1,NoS]+0.5*(1-sign(U[1,nbh]))*U[1,nbh]*(Carr[1,nbh,NoS]+deltaTime*k3[1,nbh]))/l[nbh]+(u[1,nbh]*(Carr[1,nbh-1,NoS]-(Carr[1,nbh,NoS]+deltaTime*k3[1,nbh])))/l[nbh]
  k4[nbv,nbh] = (0.5*(sign(W[nbv,nbh])+1)*W[nbv,nbh]*Carr[nbv-1,nbh,NoS]+0.5*(1-sign(W[nbv,nbh]))*W[nbv,nbh]*(Carr[nbv,nbh,NoS]+deltaTime*k3[nbv,nbh]))/h[nbv]+(w[nbv,nbh]*(Carr[nbv-1,nbh,NoS]-(Carr[nbv,nbh,NoS]+deltaTime*k3[nbv,nbh])))/h[nbv]-(w[nbv+1,nbh]*((Carr[nbv,nbh,NoS]+deltaTime*k3[nbv,nbh])-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,nbh])+1)*U[nbv,nbh]*Carr[nbv,nbh-1,NoS]+0.5*(1-sign(U[nbv,nbh]))*U[nbv,nbh]*(Carr[nbv,nbh,NoS]+deltaTime*k3[nbv,nbh]))/l[nbh]+(u[nbv,nbh]*(Carr[nbv,nbh-1,NoS]-(Carr[nbv,nbh,NoS]+deltaTime*k3[nbv,nbh])))/l[nbh]  
  
  # Equation for edge i.e. ground or walls...
  if (nbv>2)
    # for wall
    for (k in 2:(nbv-1))
    {
      k4[k,1] = (0.5*(sign(W[k,1])+1)*W[k,1]*Carr[k-1,1,NoS]+0.5*(1-sign(W[k,1]))*W[k,1]*(Carr[k,1,NoS]+deltaTime*k3[k,1]))/h[k]+(w[k,1]*(Carr[k-1,1,NoS]-(Carr[k,1,NoS]+deltaTime*k3[k,1])))/h[k]-(0.5*(sign(W[k+1,1])+1)*W[k+1,1]*(Carr[k,1,NoS]+deltaTime*k3[k,1])+0.5*(1-sign(W[k+1,1]))*W[k+1,1]*Carr[k+1,1,NoS])/h[k]-(w[k+1,1]*((Carr[k,1,NoS]+deltaTime*k3[k,1])-Carr[k+1,1,NoS]))/h[k]-(0.5*(sign(U[k,2])+1)*U[k,2]*(Carr[k,1,NoS]+deltaTime*k3[k,1])+0.5*(1-sign(U[k,2]))*U[k,2]*Carr[k,2,NoS])/l[1]-(u[k,2]*((Carr[k,1,NoS]+deltaTime*k3[k,1])-Carr[k,2,NoS]))/l[1]
      k4[k,nbh] = (0.5*(sign(W[k,nbh])+1)*W[k,nbh]*Carr[k-1,nbh,NoS]+0.5*(1-sign(W[k,nbh]))*W[k,nbh]*(Carr[k,nbh,NoS]+deltaTime*k3[k,nbh]))/h[2]+(w[k,nbh]*(Carr[k-1,nbh,NoS]-(Carr[k,nbh,NoS]+deltaTime*k3[k,nbh])))/h[k]-(0.5*(sign(W[k+1,nbh])+1)*W[k+1,nbh]*(Carr[k,nbh,NoS]+deltaTime*k3[k,nbh])+0.5*(1-sign(W[k+1,nbh]))*W[k+1,nbh]*Carr[k+1,nbh,NoS])/h[k]-(w[k+1,nbh]*((Carr[k,nbh,NoS]+deltaTime*k3[k,nbh])-Carr[k+1,nbh,NoS]))/h[k]+(0.5*(sign(U[k,nbh])+1)*U[k,nbh]*Carr[k,nbh-1,NoS]+0.5*(1-sign(U[k,nbh]))*U[k,nbh]*(Carr[k,nbh,NoS]+deltaTime*k3[k,nbh]))/l[nbh]+(u[k,nbh]*(Carr[k,nbh-1,NoS]-(Carr[k,nbh,NoS]+deltaTime*k3[k,nbh])))/l[nbh] 
    }
  if (nbh>2)  
    # for ground and rooftop
    for (i in 2:(nbh-1))
    {
      k4[1,i] = -(0.5*(sign(W[2,i])+1)*W[2,i]*(Carr[1,i,NoS]+deltaTime*k3[1,i])+0.5*(1-sign(W[2,i]))*W[2,i]*Carr[2,i,NoS])/h[1]-(w[2,i]*((Carr[1,i,NoS]+deltaTime*k3[1,i])-Carr[2,i,NoS]))/h[1]+(0.5*(sign(U[1,i])+1)*U[1,i]*Carr[1,i-1,NoS]+0.5*(1-sign(U[1,i]))*U[1,i]*(Carr[1,i,NoS]+deltaTime*k3[1,i]))/l[i]+(u[1,i]*(Carr[1,i-1,NoS]-(Carr[1,i,NoS]+deltaTime*k3[1,i])))/l[i]-(0.5*(sign(U[1,i+1])+1)*U[1,i+1]*(Carr[1,i,NoS]+deltaTime*k3[1,i])+0.5*(1-sign(U[1,i+1]))*U[1,i+1]*Carr[1,i+1,NoS])/l[i]-(u[1,i+1]*((Carr[1,i,NoS]+deltaTime*k3[1,i])-Carr[1,i+1,NoS]))/l[i]
      k4[nbv,i] = (0.5*(sign(W[nbv,i])+1)*W[nbv,i]*Carr[nbv-1,i,NoS]+0.5*(1-sign(W[nbv,i]))*W[nbv,i]*(Carr[nbv,i,NoS]+deltaTime*k3[nbv,i]))/h[nbv]+(w[nbv,i]*(Carr[nbv-1,i,NoS]-(Carr[nbv,i,NoS]+deltaTime*k3[nbv,i])))/h[nbv]-(w[nbv+1,i]*((Carr[nbv,i,NoS]+deltaTime*k3[nbv,i])-cBgarr[NoS]))/h[nbv]+(0.5*(sign(U[nbv,i])+1)*U[nbv,i]*Carr[nbv,i-1,NoS]+0.5*(1-sign(U[nbv,i]))*U[nbv,i]*(Carr[nbv,i,NoS]+deltaTime*k3[nbv,i]))/l[i]+(u[nbv,i]*(Carr[nbv,i-1,NoS]-(Carr[nbv,i,NoS]+deltaTime*k3[nbv,i])))/l[i]-(0.5*(sign(U[nbv,i+1])+1)*U[nbv,i+1]*(Carr[nbv,i,NoS]+deltaTime*k3[nbv,i])+0.5*(1-sign(U[nbv,i+1]))*U[nbv,i+1]*Carr[nbv,i+1,NoS])/l[i]-(u[nbv,i+1]*((Carr[nbv,i,NoS]+deltaTime*k3[nbv,i])-Carr[nbv,i+1,NoS]))/l[i]
    }
  
  if (nbv>2 & nbh>2)
    # Equation for in-canyon box
    for (k in 2:(nbv-1))
    {
      for (i in 2:(nbh-1))
      {
        k4[k,i] = (0.5*(sign(W[k,i])+1)*W[k,i]*Carr[k-1,i,NoS]+0.5*(1-sign(W[k,i]))*W[k,i]*(Carr[k,i,NoS]+deltaTime*k3[k,i]))/h[k]+(w[k,i]*(Carr[k-1,i,NoS]-(Carr[k,i,NoS]+deltaTime*k3[k,i])))/h[k]-(0.5*(sign(W[k+1,i])+1)*W[k+1,i]*(Carr[k,i,NoS]+deltaTime*k3[k,i])+0.5*(1-sign(W[k+1,i]))*W[k+1,i]*Carr[k+1,i,NoS])/h[k]-(w[k+1,i]*((Carr[k,i,NoS]+deltaTime*k3[k,i])-Carr[k+1,i,NoS]))/h[k]+(0.5*(sign(U[k,i])+1)*U[k,i]*Carr[k,i-1,NoS]+0.5*(1-sign(U[k,i]))*U[k,i]*(Carr[k,i,NoS]+deltaTime*k3[k,i]))/l[i]+(u[k,i]*(Carr[k,i-1,NoS]-(Carr[k,i,NoS]+deltaTime*k3[k,i])))/l[i]-(0.5*(sign(U[k,i+1])+1)*U[k,i+1]*(Carr[k,i,NoS]+deltaTime*k3[k,i])+0.5*(1-sign(U[k,i+1]))*U[k,i+1]*Carr[k,i+1,NoS])/l[i]-(u[k,i+1]*((Carr[k,i,NoS]+deltaTime*k3[k,i])-Carr[k,i+1,NoS]))/l[i]  
      }
    }
  
  # output
  for (k in 1:nbv){
    for (i in 1:nbh){
      Carr[k,i,NoS] = Carr[k,i,NoS] + deltaTime*(k1[k,i]+2*k2[k,i]+2*k3[k,i]+k4[k,i])/6
    }
  }
  }
  return(Carr)
}

## ----------- flux calculation -YQ on 29/10/2020 -------------- ##
# Note
# Ft[k,i] = (u[k,i](Carr[k,i-1,NoS]-Carr[k,i,NoS]))
# Fa[k,i] = (0.5*(sign(U[k,i])+1)*U[k,i]*Carr[k,i-1,NoS]+0.5*(1-sign(U[k,i]))*U[k,i]*Carr[k,i,NoS])
# Gt[k,i] = (w[k,i](Carr[k-1,i,NoS]-Carr[k,i,NoS]))
# Ga[k,i] = (0.5*(sign(W[k,i])+1)*W[k,i]*Carr[k-1,i,NoS]+0.5*(1-sign(W[k,i]))*W[k,i]*Carr[k,i,NoS])

# Ft = array(0,c(nbh,nbv)); Gt = array(0,c(nbh+1,nbv));
# Fa = array(0,c(nbh,nbv)); Ga = array(0,c(nbh+1,nbv));
# vertical flux
# for (k in 1:nbv)
# { for (i in 2:nbh) {
#  Ft[k,i] = (u[k,i]*(Carr[k,i-1,NoS]-Carr[k,i,NoS]))
#  Fa[k,i] = (0.5*(sign(U[k,i])+1)*U[k,i]*Carr[k,i-1,NoS]+0.5*(1-sign(U[k,i]))*U[k,i]*Carr[k,i,NoS])
# }
# }
#
# horizontal flux
# for (k in 2:nbv)
# { for (i in 1:nbh) {
#   Gt[k,i] = (w[k,i]*(Carr[k-1,i,NoS]-Carr[k,i,NoS]))
#  Ga[k,i] = (0.5*(sign(W[k,i])+1)*W[k,i]*Carr[k-1,i,NoS]+0.5*(1-sign(W[k,i]))*W[k,i]*Carr[k,i,NoS])
# }
# }
#
# for interface flux
# for (i in 1:nbh) 
# {
# Gt[nbv+1,i] = (w[nbv+1,i]*(Carr[nbv,i,NoS]-cBgarr[NoS]))
# }

# browser()
# return(list(Carr, fw=f, fu=g, we=w, wa=W, ue=u, ua=U, Amatrix=a, d=d, residual=C))
## **************************  COMPLETE!  ***************************** ##
