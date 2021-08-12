# Hybrid_RadialCheb_lsoda.R November 3 2006
# RADIAL FULL model on [0,1]

require(odesolve);

########################################################### 
#  compute D = differentiation matrix, x = Chebyshev grid   
###########################################################  
cheb.inc=function(N) {
  x = matrix(cos(pi*(0:N)/N),N+1,1); 
  c = c(2,rep(1,N-1),2)*((-1)^(0:N)); c=matrix(c,N+1,1); 
  X = matrix(x,N+1,N+1);
  dX = X-t(X);                  
  D  = (c%*%t(1/c))/(dX+diag(rep(1,N+1))) # off-diagonal entries
  D  = D - diag(apply(t(D),2,sum));       # diagonal entries
  return(list(Dmatrix=D[(N+1):1,(N+1):1],grid=rev(x)))
}

# Transform grid and matrix onto interval [lower, upper]
cheb.moved=function(N,lower=0,upper=1) {
	g=cheb.inc(N); x=g$grid; D=g$Dmatrix;
	x.moved=lower+(upper-lower)*(x+1)/2; 
	D.moved=(2/(upper-lower))*D;
	return(list(Dmatrix=D.moved,grid=x.moved))
}

# Compute advective transport term for flux f(r)  
rflux=function(f) {f/rvals + D%*%f}

###########################################################
#  Radial fungal model on r=[0,1] with QSS solution for signal 
###########################################################
RFUNGAL=function(t,y,parms){
	 ymat=matrix(y,nx-1,Ncritters); 
	 # impose zero-derivative BC at r=0 
	 ytop=ymat[1,]; 
	 yu=rbind(matrix(ytop,1,Ncritters),ymat); 
	 
# State Variables
         T=yu[,1]; P=yu[,2]; U=yu[,4]; S=yu[,3]; 
	 G=yu[,5]; W=yu[,6];
  	 feed = pmax(0,P*T);   

         # QSS signal, after transient period 
	 if(Sswitch>0){S=solve(Smatrix,feed*endzero); S[S<0]=0}  

# First derivatives
	 DT=D%*%yu[,1]; DS=D%*%S; Ucells=S*U; 

      	 dy[,1]= -q*feed;

         dy[,2]= q*feed - (mG*G + 0.5*mG*W + muP)*P - etaP*rflux(dvec*P*DT) + (etaP*parms[1])*(DXD%*%P)  

	 dy[,3]= P*T - muS*S + Ds*rflux(DS); if(Sswitch>0) dy[,3]=0;    
	
	 dy[,4]= muU*(1-U) - (deltaG+deltaW)*Ucells - etaI*rflux(U*DS)+ (etaI*parms[2])*(DXD%*%U)  

         dy[,5]= deltaG*Ucells - muI*G - etaI*rflux(G*DS) + (etaI*parms[2])*(DXD%*%G)   

         dy[,6]= deltaW*Ucells - muI*W - etaI*rflux(W*DS) + (etaI*parms[2])*(DXD%*%W)   

         dy[,7]= W;

	 # impose zero-flux BC at r=1 by finite difference 
	 dy[nx,2:6]=dy[nx-1,2:6]; 
	 
# send to the integrator derivatives at nonzero grid points 
         return(list(dy=as.vector(dy[-1,])));
}

#######	 Set up grid and derivative matrices
#######  Grid is chosen so two leftmost pts are at -eps, eps 
N=50; #number of grid points  

# function that is zero when two leftmost pts are at -eps, eps
endfunc=function(r0) {
	cn=cheb.moved(N,lower=r0,upper=1)$grid;
	return(cn[1]+cn[2])
}
# find the right grid 
cn=cheb.moved(N,lower=0,upper=1); eps=cn$grid[2];  
r0=uniroot(endfunc,lower=-eps,upper=0)$root;  
cn=cheb.moved(N,lower=r0,upper=1); 
rvals=cn$grid; nx=length(rvals); r.inv=1/rvals; 

###  Vectors to zero out endpoints of a state vector, 
###  and Spatially Varying Pathogen Advection 
endzero=c(0,rep(1,nx-2),0); dvec=(1-rvals^20)^12;

# First and second derivative matrices 
D=cn$Dmatrix; D2=D%*%D; 

###### Spatially varying "artificial diffusion" matrix
dcoef=pmax(0,1-3.8*(rvals-0.5)^2)
### 1-dimensional model: 
# DXD=D%*%diag(dcoef)%*%D; 
### In polar coordinates: 
DXD=diag(r.inv)%*%D%*%diag(rvals*dcoef)%*%D; 

## Ncritters=Number of State Variables
Ncritters=7; dy=matrix(0,nx,Ncritters);   

# Initial conditions 
Po=0.5*pmax(0,(1-(12*rvals)^6)^3); P0=Po[2:nx];
T0=1-P0; S0=rep(0,nx-1); 
U0=rep(1,nx-1); G0=rep(0,nx-1); W0=rep(0,nx-1); B0=rep(0,nx-1);
y0=c(T0,P0,S0,U0,G0,W0,B0); 


############### END SETUP 

############## Run the model 
### Default Parameter values 
q=.05; muP = 0.1*q; etaP = 0.0005;  		# Pathogen 
Ds=0.1;   muS=.1;  				# Signal 
etaI=1;   muU=0.1; 				# Undifferentiated 
deltaG=1; deltaW=1; muI=0.1; mG=0.02;		# Gobblers and Wallers

### second set of default params
# q=0.02; muP = 0.1*q; mG=0.002; 
	 		
# matrix for calculating QSS signal distribution 
Smatrix= - Ds*diag(r.inv)%*%D%*%diag(rvals)%*%D; 
diag(Smatrix)=diag(Smatrix)+muS; 
Smatrix[nx,]=D[nx,]; 
Smatrix[1,]=c(1,-1,rep(0,nx-2)); 

### Integrate 
dt=1; tvals=seq(from=0,to=180,by=dt); nt=length(tvals); tS=500; 
yt=matrix(0,nt,length(y0)); yt[1,]=y0; 
for(j in 2:nt) {
	Sswitch=ifelse((j*dt)>tS,1,0)  
	yt.n1=lsoda(yt[j-1,],times=tvals[(j-1):j],func=RFUNGAL,parms=c(0.1,0.05),hmax=.1);
	yt.n1=yt.n1[2,-1]; 
	if((dt*(j-1))%%2==1) cat((j-1)*dt,sum(is.na(yt.n1)),"\n"); 
	yt[j,]=yt.n1; 
}

### Compute metrics
jrow=length(tvals); Y=matrix(yt[jrow,],nx-1,Ncritters);
if(jrow>tS) {
  feed=Y[,1]*Y[,2]; 
  S=solve(Smatrix,c(0,feed)*endzero); S[S<0]=0
  Y[,3]=S[-1]; 
} 
xvals=c(-rvals[nx:2],rvals[2:nx]);   
Y=rbind(Y[(nx-1):1,],Y[1:(nx-1),]); 

# interpolate and integrate to get total P and T 	
Pfun=approxfun(xvals,abs(xvals)*Y[,2]); Pfinal=2*integrate(Pfun,0,1,subdivisions=500)$value	
Tfun=approxfun(xvals,abs(xvals)*Y[,1]); Tfinal=2*integrate(Tfun,0,1,subdivisions=500)$value
Bfun=approxfun(xvals,abs(xvals)*Y[,7]); Bfinal=2*integrate(Bfun,0,1,subdivisions=500)$value

############ END RUN 

### Plotting 
win.graph(); par(mfrow=c(3,2));
plotcolors=c("green","red","cyan","black","darkgreen")
plotcols=c(1,2,3,4,5); plottimes=c(0,10,30,60,120,180); 
for(j in 1:6){
   jrow=which(tvals==plottimes[j]); 
   Y=matrix(yt[jrow,],nx-1,Ncritters);
   if(plottimes[j]>tS) {
	feed=Y[,1]*Y[,2]; 
	S=solve(Smatrix,c(0,feed)*endzero); S[S<0]=0
	Y[,3]=S[-1]; 
   } 
   xvals=c(-rvals[nx:2],rvals[2:nx]);   
   Y=rbind(Y[(nx-1):1,],Y[1:(nx-1),]); 

   #### interpolate and integrate to get total P and T 	
   Pfun=approxfun(xvals,abs(xvals)*Y[,2]); totalP=2*integrate(Pfun,0,1,subdivisions=500)$value	
   Tfun=approxfun(xvals,abs(xvals)*Y[,1]); totalT=2*integrate(Tfun,0,1,subdivisions=500)$value

   matplot(xvals,Y[,plotcols],type="l",ylim=range(c(0,1,as.vector(Y[,plotcols]))), 
   lwd=2, xlab="Radial distance r", ylab=" ",col=plotcolors[plotcols], lty=1)
   title(paste("t=",tvals[jrow]," Total P=",round(totalP,digits=3)," Total T=",round(totalT,digits=3)  )); 
}

