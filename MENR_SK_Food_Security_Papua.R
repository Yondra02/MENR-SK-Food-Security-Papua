library(pracma) 
# Read the CSV file
data <- read.csv(file.choose(), header=TRUE, sep=";", dec=".")

# Assigning the variables
# Response variable
y <- data$Food.Security.Index
# Predictor variables
x1 <- data$Mean.Years.of.Schooling
x2 <- data$Percentage.of.Poor.Population
x3 <- data$Percentage.of.Households.with.Access.to.Proper.Sanitation
x4 <- data$Stunting.Prevalence

# Membuat data frame
df <- data.frame(y, x1, x2, x3, x4)

# Model regresi linier
model_lm <- lm(y ~ x1 + x2 + x3 + x4, data = df)

# Load package
library(car)

# Hitung VIF
vif(model_lm)


library(ggplot2)
# X1: Mean Years of Schooling vs Food Security Index
ggplot(df, aes(x = x1, y = y)) +
  geom_point() +
  ggtitle("Mean Years of Schooling vs Food Security Index") +
  xlab("Mean Years of Schooling") +
  ylab("Food Security Index")


# X2: Percentage of Poor Population vs Food Security Index
ggplot(df, aes(x = x2, y = y)) +
  geom_point() +
  ggtitle("Percentage of Poor Population vs Food Security Index") +
  xlab("Percentage of Poor Population") +
  ylab("Food Security Index")


# X3: Access to Proper Sanitation vs Food Security Index
ggplot(df, aes(x = x3, y = y)) +
  geom_point() +
  ggtitle("Access to Proper Sanitation vs Food Security Index") +
  xlab("Percentage of Households with Access to Proper Sanitation") +
  ylab("Food Security Index")


# X4: Stunting Prevalence vs Food Security Index
ggplot(df, aes(x = x4, y = y)) +
  geom_point() +
  ggtitle("Stunting Prevalence vs Food Security Index") +
  xlab("Stunting Prevalence") +
  ylab("Food Security Index")

#=========================================================================================
library(pracma) 

data
y=data[,3] #variabel y
x=as.matrix(data[c(4:7)]) 
xk=as.matrix(x[,c(1)]) #variabel kernel
xs=as.matrix(x[,c(2,3,4)]) #variabel spline 

kn=1 #Titik Knot Dicobakan 1 -3 

n=length(y) #jumlah pengamatan 
pk=ncol(xk) #jumlah variabel kernel 
ps=ncol(xs) #jumlah variabel spline

int.kr=30 #jumlah pembagi titik bandwidth yang diinginkan 
int.sp=30 #jumlah pembagi titik knot yang diinginkan 
alpha=0.05

#matrix
m1.nn=matrix(1, nrow=n, ncol=n)	#matriks 1 nxn 
m1.n1=matrix(1, nrow=n)		#matriks 1 nx1
mi.nn=diag(1,n,n)		#matriks identitas nxn

#penentuan titik knot 
knot=matrix(0,int.sp,ps) 
for (i in 1:ps)
{
  knot[,i]=seq(min(xs[,i]),max(xs[,i]),length.out=int.sp)
}
knot=as.matrix(knot[2:(int.sp-1),]) 
nknot=nrow(knot)

if (kn==1){ 
  knot=as.matrix(knot)
}else if (kn==2)
{
  #knot2 
  nkomb=(nknot*(nknot-1)/2)
  knot2=matrix(0,nkomb,kn*ps) 
  v=1
  for (i in 1:(nknot-1))
  {
    for (j in (i+1):nknot)
    {
      kk=0
      for (l in 1:ps)
      {
        a=cbind(knot[i,l],knot[j,l]) 
        kk=cbind(kk,a)
      }
      knot2[v,]=kk[1,2:ncol(kk)]
      v=v+1
    }
  }
  knot=as.matrix(knot2) 
  nknot=nrow(knot)
}else
{
  #knot3
  nkomb=(nknot*(nknot-1)*(nknot-2)/6) 
  knot3=matrix(0,nkomb,kn*ps)
  v=1
  for (i in 1:(nknot-2))
  {
    for (j in (i+1):(nknot-1))
    {
      for (k in (j+1):nknot)
      {
        kk=0
        for (l in 1:ps)
        {
          a=cbind(knot[i,l],knot[j,l],knot[k,l]) 
          kk=cbind(kk,a)
        }
        knot3[v,]=kk[1,2:ncol(kk)] 
        v=v+1
      }
    }
  }
  knot=as.matrix(knot3) 
  nknot=nrow(knot)
}

#penentuan bandwidth 
bw=matrix(0,int.kr,pk) 
for (i in 1:pk)
{
  bw[,i]=seq(0,(max(xk[,i])-min(xk[,i])),length.out=int.kr)
}
bw=as.matrix(bw[2:(int.kr-1),]) 
nband=nrow(bw)

#desain matriks X(k) pada spline 
MSE=matrix(0,nband*nknot) 
GCV=matrix(0,nband*nknot) 
code=matrix(0,nband*nknot,kn*ps+pk) 
o=1
for (i in 1:nknot)
{
  for (j in 1:nband)
  {
    #matrik spline 
    Z=cbind(1,xs) 
    a=1
    for (k in 1:ps)
    {
      for (l in 1:kn)
      {
        Z=cbind(Z,(pmax(0,xs[,k]-knot[i,a]))) 
        a=a+1
      }
    }
    sum.v.phi=0  
    for (k in 1:pk)
    {
      v.diag=diag(xk[,k]) 
      V=m1.nn %*%v.diag 
      z=(t(V)-V)/bw[j,k]
      K=1/sqrt(2*pi)*exp(-1/2*z^2) #fungsi kernel gaussian
      K.Z=(1/bw[j,k])*K
      W.penyebut=diag(c(1/n*K.Z%*%m1.n1))%*%m1.nn 
      V.phi=1/n*K.Z/W.penyebut
      
      #penimbang V(phi).1
      sum.v.phi=sum.v.phi+V.phi	#nilai kernel untuk setiap variabel
    }
    # penimbang kernel gabungan 
    V.phi=sum.v.phi/pk 		#nilai kernel rata-rata
    
    #estimasi parameter 
    beta=0
    C=pinv(t(Z)%*%Z)%*%t(Z)%*%(mi.nn-V.phi) 
    beta=C%*%y
    A=Z%*%C
    B=A+V.phi 
    yhat=B%*%y 
    error=y-yhat
    MSE[o]=n^-1*(t(error)%*%error)
    db=(sum(diag(mi.nn-A-V.phi))/(n))^2 
    GCV[o]=MSE[i]/(1-db)
    code[o,]=c(knot[i,],bw[j,]) 
    o=o+1
  }
}
optimum=cbind(code,MSE,GCV)
optimum
GCVmin=optimum[order(optimum[,(kn*ps+pk+2)]),] 
GCVmin

#mengurutkan nilai GCV minimum 
knot.opt=GCVmin[1,1:(kn*ps)] 
band.opt=GCVmin[1,(kn*ps)+1:ncol(xk)] 
gcv.opt=GCVmin[1,ncol(GCVmin)]
#validasi nilai GCV terkecil 
#matrik  spline 
Z=cbind(1,xs)
a=1
for (k in 1:ps)
{
  for (l in 1:kn)
  {
    Z=cbind(Z,(pmax(0,xs[,k]-knot.opt[a]))) 
    a=a+1
  }
}
sum.v.phi=0  
for (k in 1:pk)
{
  v.diag=diag(xk[,k]) 
  V=m1.nn %*%v.diag 
  z=(t(V)-V)/band.opt[k]
  K=1/sqrt(2*pi)*exp(-1/2*z^2) #fungsi kernel gaussian
  K.Z=(1/band.opt[k])*K 
  W.penyebut=diag(c(1/n*K.Z%*%m1.n1))%*%m1.nn 
  V.phi=1/n*K.Z/W.penyebut
  #penimbang V(phi).1
  sum.v.phi=sum.v.phi+V.phi #nilai kernel rata-rata setiap variabel
}

# penimbang kernel gabungan
V.phi=sum.v.phi/pk	#nilai kernel rata-rata


#estimasi parameter 
beta=0
C=pinv(t(Z)%*%Z)%*%t(Z)%*%(mi.nn-V.phi) 
beta=C%*%y
A=Z%*%C
B=A+V.phi 
yhat=B%*%y 
error=y-yhat
db=matrix(NA,nrow=3) 
SS=matrix(NA,nrow=3) 
MS=matrix(NA,nrow=3)
PValue=matrix(NA,nrow=3)
deci=matrix(0,nrow=(ps*(kn+1)+1)) 
Fhitung=matrix(NA,nrow=3) 
db[1]=ncol(Z)
db[2]=n-db[1]-1 
db[3]=n-1
SS[1]=sum((yhat-mean(y))^2)
SS[2]=sum((y-yhat)^2)
SS[3]=sum((y-mean(y))^2) 
MS[1]=SS[1]/db[1]
MS[2]=SS[2]/db[2] 
R2=(SS[1]/(SS[1]+SS[2]))*100
MAPE=((sum((abs(y-yhat))/y))/n)*100

#Uji F (Uji Serentak) 
Fhitung[1]=MS[1]/MS[2]
PValue[1]=pf(Fhitung[1],db[1],db[2],lower.tail=FALSE)
ANOVA=cbind(db,SS,MS,Fhitung,PValue) 
colnames(ANOVA)=c("db","SS","MS","Fhitung","P-Value")
rownames(ANOVA)=c("Regresi","Error","Total") 
Ftabel=qf(0.95,db[1],db[2])
if (Fhitung[1]>Ftabel)
{dec='H0 ditolak'
}else
  dec='H0 gagal ditolak'
ANOVA

library(MASS)

#Uji Parsial
n1=nrow(beta)
thit=rep(NA,n1)
pvaluee=rep(NA,n1)
bag1=(mi.nn-t(mi.nn-V.phi)%*%Z%*%pinv(t(Z)%*%Z)%*%t(Z))
bag2=(mi.nn-(Z%*%pinv(t(Z)%*%Z)%*%t(Z)%*%(mi.nn-V.phi)))
MSE=(t(y)%*%bag1%*%bag2%*%y)/n
A=pinv(t(Z)%*%Z)%*%t(Z)%*%(mi.nn-V.phi)
baru = (A%*%t(A)) / MSE[1,1] #(MSE)
SE=sqrt(diag(baru))
for (i in 1:n1)
{
  thit[i]=beta[i,1]/SE[i]
  pvaluee[i]=round(2*(pt(abs(thit[i]),(n-n1),lower.tail=FALSE)),5)
}
thit=as.matrix(thit)
colnames(thit)<-"thitung"
colnames(beta)<-"parameter untuk Beta"

tg1=cbind(beta,thit,pvaluee)
tg1

knot.opt
band.opt
gcv.opt
R2
MAPE
MSE = mean(SS[2])
beta

# ============================================================
# ACTUAL VS PREDICTED FOOD SECURITY INDEX
# ============================================================

library(ggplot2)

# Create data frame for plotting
data_plot <- data.frame(
  Observation = 1:length(y),
  Actual = y,
  Predicted = yhat
)

# Plot Actual vs Predicted
ggplot(data_plot, aes(x = Observation)) +
  geom_line(
    aes(y = Actual, color = "Actual"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = Predicted, color = "Predicted"),
    linewidth = 1
  ) +
  labs(
    title = "Actual vs Predicted Food Security Index",
    x = "Observation",
    y = "Food Security Index",
    color = NULL
  ) +
  scale_color_manual(
    values = c(
      "Actual" = "black",
      "Predicted" = "red"
    )
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank()
  )