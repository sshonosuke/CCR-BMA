source("CCR-function-new.R")
library(mada)

# Data
data(AuditC)
dd=madad(AuditC)

logit=function(x){ log(x/(1-x)) }
invlogit=function(x){ exp(x)/(1+exp(x)) }

pp=cbind(dd$sens$sens,1-dd$spec$spec)
y=logit(pp)
ni=cbind(AuditC$TP+AuditC$FN,AuditC$FP+AuditC$TN)
S=1/(pp*(1-pp)*ni)


# Proposed confidence region
fit=CCR(y,S,alpha=0.05)

h.ast=fit[[2]]
hbeta=fit[[3]]
IH=fit[[4]]
xx=fit[[5]]

Theta=seq(0,2*pi,length=200)
Lam=eigen(IH)$values
HH=t(eigen(IH)$vectors)
t(HH)%*%diag(Lam)%*%HH

c=xx*(1+h.ast)
tx=sqrt(c/Lam)*rbind(cos(Theta),sin(Theta))
CCR=t(hbeta+t(HH)%*%tx)



# Reitsuma's method
fits=reitsma(AuditC,method="reml")

QQ <- apply(diag(fits$Psi)*t(1/S), 1, sum)
m <- dim(ni)[1]
(QQ-(m-1))/QQ*100

# load the result of Sugasawa and Noma (2019)
load("AuditC-result.RData")

# Plot
postscript("AuditC.eps",height=8,width=8,pointsize=16,horizontal=F)
plot(fits,sroclwd=2,main ="",
     ylim=c(0.5,1),xlim=c(0,0.52))
points(invlogit(CCR)[,c(2,1)],type="l",lty=2)
points(invlogit(sECR)[,c(2,1)],type="l",lty=3)
points(invlogit(y)[,c(2,1)],pch=2)
legend("bottomright",c("Data point","Summary estimate","SROC","Approximate CR","Corrected CR","Exact CR"),
       pch=c(2,1,NA,NA,NA,NA),lwd=c(NA,NA,2,1,1,1),lty=c(NA,NA,1,1,2,3))
dev.off()


pdf("AuditC.pdf",height=8,width=8,pointsize=16)
plot(fits,sroclwd=2,main ="",
     ylim=c(0.5,1),xlim=c(0,0.52))
points(invlogit(CCR)[,c(2,1)],type="l",lty=2)
points(invlogit(sECR)[,c(2,1)],type="l",lty=3)
points(invlogit(y)[,c(2,1)],pch=2)
legend("bottomright",c("Data point","Summary estimate","SROC","Approximate CR","Corrected CR","Exact CR"),
       pch=c(2,1,NA,NA,NA,NA),lwd=c(NA,NA,2,1,1,1),lty=c(NA,NA,1,1,2,3))
dev.off()
