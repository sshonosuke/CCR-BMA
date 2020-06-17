source("CCR-BMA-function.R")
library(mada)

# Data
data(AuditC)
dd <- madad(AuditC)

logit <- function(x){ log(x/(1-x)) }
invlogit <- function(x){ exp(x)/(1+exp(x)) }

pp <- cbind(dd$sens$sens, 1-dd$spec$spec)
Y <- logit(pp)
ni <- cbind(AuditC$TP+AuditC$FN, AuditC$FP+AuditC$TN)
S <- 1/(pp*(1-pp)*ni)
n <- dim(Y)[1]


# Confidence region
CR <- CCR.BMA(Y, S, alpha=0.05)

# Reitsuma's method
fits <- reitsma(AuditC, method="reml")

# heterogeneity measure
ww <- 1/S
QQ <- (n-1)*apply(ww, 2, sum) / (apply(ww, 2, sum)^2 - apply(ww^2, 2, sum))
hPsi <- diag(fits$Psi)
hPsi/(QQ+hPsi)*100


# Plot
plot(fits, sroclwd=2, ylim=c(0.5,1), xlim=c(0,0.52))
points(invlogit(CR)[,c(2,1)], type="l", lty=2)
points(invlogit(Y)[,c(2,1)], pch=2)
meth <- c("Data point","Summary estimate","SROC","Approximate CR","Corrected CR")
legend("bottomright", legend=meth, pch=c(2,1,NA,NA,NA), lwd=c(NA,NA,2,1,1), lty=c(NA,NA,1,1,2))
