## program for analysis of banteng population projection

## ----------------------------------------
## base model - rates from Choquenot (1993)
## giving stability
## ----------------------------------------

## Removes everything in memory
rm(list = ls())

## Add packages 'base'

library(base)

## define max age
age.max<-17

## define survival - fitted function (Choquenot 1993)
## s0<-0.7402 ## fitted from Choquenot (1993)
## s.vec <- c(0.55,0.9278,0.9346,0.9370,0.9382,0.9389,0.9394,0.9397,0.9400,0.9402,0.9403,0.9405,0.9406,0.9407,0.9408,0.9408,0.9409,0.9409)

## Mortality vector from Choquenot (1993)
q.vec <- c(0.26,0.08,0.075,0.07,0.075,0.06,0.06,0.06,0.06,0.05,0.05,0.05,0.07)

## Fit survival values
library(stats)
library(base)

## Create data frame
age.vec.q <- seq(-1,length(q.vec)-2)
age.vec.q[1] <- 0
age.vec.q[2:(length(age.vec.q))] <- age.vec.q[2:(length(age.vec.q))] + 0.5

q.fit.data <- data.frame(age.vec.q,q.vec)

## model formula (exponential rise to maximum): yhat = a*exp(b/(x+c))

fit.q.vec <- nls(q.fit.data$q.vec ~ a.coeff * exp(b.coeff/(q.fit.data$age.vec.q + c.coeff)),
		data = q.fit.data,
		start = list(a.coeff = 0.1, b.coeff = 1, c.coeff = 1),
		trace = TRUE)
sum.fit.q.vec <- summary(fit.q.vec)

## Coefficients from fit
coeff.fit.q.vec <- as.numeric(sum.fit.q.vec$parameters)

## age vector
age.vec<-rep(0,age.max+1)
	for (j in 1:(age.max)) {
	    age.vec[j+1] <- j
	}

## Predict new q vector with integer age inputs
qf.new.vec <- coeff.fit.q.vec[1] * exp(coeff.fit.q.vec[2]/(coeff.fit.q.vec[3]+age.vec))

row <- 1
col <- 1
par(mfrow=c(row,col))
plot(age.vec.q,q.vec,xlab="age (years)",ylab="qx")
title(main="Fitted Survival")
lines(age.vec,qf.new.vec)
par(mfrow=c(1,2))

s.vec <- 1 - qf.new.vec

## change s0 in s.vec to give population stability
##s.vec[1] <- 0.62

## fitted fecundity parameters
## m.vec <- c(0.0000,0.0001,0.0133,0.0947,0.2303,0.3440,0.3964,0.3921,0.3525,0.2979,0.2414,0.1903,0.1471,0.1122,0.0849,0.0638,0.0478,0.0358)
m.vec.obs <- c(0,0,0,0,0.22,0.351,0.5,0.45,0.2,0.325,0.125,0.225) ## Choquenot (1993)

## Create data frame
age.vec.m <- seq(0,11)

m.fit.data <- data.frame(age.vec.m,m.vec.obs)

## model formula (log-normal): yhat=a*exp(-0.5*(ln(x/x0)/b)^2) 

fit.m.vec <- nls(m.fit.data$m.vec.obs ~ a.coeff.m * exp(-0.5*(log(m.fit.data$age.vec.m/x0.coeff.m)/b.coeff.m)^2),
		data = q.fit.data,
		start = list(a.coeff.m = 0.4, x0.coeff.m = 6, b.coeff.m = 0.4),
		trace = TRUE)
sum.fit.m.vec <- summary(fit.m.vec)

## Coefficients from fit
coeff.fit.m.vec <- as.numeric(sum.fit.m.vec$parameters)

## Predict new q vector with integer age inputs
mf.new.vec <- coeff.fit.m.vec[1] * exp(-0.5*(log(age.vec/coeff.fit.m.vec[2])/coeff.fit.m.vec[3])^2)

## re-scale m.vec.obs for plotting
m.vec.obs.rs <- rep(0,age.max+1)
m.vec.obs.rs[1:12] <- m.vec.obs

row <- 1
col <- 1
par(mfrow=c(row,col))
plot(age.vec,m.vec.obs.rs,xlab="age (years)",ylab="mx",xlim=range(c(0,age.max)),ylim=range(c(0,1)))
title(main="Fitted Fertility")

lines(age.vec,mf.new.vec)
par(mfrow=c(1,2))

m.vec <- mf.new.vec

primi<-4 ## age at first reproduction (females)

## initial population sizes
## Nmin<-4000 ## J. Christopherson (2004)
## Nmax<-5000 ## J. Christopherson (2004)
Nmin<-7000 ## K. Saalfeld (2004)
Nmax<-8000 ## K. Saalfeld (2004)
Navg <- mean(c(Nmin, Nmax))

## Use 2 scenarios - pop.init at 4000 and 8000
N1 <- 5000
N2 <- 9000

## Density estimates
Area.Cob <- 220000 ## Area of Cobourg Peninsula in hectares
D.avg1 <- N1/(Area.Cob*0.01)
D.avg2 <- N2/(Area.Cob*0.01)

## sex ratios
sr <- 0.5 ## overall sex ratio
x <- 0.5 ## foetal sex ratio (Choquenot 1993)

##total females in population
f1<-N1*sr
f2<-N2*sr

##total males in population
mal1<-N1*(1-sr)
mal2<-N2*(1-sr)

##Initial population size vector
N<-N1

## Matrix size
k <- 34

## Create matrix shell
a <- matrix(data<-0,nrow<-k,ncol<-k)

## Add survival vectors to matrix
diag(a[2:age.max,1:(age.max-1)]) <- s.vec[2:age.max] ## adds s.vec to female quadrant
diag(a[(age.max+2):k,(age.max+1):(k-1)]) <- s.vec[2:age.max] ## adds s.vec to male quadrant
a[1,1:(age.max)] <- s.vec[1] * m.vec[1:(age.max)] ## adds female fecundity
a[(age.max+1),1:(age.max)] <- s.vec[1] * m.vec[1:(age.max)] ## adds male fecundity

## Maximum lambda function
max.lambda <- function(x) Re((eigen(x)$values)[1]) ## where 'x' is a Leslie matrix

## Maximum r function
max.r <- function(x) log(Re((eigen(x)$values)[1])) ## where 'x' is a Leslie matrix

## Stable stage distribution
stable.stage.dist <- function(x) ((x %*% (Re((eigen(x)$vectors)[,1])))/(sum((x %*% (Re((eigen(x)$vectors)[,1]))))))[,1]

## Generation length function

R.val <- function(X,age.max) ## reproductive value (R0) where X = Leslie matrix; age.max = maximum age of females
{		
		## define the transition matrix
		T <- X[1:age.max,1:age.max]
		T[1,1:(age.max)] <- 0

		## define the fertility matrix
		F <- X[1:age.max,1:age.max]
		diag(F[2:age.max,1:(age.max-1)]) <- 0

		## define the identity matrix
		I <- matrix(data<-0,nrow<-age.max,ncol<-age.max)
		diag(I) <- 1

		## define the fundamental matrix
		library(MASS)
		N.fund <- ginv(I-T)

		## define the reproductive matrix
		R <- F %*% N.fund

		## define R0 (number of female offspring produced per female during lifetime)
		R0 <- Re((eigen(R)$values)[1])
		
		## output
		print("number of female offspring produced per female during its lifetime")
		print("_________________________________________________________________")
		print(R0)

}

## Mean generation time function
G.val <- function (X,age.max) ## where X is a Leslie Matrix
{	
		G <- (log(R.val(X,age.max)))/(log(Re((eigen(X)$values)[1])))
		print("mean generation time")
		print("____________________")
		print(G)
}

## max lambda
maxlambda<-max.lambda(a)

## r
r <- max.r(a)

## Reproductive value
R0 <- R.val(a,age.max)

## Generation time
G <- G.val(a,age.max)

## Stable stage distribution
ssd <- stable.stage.dist(a)

## ssd classes
    ## female
    ssd.juvf <- sum(ssd[1:primi-1])
    ssd.adf <- sum(ssd[primi:age.max])

    ## male
    ssd.juvm <- sum(ssd[(age.max+1):(age.max+primi-1)])
    ssd.adm <- sum(ssd[(age.max+primi):k])

## Average age of mothers at stability
G.age.mean <- weighted.mean(age.vec[2:(age.max+1)],ssd[1:age.max])

## pick initial vectors
n<-ssd*N

##Calculate Quasi-extinction times
thresh <- 50 ## Define quasi-exinction threshold

Q<-(log(thresh/sum(n)))/log(maxlambda)
	if (Q < 0) {
	    Q <- "infinity"
	}
	Q <- Q

## do a simulation
## first specify the initial condition and length of simulation
tlimit <- 30

##set population size year step vector
pop.vec <- rep(1,tlimit+1)
pop.vec[1] <- sum(n)

##set year step vector
yr.vec <- rep(1,tlimit+1)
yr.vec[1] <- 0

	for (j in 1:tlimit) {
    		yr.vec[j+1] <- j
	}

##then iterate

	for (ii in 1:tlimit) { 
	   n <- a %*% n 
	   pop.vec[ii+1] <- sum(n)
	} 

log.pop.vec <- log10(pop.vec)

##total population size after 'tlimit' years
pop.st <- N
pop.end <- sum(n)
tlimit
maxlambda
r
Q
R0
G

## continue displays
ssd.juv <- ssd.juvf + ssd.juvm
ssd.ad <- ssd.adf + ssd.adm

##Make density independent plots
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yr.vec,pop.vec,xlab="year",ylab="N",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(pop.vec)+(0.25*(max(pop.vec))))),type='l')
plot(yr.vec,log.pop.vec,xlab="year",ylab="log N",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(log.pop.vec)+(0.25*(max(log.pop.vec))))),type='l')

##Make survival and fecundity plots
plot(age.vec,s.vec,xlab="age (years)",ylab="P(survival)",xlim=range(0.5:(age.max+0.5)),ylim=range(0:1),type='l')
plot(age.vec,m.vec,xlab="age (years)",ylab="m (fertility)",xlim=range(0.5:(age.max+0.5)),ylim=range(0:1),type='l')
par(mfrow=c(1,2))




## *******************************************
## Density-dependence in survival & fertility
## Set kill rates
## *******************************************

## ***********************************************************************************************
##K population size vector
Nd <- 9000

## Set time limit for projection
tlimit <- 30

## Set kill type (2 = adult females (>2 years) only; 1 = adult males (>5 years) only; 0 = adults (fem >2, mal >3) only)
kill.type <- 0

if (kill.type == 0) scenario <- "adult kill"
if (kill.type == 1) scenario <- "adult male (>5 years) kill"
if (kill.type == 2) scenario <- "adult female (>2 years) kill"

## Set kill range
min.kill <- 2 ## Set minimum kill rate
max.kill <- 1000 ## Set maximum kill rate

## Set kill interval
kill.interval <- 50 ## (calculate kill in increments of x)

## Set target per cent reduction of original population size
thresh <- 0

## Stochasticity in rainfall
## 0 = deterministic; 1 = stochastic
stoch.select <- 1
rain.cv = 24 ## coefficient of variation in rainfall (Choquenot 1993)

## Set number of iterations to estimate mean & confidence interval for projections based on flooding stochasticity
perm <- 1000 ## number of permutations
if (stoch.select == 0) perm <- 1 else perm <- perm

## Tertiary sex ratio (0.83 = balanced in polygynous ungulates)
tsr.max <- 0.95 ## Tertiary sex ratio (breeding females:breeding adults)

######################################
## Economic analysis
######################################

## Exchange rate (AU$ to Euros); updated live from server
##exch.aud.euro <- 0.587673

exch.rates <- download.file("http://www.securetrading.net/download/currencyrates/curr-gbp-ns-today.txt","C:\\Documents and Settings\\c_bradshaw\\My Documents\\Biology\\Banteng\\R\\exch\\exch.txt", "internal", quiet = FALSE, mode="w")
exch.rat.dat <- scan(file="C:\\Documents and Settings\\c_bradshaw\\My Documents\\Biology\\Banteng\\R\\exch\\exch.txt",what="character",skip=3,nlines=1)
exch.scan.aud <- scan(file="C:\\Documents and Settings\\c_bradshaw\\My Documents\\Biology\\Banteng\\R\\exch\\exch.txt",what="character",skip=7,nlines=1)
##exch.rates <- download.file("http://www.securetrading.net/download/currencyrates/curr-gbp-ns-today.txt","F:\\Banteng\\R\\exch\\exch.txt", "internal", quiet = FALSE, mode="w")
##exch.rat.dat <- scan(file="F:\\Banteng\\R\\exch\\exch.txt",what="character",skip=3,nlines=1)
##exch.scan.aud <- scan(file="F:\\Banteng\\R\\exch\\exch.txt",what="character",skip=7,nlines=1)

aud.gbp <- 1/((as.numeric(as.character((strsplit(exch.scan.aud," ")[2]))))/10000000)

exch.scan.eur <- scan(file="C:\\Documents and Settings\\c_bradshaw\\My Documents\\Biology\\Banteng\\R\\exch\\exch.txt",what="character",skip=11,nlines=1)
##exch.scan.eur <- scan(file="F:\\Banteng\\R\\exch\\exch.txt",what="character",skip=11,nlines=1)

eur.gbp <- 1/((as.numeric(as.character((strsplit(exch.scan.eur," ")[2]))))/10000000)
print(aud.gbp); print(eur.gbp)
exch.aud.euro <- aud.gbp/eur.gbp
print(exch.aud.euro)
print(exch.rat.dat)

## Eastern Young Cattle Indicator (EYCI) market data (AU$/kilo meat)
eyci <- download.file("http://www.mla.com.au/content.cfm?sid=1173","C:\\Documents and Settings\\c_bradshaw\\My Documents\\Biology\\Banteng\\R\\eyci\\eyci.txt", "internal", quiet = FALSE, mode="w")
eyci.lat <- scan(file="C:\\Documents and Settings\\c_bradshaw\\My Documents\\Biology\\Banteng\\R\\eyci\\eyci.txt",what="character",skip=641,nlines=1)
##eyci <- download.file("http://www.mla.com.au/content.cfm?sid=1173","F:\\Banteng\\R\\eyci\\eyci.txt", "internal", quiet = FALSE, mode="w")
##eyci.lat <- scan(file="F:\\Banteng\\R\\eyci\\eyci.txt",what="character",skip=643,nlines=1)

eyci.now <- (as.numeric(substr((as.character(eyci.lat[2])),23,28)))/100
print(eyci.now)

## Read in full range of EYCI since 1996
eyci.dat <- read.table("C:\\Documents and Settings\\c_bradshaw\\My Documents\\Biology\\Banteng\\R\\eyci\\eycidat.csv", sep=",", header=TRUE)
##hist(eyci.dat$adjeyci)
##plot(eyci.dat$year,eyci.dat$adjeyci,type="l")
yrf.rnd <- factor(floor(eyci.dat$year))
yrn.rnd <- (floor(eyci.dat$year))
eycidat <- data.frame(eyci.dat$adjeyci,yrf.rnd,yrn.rnd)

## Summary table
attach(eycidat)
eyci.tab <- table(yrf.rnd)
eyci.xtab <- xtabs(eyci.dat$adjeyci ~ yrf.rnd)
detach(eycidat)
eyci.yr.mean <- eyci.xtab/eyci.tab/100

##plot(eyci.yr.mean,type="l")

## Min and max eyci (since 1996, set in 2004 AU$)
eyci.min <- min(eyci.dat$adjeyci)/100 ## 1.8596
eyci.max <- max(eyci.dat$adjeyci)/100 ## 4.3598

## Temporal autocorrelation
lag1 <- as.numeric(eyci.yr.mean[2:(nrow(eyci.yr.mean))])
base <- as.numeric(eyci.yr.mean[1:(nrow(eyci.yr.mean)-1)])

eyci.fit <- lm(lag1 ~ base)

##plot(base,lag1)
##abline(eyci.fit,col="red")

## Mean & SD
eyci.mean <- as.numeric(mean(eyci.yr.mean))
eyci.var <- as.numeric(var(eyci.yr.mean))

## User-defined values
ppkg.meat <- 10 ## (price in AU$ to obtain 1 kg banteng-equivalent meat); T. Griffiths, pers. comm.

## Meat eaten per day per person
adult.meat.day <- 67 ## g/day
adult.meat.yr <- (adult.meat.day*365)/1000 ## kg/year/adult
minjil.pop <- 254 ## 2001 census data
minjil.juv <- 39 + 53
minjil.ad <- 35 + 69 + 42 + 16
minjil.meat <- (adult.meat.yr*minjil.ad) + ((adult.meat.yr*minjil.juv)/2)

## Proportion of full animal producing marketable meat (dressing proportion)
meat.prop <- 0.526 ## for banteng (Kirby 1979)

## Price per mature male harvested passed to Traditional owners
safari.to <- 2500 ## AU$

## Costing 
muster <- 10 ## per head
transport <- 10 ## per head (http://www.agric.nsw.gov.au/reader/an-transport/6980)

##cost.tot <- muster + transport ## set individually above OR,
cost.tot <- 150 ## set total cost in AU$ per head (Alister Trier, Proj. Mgr - NT Indigenous Pastoral Project

## Cost function
## Increasing costs as a function of a reduction in population size
cost.max <- 500 ## max cost (AU$) at 0.1*K
cost.min <- 150 ## min cost (AU$) at K

## Profit (http://www.abs.gov.au/websitedbs/c311215.NSF/0/641c9a7c2c02c6d9ca256d480012c608?OpenDocument)
profit <- 29 ## average % profit per $ invested (Australia average)

## Demand function
## Set the the probability (p.dem) of attracting max n safari hunters (n.saf) to kill trophy males
n.saf <- 40
p.dem <- 0.5


## *********************************************************************************
## ***************
## End user input
## ***************
## *********************************************************************************

## Calculate demand function (3-parameter Hill logistic decay)
## yhat = ((a * x^b)/(c^b + x^b))

## Create data.frame
p.dem.vec <- c(1,0.99,0.98,0.95,0.90,0.80,0.64,p.dem)
n.saf.vec <- c(1,n.saf/4,n.saf/2,n.saf*0.75,n.saf*0.85,n.saf*0.92,n.saf*0.95,n.saf)
dem.data <- data.frame(n.saf.vec,p.dem.vec)

fit.dem <- nls(dem.data$p.dem.vec ~ ((a.dem * ((dem.data$n.saf.vec)^b.dem))/((c.dem^b.dem)+(dem.data$n.saf.vec^b.dem))),
		data = dem.data,
		start = list(a.dem = 1, b.dem = -(n.saf/2), c.dem = n.saf),
		trace = TRUE)
sum.fit.dem <- summary(fit.dem)
coeff.fit.dem <- as.numeric(sum.fit.dem$parameters)

dem.n.vec <- seq(1,n.saf*2,1)
dem.pred <- ((coeff.fit.dem[1] * ((dem.n.vec)^coeff.fit.dem[2]))/((coeff.fit.dem[3]^coeff.fit.dem[2])+(dem.n.vec^coeff.fit.dem[2])))

## Plot demand function
row <- 1
col <- 1
par(mfrow=c(row,col))
plot(dem.n.vec,dem.pred,xlab="number of male trophy tags",ylab="demand probability",ylim=range(c(0,1)),type="l")
title(main="Trophy Demand Function")
par(mfrow=c(1,2))

## Translate demand to numbers of harvest animals
demf.pred <- (1 - (0.97^dem.n.vec))
demf <- demf.pred*n.saf

row <- 1
col <- 1
par(mfrow=c(row,col))
plot(dem.n.vec,demf.pred,xlab="number of male trophy tags",ylab="demand n probability",type="l")
title(main="Trophy n Demand Function")
par(mfrow=c(1,2))

## Cost function
cost.dif <- cost.max - cost.min

kprop.vec <- (seq(0.1,1,0.1))*Nd
cost.vec <- c(cost.max,0.833*cost.max,(cost.min+cost.max)/2,0.66*cost.max,1.25*cost.min,1.2*cost.min,1.15*cost.min,1.1*cost.min,1.05*cost.min,cost.min)

## Fit function (exponential decay) f=y0+a*exp(-b*x)
cost.data <- data.frame(kprop.vec,cost.vec)

fit.cost <- nls(cost.data$cost.vec ~ y0.cost + a.cost*exp(-b.cost*cost.data$kprop.vec),
		data = cost.data,
		start = list(a.cost = (0.66*cost.max), b.cost = 0.0004, y0.cost = cost.min),
		trace = TRUE)

sum.fit.cost <- summary(fit.cost)
coeff.fit.cost <- as.numeric(sum.fit.cost$parameters)

cost.n.vec <- seq(1,Nd,100)
cost.pred <- coeff.fit.cost[3] + coeff.fit.cost[1]*exp(-coeff.fit.cost[2]*cost.n.vec)

	for (c in 1:(length(cost.n.vec))) {
		if (cost.pred[c] > cost.max) cost.pred[c] <- cost.max
		if (cost.pred[c] < cost.min) cost.pred[c] <- cost.min
	}

## Plot cost function
row <- 1
col <- 1
par(mfrow=c(row,col))
plot(cost.n.vec,cost.pred,xlab="Population Size",ylab="Cost per head ($AU)",type="l")
title(main="Cost Function")
par(mfrow=c(1,2))


## Growth functions
f.wt.vec <- c(17,82,300); m.wt.vec <- c(17,100,600)
f.age.vec <- c(0,1,4); m.age.vec <- c(0,1,6)

## Create data frame
gr.data <- data.frame(f.wt.vec,f.age.vec,m.wt.vec,m.age.vec)

## model formula
## Von Bertalanffy growth function: M(t) = Mmax - (Mmax - M0) exp(-kt)
## M0 = mean birth mass
## Mmax = mean maximum mass
## k = rate constant per year

fit.gr.f <- nls(gr.data$f.wt.vec ~ gr.data$f.wt.vec[3] - (gr.data$f.wt.vec[3] - gr.data$f.wt.vec[1]) * exp(-k.coeff*gr.data$f.age.vec),
		data = gr.data,
		start = list(k.coeff = 0.2),
		trace = TRUE)
sum.fit.gr.f <- summary(fit.gr.f)

fit.gr.m <- nls(gr.data$m.wt.vec ~ gr.data$m.wt.vec[3] - (gr.data$m.wt.vec[3] - gr.data$m.wt.vec[1]) * exp(-k.coeff*gr.data$m.age.vec),
		data = gr.data,
		start = list(k.coeff = 0.2),
		trace = TRUE)
sum.fit.gr.m <- summary(fit.gr.m)

## Coefficients from fit
coeff.fit.gr.f <- as.numeric(sum.fit.gr.f$parameters)
coeff.fit.gr.m <- as.numeric(sum.fit.gr.m$parameters)

## Predict new q vector with integer age inputs
grf <- gr.data$f.wt.vec[3] - (gr.data$f.wt.vec[3] - gr.data$f.wt.vec[1]) * exp(-coeff.fit.gr.f[1]*age.vec)
grm <- gr.data$m.wt.vec[3] - (gr.data$m.wt.vec[3] - gr.data$m.wt.vec[1]) * exp(-coeff.fit.gr.m[1]*age.vec)

row <- 1
col <- 1
par(mfrow=c(row,col))
plot(age.vec,grf,xlab="age (years)",ylab="mass (kg)",ylim=range(c(0,max(m.wt.vec))),type="l")
lines(age.vec,grm,col="red")
title(main="von Bertalanffy Growth Functions",sub="(male vs. female)")
par(mfrow=c(1,2))

## Derivate growth expressions at "age.input"
vb.exp <- expression(Mmax - (Mmax - M0) * exp((-gr.k)*age))

## For a male
Mmax <- gr.data$m.wt.vec[3]
M0 <- gr.data$m.wt.vec[1]
gr.k <- coeff.fit.gr.m[1]
age1 <- 16/12
age2 <- age1 + (15/12)
age3 <- 10

age <- age1
dvb.exp1 <- D(vb.exp, c("age"))
m1.gr.yr <- eval(dvb.exp1)
m1.yr <- eval(vb.exp)
m1.gr.dy <- m1.gr.yr/365

age <- age2
dvb.exp2 <- D(vb.exp, c("age"))
m2.gr.yr <- eval(dvb.exp2)
m2.yr <- eval(vb.exp)
m2.gr.dy <- m2.gr.yr/365

age <- age3
dvb.exp3 <- D(vb.exp, c("age"))
m3.gr.yr <- eval(dvb.exp3)
m3.yr <- eval(vb.exp)
m3.gr.dy <- m3.gr.yr/365

## Average growth rate (per yr & per day) between
## the ages of 22 and 37 months (Moran 1973)
## Moran's value = 0.22 kg/day

m.sec.yr <- (m2.yr - m1.yr)/(age2 - age1)
m.sec.dy <- m.sec.yr/365 ## predicted kg/day

m.intcpt <- ((m1.yr + m2.yr) - ((age1+age2)*m.sec.yr))/2

row <- 1
col <- 1
par(mfrow=c(row,col))
plot(age.vec,grf,xlab="age (years)",ylab="mass (kg)",ylim=range(c(0,max(m.wt.vec))),type="l")
lines(age.vec,grm,col="red")
abline(v=age1,col="green")
abline(h=m1.yr,col="green")
abline(m.intcpt,m.sec.yr)
abline(v=age2,col="green")
abline(h=m2.yr,col="green")
##abline(v=age3,col="green")
##abline(h=m3.yr,col="green")
points(age1,m1.yr,col="black",pch=19)
points(age2,m2.yr,col="black",pch=19)
##points(age3,m3.yr,col="black",pch=19)
title(main="von Bertalanffy Growth Functions",sub="(male vs. female)")
par(mfrow=c(1,2))



## ************************************************************************************************
## Re-calculate logistic functions for survival & fecundity based on initial population vector

library(stats)
library(base)

pop.k <- Nd
pop.min <- 20
pop.step <- 50
xmid <- (Nd + pop.min)/2
num.vec <- seq(pop.min,pop.k,((pop.k-pop.min)/pop.step))

## Set up base logistic form
quant <- seq(0,1,(1/pop.step)) ## sets up quantile vector

## logistic.dens <- dlogis(quant, location = 0.01, scale = 10, log = FALSE) ## estimates logistic density function with long left tail
logistic.dens <- 0 + ((0.997)/(1+(quant/1.8217)^3.9151))
log.dens <- (logistic.dens - min(logistic.dens)) / max(logistic.dens - min(logistic.dens))

## Set min and max s0 values
s0.min <- 0.92
s0.max <- 0.62

## Calculate logistic function between min and max s0
s0.vec <- (seq(s0.min,s0.max,-((s0.min-s0.max)/pop.step)))
range.s0.vec <- range(s0.vec)[2] - range(s0.vec)[1]
s0.vec.logis <- (range.s0.vec * log.dens) + min(s0.vec)

## Set min and max s1+ values
s1p.min <- 0.99
s1p.max <- s0.min

## Calculate logistic function between min and max s1p
s1p.vec <- (seq(s1p.min,s1p.max,-((s1p.min-s1p.max)/pop.step)))
range.s1p.vec <- range(s1p.vec)[2] - range(s1p.vec)[1]
s1p.vec.logis <- (range.s1p.vec * log.dens) + min(s1p.vec)

## Set min and max fecundity a coefficient values
af.min <- 0.62
af.max <- coeff.fit.m.vec[1]

## Calculate logistic function between min and max
af.vec <- (seq(af.min,af.max,-((af.min-af.max)/pop.step)))
range.af.vec <- range(af.vec)[2] - range(af.vec)[1]
af.vec.logis <- (range.af.vec * log.dens) + min(af.vec)

## Plot dd functions
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(num.vec,s0.vec.logis,xlab="N",ylab="s0",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
title(main="Logistic DD Functions")
plot(num.vec,s1p.vec.logis,xlab="N",ylab="s1+",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec,af.vec.logis,xlab="N",ylab="a.f",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
par(mfrow=c(1,2))

log.func <- data.frame(num.vec,s0.vec.logis,s1p.vec.logis,af.vec.logis) ## put new vectors into data.frame

## Estimate 4-parameter logistic model coefficients
s0.param <- getInitial(s0.vec.logis ~ SSfpl(num.vec, a.s0, b.s0, xmid.s0, scal.s0), data=log.func)
s1p.param <- getInitial(s1p.vec.logis ~ SSfpl(num.vec, a.s1p, b.s1p, xmid.s1p, scal.s1p), data=log.func)
af.param <- getInitial(af.vec.logis ~ SSfpl(num.vec, a.af, b.af, xmid.af, scal.af), data=log.func)

## ************************************************************************************************

## Stochastic catastrophe scenario modelled from 0.14G (Reed et al. 2004)
p.cat <- 0.147/G.age.mean
sev.cat.vec <- seq(0,100,10)
sev.cat.prob <- (2510*(0.9376^(sev.cat.vec)))/100

sev.cat.dat <- data.frame(sev.cat.vec,sev.cat.prob)
sev.cat.fit <- nls(sev.cat.vec ~ (a.cat*(sev.cat.prob^b.cat)),
		data = sev.cat.dat,
		start = list(a.cat = 50, b.cat = -0.2),
		trace = TRUE)
sum.sev.cat.fit <- summary(sev.cat.fit)
sev.cat.coeff <- as.numeric(sum.sev.cat.fit$parameters)

karray <- seq(min.kill,max.kill,kill.interval) ## kill rate array
lkarray <- length(karray)

## Rainfall setup for stochastic projections
rain.sum.mean <- round(sum(c(275,240,300,140,20,5,2,0,5,20,120,230)),-1) ## mean monthly rainfall (1921-1978 - Choquenot 1993)
rain.var <- (rain.cv*rain.sum.mean)/100
rain.mult.sd <- (sqrt(rain.cv))/100
rain.mult.min <- 1-(2*rain.mult.sd)
rain.mult.max <- 1+(2*rain.mult.sd)

## Stochastic loop
## Set storage vectors for mean and confidence intervals for population sizes
min.popd.perm <- matrix(0,nrow=perm,ncol=lkarray)
minmal6p.perm <- matrix(0,nrow=perm,ncol=lkarray)
meat.dol.perm <- matrix(0,nrow=perm,ncol=lkarray)
cost.dol.perm <- matrix(0,nrow=perm,ncol=lkarray)
proft.dol.perm <- matrix(0,nrow=perm,ncol=lkarray)
profm.dol.perm <- matrix(0,nrow=perm,ncol=lkarray)
meat.mass.perm <- matrix(0,nrow=perm,ncol=lkarray)
safari.dol.perm <- matrix(0,nrow=perm,ncol=lkarray)

## Start stochastic loop

for (p in 1:perm) {

	min.popd.vec <- rep(0,length(karray)) ## Minimum population size vector
	mald.vec <- rep(0,length(karray)) ## Number of males vector
	minmal6p.vec <- rep(0,length(karray)) ## Minimum male (>5 years) vector
	minpmal6p.vec <- rep(0,length(karray)) ## Minimum proportion of males (>5 years) per males vector
	minfem3p.vec <- rep(0,length(karray)) ## Minimum female (>2 years) vector
	minpfem3p.vec <- rep(0,length(karray)) ## Minimum proportion of females (>2 years) per females vector
	lambdad.vec <- rep(0,length(karray)) ## mean population rate of change (r)
	mal6p.kill <- rep(0,length(karray)) ## mature males killed vector

	## Economic vectors
	meat.mass.tot <- rep(0,length(karray)) ## Total mass meat harvested
	meat.dol.tot <- rep(0,length(karray)) ## total AU$ of harvested meat
	safari.dol.tot <- rep(0,length(karray)) ## total AU$ of harvested meat
	cost.dol.tot <- rep(0,length(karray)) ## total AU$ costs
	meat.prop.worth <- rep(0,length(karray)) ## proportional worth of harvested meat
	safari.prop.worth <- rep(0,length(karray)) ## proporational worth of safari

	## Start simulation over kill vector

	for (z in 1:lkarray) {

	## Stochastic measure (variation in rainfall)
	rain.mult <- runif(1,min=rain.mult.min,max=rain.mult.max)

	## define start survival

	        s0.d <- rain.mult*(as.numeric((SSfpl(sum(Nd),s0.param[1],s0.param[2],s0.param[3],s0.param[4]))))
	            if (s0.d < s0.max) s0.d <- s0.max else s0.d <- s0.d
		    if (s0.d > s0.min) s0.d <- s0.min else s0.d <- s0.d

	        s1p.d <- rain.mult*(as.numeric((SSfpl(sum(Nd),s1p.param[1],s1p.param[2],s1p.param[3],s1p.param[4]))))
	            if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
		    if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d
        
	## define start fecundity
	        a.md <- as.numeric((SSfpl(sum(Nd),af.param[1],af.param[2],af.param[3],af.param[4])))
	            if (a.md < af.max) a.md <- af.max else a.md <- a.md
	            if (a.md > af.min) a.md <- af.min else a.md <- a.md
	            
	        ## Re-define fertility vector
	        md.vec <- rep(0,age.max); b.fert <- coeff.fit.m.vec[3]; x0.fert <- coeff.fit.m.vec[2]
	        for (j in (1:age.max)) {
	            md<-a.md*exp(-0.5*(log(j/x0.fert)/b.fert)^2)
	            md.vec[j] <- md
	        }

	##total females in population
	f<-Nd*sr
	
	##total males in population
	mal<-Nd*(1-sr)
	
	## The normal matrix
	ad <- matrix(data<-0,nrow<-k,ncol<-k)
	
	## Fill matrix with survival & fecundity vectors
	diag(ad[2:age.max,1:(age.max-1)]) <- s1p.d ## adds s1p.d to female quadrant
	diag(ad[(age.max+2):k,(age.max+1):(k-1)]) <- s1p.d ## adds s1p.d to male quadrant
	ad[1,2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds female fecundity
	ad[(age.max+1),2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds male fecundity

	## Stable stage distribution
	ssdd <- stable.stage.dist(ad)

	## pick initial vectors
	nd<-ssdd*Nd

	## max r
	rd <- max.r(ad)

	##set population size year step vector
	popd.vec <- rep(1,tlimit+1)
	popd.vec[1] <- sum(nd)

	##set female population size year step vector
	femd.vec<-rep(1,tlimit+1)
	femd.vec[1]<-sum(nd[1:age.max])

	femda.vec<-rep(1,tlimit+1)
	femda.vec[1]<-sum(nd[3:age.max])

	##set male population size year step vector
	mald.vec<-rep(1,tlimit+1)
	mald.vec[1]<-sum(nd[(age.max+1):k])

	malda.vec<-rep(1,tlimit+1)
	malda.vec[1]<-sum(nd[21:k])

	## set large male year step vector
	male6p.vec <- rep(0,tlimit+1)
	male6p.vec[1]<-sum(nd[23:k])

	##set year step vector
	yrd.vec <- rep(1,tlimit+1)
	yrd.vec[1] <- 0

	for (j in 1:tlimit) {
    		yrd.vec[j+1] <- j
	}

	## Mature male kill vector
	mal6p.kill.vec <- rep(0,tlimit+1)

	## Economic vectors
	meat.tot.vec <- rep(0,tlimit)
	meat.dol.vec <- rep(0,tlimit)
	cost.tot.vec <- rep(0,tlimit)
	safari.tot.vec <- rep(0,tlimit)
	meat.worth.vec <- rep(0,tlimit)
	safari.worth.vec <- rep(0,tlimit)
	meat.prop.vec <- rep(0,tlimit)
	safari.prop.vec <- rep(0,tlimit)

	eyci <- runif(1,eyci.min,eyci.max)

	##then iterate
		for (ii in 1:tlimit) {
			nd <- ad %*% nd

		## Calculate severity of a potential catastrophe
		sev.cat.input <- runif(1,0,(max(sev.cat.prob)))
		sev.cat <- ((sev.cat.coeff[1]*(sev.cat.input^sev.cat.coeff[2])))/100
		##1-sev.cat
		if (sev.cat > 1) sev.cat <- 0.99

		## calculate whether a catastrophe occurs
		if (runif(1,0,1) <= p.cat) nd <- nd*(1-sev.cat) else nd <- nd
		
		## Set any negative values in nd to 0
		zero.sub <- which(nd<0)
		nd[zero.sub] <- 0

		na.sub <- which(is.na(nd))
		nd[na.sub] <- 0

		if (sum(nd) == 0) next

		## Stochastic measure (variation in rainfall)
		rain.multp <- runif(1,min=rain.mult.min,max=rain.mult.max)
		## rain.multp <- ifelse(stoch.select == 1, rain.multp <- rain.multp, rain.mult <- 1)
		##rain.multp <- 1
		
	        ## Set negative density feedback function for the matrix
	
	        ## Redefine survival probabilities
	        s0.d <- rain.multp*(as.numeric((SSfpl(sum(nd),s0.param[1],s0.param[2],s0.param[3],s0.param[4]))))
	            if (s0.d < s0.max) s0.d <- s0.max else s0.d <- s0.d
		    if (s0.d > s0.min) s0.d <- s0.min else s0.d <- s0.d

	        s1p.d <- rain.multp*(as.numeric((SSfpl(sum(nd),s1p.param[1],s1p.param[2],s1p.param[3],s1p.param[4]))))
	            if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
		    if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d

		## re-define fecundity a coefficient
	        a.md <- as.numeric(SSfpl(sum(nd),af.param[1],af.param[2],af.param[3],af.param[4]))
	            if (a.md < af.max) a.md <- af.max else a.md <- a.md
	            if (a.md > af.min) a.md <- af.min else a.md <- a.md

	        ## Re-define fertility vector
	        md.vec <- rep(0,age.max)
		for (j in (1:age.max)) {
	            md <- a.md*exp(-0.5*(log(j/x0.fert)/b.fert)^2)
	            md.vec[j] <- md
	        }

		## Adjust for demographic stochasticity (binomial & poisson resampling)
		n.mf <- 2*md.vec[1:(age.max)]*nd[1:age.max]
		n.f <- sum(nd[1:age.max])
		sum.n.mf <- sum(n.mf)
		if (sum.n.mf == 0) sum.n.mf <- 1
		if ((is.na(sum.n.mf)) == TRUE) sum.n.mf <- 1
		mds.vec <- sum(md.vec)*(rpois(age.max,(n.mf*n.f))/(sum.n.mf*n.f))

			for (mm in 1:(length(md.vec))) {
				if (md.vec[mm] > af.min) md.vec[mm] <- af.min else md.vec[mm] <- md.vec[mm]
			} ## caps m at 0.62

	        ## Tertiary sex ratio
	        	if (sum(nd[1:age.max] / (1+sum(nd))) > tsr.max) md.vec <- rep(0,age.max) else md.vec <- md.vec

		n.1 <- round(sum(nd))
		if (n.1 == 0) n.1 <- 1
		s1p.d <- rbinom(1,n.1,s1p.d)/n.1		
		if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d

		n.0 <- round(sum(md.vec*n[1:age.max]))
		if (n.0 == 0) n.0 <- 1
		s0.d <- rbinom(1,n.0,s0.d)/n.0

		## Fill matrix with adjusted survival & fecundity vectors
		diag(ad[2:age.max,1:(age.max-1)]) <- s1p.d ## adds s1p.d to female quadrant
		diag(ad[(age.max+2):k,(age.max+1):(k-1)]) <- s1p.d ## adds s1p.d to male quadrant
		ad[1,2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds female fecundity
		ad[(age.max+1),2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds male fecundity
	
		## Fluctuating market values for beef sales (EYCI)
		## random uniform distribution between 1996-2004 min and max EYCI (in 2004 AU$)
		## Autocorrelated from last year
		eyci <- as.numeric(eyci.fit$coefficients[1] + (eyci.fit$coefficients[2])*eyci)
		eyci <- as.numeric(rnorm(1, mean=eyci,sd=(sqrt(eyci.var))))

		## if (stoch.select == 0) eyci <- eyci.now else eyci <- eyci
		##eyci <- ppkg.meat
		##print(eyci)


		## Cost function
		cost.now <- coeff.fit.cost[3] + coeff.fit.cost[1]*exp(-coeff.fit.cost[2]*(sum(nd)))
		if (cost.now > cost.max) cost.now <- cost.max
		if (cost.now < cost.min) cost.now <- cost.min
		##cost.now <- 0
	
	        ## *****************************************************************
	        ## Kill adult females (3-17 years)
	        ## *****************************************************************

		fem.kill <- ifelse(kill.type == 2, karray[z], 0)
	
	        nf.kill <- fem.kill/15 ## uniform distribution of kill in ages > 2
	        kill.f.vec <- rep(0,k)
	        kill.f.vec[3:age.max] <- nf.kill
	            if (sum(nd) < fem.kill) {
	                nd <- nd 
	            }      
	                nd <- nd - kill.f.vec ## Adjust nd vector for kill

		f.mass.tot <- sum(kill.f.vec[3:age.max] * grf[4:(age.max+1)])
		f.meat.tot <- f.mass.tot*meat.prop
		
		f.cost.tot <- sum(kill.f.vec[3:age.max]) * cost.now
		f.meat.dol <- f.meat.tot*eyci
	        
	        ## *****************************************************************

	        ## *****************************************************************
	        ## Kill large males (males 6-17 years)
	        ## *****************************************************************

		male.kill <- ifelse(kill.type == 1, karray[z], 0)

	        nm.kill <- male.kill/12 ## uniform distribution of kill in ages > 5
	        kill.m.vec <- rep(0,k)
	        kill.m.vec[23:k] <- nm.kill
	            if (sum(nd) < male.kill) {
	                nd <- nd 
	            }      
	                nd <- nd - kill.m.vec ## Adjust nd vector for kill

		m.mal6p.kill <- sum(kill.m.vec[23:k])
		m.mass.tot <- sum(kill.m.vec[23:k] * grm[7:(age.max+1)])

		m.cost.tot <- sum(kill.m.vec[23:k]) * cost.now
		m.meat.tot <- m.mass.tot*meat.prop
		m.meat.dol <- m.meat.tot*eyci

		m.m.kill <- m.mal6p.kill
		
		dem.sub <- ifelse((m.mal6p.kill > length(dem.pred)), length(dem.pred), ceiling(m.mal6p.kill))
		dem.kill <- (demf[dem.sub])
		
		m.safari.tot <- dem.kill * safari.to
		if (sum(kill.m.vec) == 0) m.safari.tot <- 0 else m.safari.tot <- m.safari.tot

	       
	        ## *****************************************************************

	        ## *****************************************************************
	        ## Kill adults (females >2 years; males >3)
	        ## *****************************************************************

		all.kill <- ifelse(kill.type == 0, karray[z], 0)

	        naf.kill <- (all.kill/2)/15 ## uniform distribution of females killed (> 2 years)
		nam.kill <- (all.kill/2)/14 ## uniform distribution of males killed (> 3 years)
	
	        kill.all.vec <- rep(0,k)
		kill.all.vec[3:age.max] <- naf.kill
		kill.all.vec[21:k] <- nam.kill
	            if (sum(nd) < all.kill) {
	                nd <- nd 
	            }      
	                nd <- nd - kill.all.vec ## Adjust nd vector for kill
	
		a.mal6p.kill <- sum(kill.all.vec[23:k])
		mass.tot.m6 <- sum(kill.all.vec[23:k] * grm[7:(age.max+1)])
		mass.tot.m45 <- sum(kill.all.vec[21:22] * grm[5:6])
		mass.tot.f <- sum(kill.all.vec[3:age.max] * grf[4:(age.max+1)])

		cost.tot.m <- sum(kill.all.vec[23:k]) * cost.now
		cost.tot.f <- sum(kill.all.vec[3:age.max]) * cost.now
		a.mass.tot <- mass.tot.m6 + mass.tot.m45 + mass.tot.f
		a.cost.tot <- cost.tot.m + cost.tot.f
		a.meat.tot <- a.mass.tot*meat.prop
		a.meat.dol <- a.meat.tot*eyci

		a.m.kill <- a.mal6p.kill
		dem.sub <- ifelse((a.mal6p.kill > length(dem.pred)), length(dem.pred), ceiling(a.mal6p.kill))
		dem.kill <- (demf[dem.sub])
		
		a.safari.tot <- dem.kill * safari.to
		if (sum(kill.all.vec) == 0) a.safari.tot <- 0 else a.safari.tot <- a.safari.tot
	       
	        ## *****************************************************************

		## Potential economic value of entire remaining herd
		fem.meat.worth <- sum(nd[1:age.max]*grf[2:(age.max+1)]*meat.prop*eyci)
		mal.meat.worth <- sum(nd[(age.max+1):k]*grm[2:(age.max+1)]*meat.prop*eyci)
		mal.saf.worth <- sum(nd[23:k])*safari.to

		     popd.vec[ii+1]<-(sum(nd))    
		     femd.vec[ii+1]<-sum(nd[1:age.max]) ## all females
		     femda.vec[ii+1]<-sum(nd[3:age.max]) ## adult (>2 years) females
		     mald.vec[ii+1]<-sum(nd[(age.max+1):k]) ## all males
		     malda.vec[ii+1]<-sum(nd[21:k]) ## adult (>3 years) males
		     male6p.vec[ii+1]<-sum(nd[23:k]) ## males >5 years
		     mal6p.kill.vec[ii+1]<-a.mal6p.kill + m.mal6p.kill

		     meat.tot.vec[ii] <- f.meat.tot + m.meat.tot + a.meat.tot
		     meat.dol.vec[ii] <- f.meat.dol + m.meat.dol + a.meat.dol
		     safari.tot.vec[ii] <- m.safari.tot + a.safari.tot
		     cost.tot.vec[ii] <- f.cost.tot + m.cost.tot + a.cost.tot
		     meat.worth.vec[ii] <- fem.meat.worth + mal.meat.worth
		     safari.worth.vec[ii] <- mal.saf.worth
		     meat.prop.vec[ii] <- (f.meat.tot + m.meat.tot + a.meat.tot)/(fem.meat.worth + mal.meat.worth)
		     safari.prop.vec[ii] <- (m.safari.tot + a.safari.tot)/mal.saf.worth
	     
		} ## end ii loop

	log.popd.vec<-log10(popd.vec)

	## Place data in z vectors
	min.popd.vec[z] <- min(ifelse(popd.vec>0,popd.vec,0))

	minmal6p.vec[z] <- min(male6p.vec)
	min.male6p.sub <- which(male6p.vec == min(male6p.vec))
	min.mald <- mald.vec[min.male6p.sub]
	minpmal6p.vec[z] <- min(male6p.vec) / min.mald[1]

	minfem3p.vec[z] <- min(femda.vec)
	min.female3p.sub <- which(femda.vec == min(femda.vec))
	min.femd <- femd.vec[min.female3p.sub]
	minpfem3p.vec[z] <- min(femda.vec) / min.femd[1]

	mal6p.kill[z] <- sum(mal6p.kill.vec) / tlimit ## average annual mature male kill

	## Economic vectors
	meat.mass.tot[z] <- sum(meat.tot.vec) / tlimit ## average per annum value
	meat.dol.tot[z] <- sum(meat.dol.vec) / tlimit ## average per annum value
	safari.dol.tot[z] <- sum(safari.tot.vec) / tlimit ## average per annum value
	cost.dol.tot[z] <- sum(cost.tot.vec) / tlimit
	meat.prop.worth[z] <- sum(meat.prop.vec) / tlimit ## average proportional worth of meat
	safari.prop.worth[z] <- sum(safari.prop.vec) / tlimit ## average proportional worth of meat

	## Calculate stochastic r
	popd1.vec<-rep(0,length(popd.vec)+1)
	popd1.vec[2:(length(popd.vec)+1)] <- popd.vec
	lambda.m <- mean((popd.vec[2:(length(popd.vec))]) / (popd1.vec[2:(length(popd.vec))]))
	lambdad.vec[z] <- lambda.m

	if (stoch.select == 1) next

	## Single projection plots
	row <- 2
	col <- 2
	par(mfrow=c(row,col))
	plot(yrd.vec,popd.vec,xlab="year",ylab="N",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(popd.vec)+(0.25*(max(popd.vec))))),type='l')
	title(main = "Kill rate")	
	mtext("=", line=0, at=8)
	mtext(z, line=0, at=14)
	mtext("per annum", line=0, at=25)
	plot(yrd.vec,femda.vec,xlab="year",ylab="N adult females",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(femda.vec)+(0.25*(max(femda.vec))))),type='l')
	title(main=scenario)
	plot(yrd.vec,malda.vec,xlab="year",ylab="N adult males",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(malda.vec)+(0.25*(max(mald.vec))))),type='l')
	plot(yrd.vec,male6p.vec,xlab="year",ylab="N males >5 years",xlim=range(-0.5:(tlimit+1)),ylim=range(0:(max(male6p.vec)+(0.25*(max(male6p.vec))))),type='l')
	par(mfrow=c(1,2))

	} ## end z loop


	## Final economics post-z loop
	profm.dol.tot <- meat.dol.tot - cost.dol.tot
	proft.dol.tot <- (meat.dol.tot + safari.dol.tot) - cost.dol.tot

	## Fix min.popd.vec
	reals1 <- which(min.popd.vec>0)
	reals1.sub <- which(diff(reals1)>1)
	reals1a.sub <- ifelse((is.na(summary(reals1.sub)[1] > 0) == TRUE), length(min.popd.vec), reals1.sub)
	    if (reals1a.sub[1] != length(min.popd.vec)) {
		            min.popd.vec[((reals1a.sub[1]+1)):(length(min.popd.vec))] <- 0
	        }
	            min.popd.vec <- min.popd.vec

	## Fix minmal6p
	reals2 <- which(minmal6p.vec>0)
	reals2.sub <- which(diff(reals2)>1)
	reals2a.sub <- ifelse((is.na(summary(reals2.sub)[1] > 0) == TRUE), length(minmal6p.vec), reals2.sub)
	    if (reals2a.sub[1] != length(minmal6p.vec)) {
	            minmal6p.vec[((reals2a.sub[1]+1)):(length(minmal6p.vec))] <- 0
	        }
	            minmal6p.vec <- minmal6p.vec

	## Fix minfem3p
	reals3 <- which(minfem3p.vec>0)
	reals3.sub <- which(diff(reals3)>1)
	reals3a.sub <- ifelse((is.na(summary(reals3.sub)[1] > 0) == TRUE), length(minfem3p.vec), reals3.sub)
	    if (reals3a.sub[1] != length(minfem3p.vec)) {
	            minfem3p.vec[((reals3a.sub[1]+1)):(length(minfem3p.vec))] <- 0
	        }
	            minfem3p.vec <- minfem3p.vec

		## Set any negative values to 0
		zero.sub3 <- which(minfem3p.vec<0)
		minfem3p.vec[zero.sub3] <- 0


	## Fix minpfem3p
	reals4 <- which(minpfem3p.vec>0)
	reals4.sub <- which(diff(reals4)>1)
	reals4a.sub <- ifelse((is.na(summary(reals4.sub)[1] > 0) == TRUE), length(minpfem3p.vec), reals4.sub)
	    if (reals4a.sub[1] != length(minpfem3p.vec)) {
 	           minpfem3p.vec[((reals4a.sub[1]+1)):(length(minpfem3p.vec))] <- 0
 	       }
 	           minpfem3p.vec <- minpfem3p.vec

	## Fix minpmal6p
	reals5 <- which(minpmal6p.vec>0)
	reals5.sub <- which(diff(reals5)>1)
	reals5a.sub <- ifelse((is.na(summary(reals5.sub)[1] > 0) == TRUE), length(minpmal6p.vec), reals5.sub)
	    if (reals5a.sub[1] != length(minpmal6p.vec)) {
	            minpmal6p.vec[((reals5a.sub[1]+1)):(length(minpmal6p.vec))] <- 0
	        }
	            minpmal6p.vec <- minpmal6p.vec

	## Fix stochastic r for real numbers only
	reals6 <- which(lambdad.vec>0)
	reals6.sub <- which(diff(reals6)>1)
	reals6a.sub <- ifelse((is.na(summary(reals6.sub)[1] > 0) == TRUE), length(lambdad.vec), reals6.sub)
	    if (reals6a.sub[1] != length(lambdad.vec)) {
	            lambdad.vec[((reals6a.sub[1]+1)):(length(lambdad.vec))] <- 0
	        }
 	           lambdad.vec <- lambdad.vec

	## If population goes to 0, make lambdad.vec 0
	zero.popd.sub <- which(min.popd.vec==0)
	zero.popd.sub1 <- zero.popd.sub[1]
	if (is.na(zero.popd.sub1) == "TRUE") zero.popd.sub1 <- 0 else zero.popd.sub1 <- zero.popd.sub1
	if (zero.popd.sub1 > 0) lambdad.vec[zero.popd.sub1:(length(lambdad.vec))] <- 0 else lambdad.vec <- lambdad.vec
	
	## If population goes to 0, proportion male/female vectors go to 0
	if ((zero.popd.sub1 > 0) & (minpmal6p.vec[zero.popd.sub1-1] == 1)) minpmal6p.vec[zero.popd.sub1:(length(minpmal6p.vec))] <- 0 else minpmal6p.vec <- minpmal6p.vec
	if ((zero.popd.sub1 > 0) & (minpfem3p.vec[zero.popd.sub1-1] == 1)) minpfem3p.vec[zero.popd.sub1:(length(minpfem3p.vec))] <- 0 else minpfem3p.vec <- minpfem3p.vec

	## If population goes to 0, make profit vectors go to 0
	if (zero.popd.sub1 > 0) meat.dol.tot[zero.popd.sub1:(length(meat.dol.tot))] <- 0 else meat.dol.tot <- meat.dol.tot
	if (zero.popd.sub1 > 0) meat.mass.tot[zero.popd.sub1:(length(meat.mass.tot))] <- 0 else meat.mass.tot <- meat.mass.tot
##	if (zero.popd.sub1 > 0) cost.dol.tot[zero.popd.sub1:(length(cost.dol.tot))] <- 0 else cost.dol.tot <- cost.dol.tot
	if (zero.popd.sub1 > 0) proft.dol.tot[zero.popd.sub1:(length(proft.dol.tot))] <- 0 else proft.dol.tot <- proft.dol.tot
	if (zero.popd.sub1 > 0) profm.dol.tot[zero.popd.sub1:(length(profm.dol.tot))] <- 0 else profm.dol.tot <- profm.dol.tot


	## Find kill rate required to reduce minimum population by x %
	min.red.pc <- round(((Nd - min.popd.vec)/Nd)*100)
	thresh.sub <- (which(min.red.pc == thresh))[1]
	thresh.kill <- karray[thresh.sub]

	print(p)

	## storage vectors for stochastic simulation
	min.popd.perm[p,] <- min.popd.vec; minmal6p.perm[p,] <- minmal6p.vec
	meat.dol.perm[p,] <- meat.dol.tot; proft.dol.perm[p,] <- proft.dol.tot
	meat.mass.perm[p,] <- meat.mass.tot; cost.dol.perm[p,] <- cost.dol.tot
	safari.dol.perm[p,] <- safari.dol.tot; profm.dol.perm[p,] <- profm.dol.tot

	if (stoch.select == 1) next

	## Make density-dependent plots
	row <- 3
	col <- 2
	par(mfrow=c(row,col))

	plot(karray,minmal6p.vec,xlab="animals killed/year",ylab="Min N mature males",xlim=range(-0.5:(max(karray))),ylim=range(0:(max(minmal6p.vec)+(0.25*(max(minmal6p.vec))))),type='l')
	abline(v=thresh.kill,col="red")
	title(main = "Scenario = ")
	plot(karray,minpmal6p.vec,xlab="animals killed/year",ylab="Min prop mature males",xlim=range(-0.5:(max(karray))),ylim=range(c(0,1)),type='l')
	title(main = scenario)
	abline(v=thresh.kill,col="red")

	plot(karray,minfem3p.vec,xlab="animals killed/year",ylab="Min N adult females",xlim=range(-0.5:(max(karray))),ylim=range(0:(max(minfem3p.vec)+(0.25*(max(minfem3p.vec))))),type='l')
	abline(v=thresh.kill,col="red")
	plot(karray,minpfem3p.vec,xlab="animals killed/year",ylab="Min prop adult females",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(minpfem3p.vec)+(0.25*(max(minpfem3p.vec))))),type='l')
	abline(v=thresh.kill,col="red")

	plot(karray,min.popd.vec,xlab="animals killed/year",ylab="Min N",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(min.popd.vec)+(0.25*(max(min.popd.vec))))),type='l')
	abline(v=thresh.kill,col="red")
	plot(karray,lambdad.vec,xlab="animals killed/year",ylab="stochastic lambda",xlim=range(-0.5:(max(karray))),ylim=c(0.95,1.05),type='l')
	lines(karray,(rep(1,length(karray))),type='l',col="red")
	abline(v=thresh.kill,col="red")

	par(mfrow=c(1,2))

	## Summary
	scenario ## kill scenario
	z ## final kill rate per annum
	Nd ## initial population size
	min.popd.vec[z] ## lowest population size
	((Nd - min.popd.vec[z])/Nd)*100 ## per cent reduction at final kill rate
	thresh ## target per cent reduction of original population size
	thresh.kill ## kill rate required to reduce population by x %
	lambdad.vec[z] ## stochastic lambda

	## Economic plots
	## $ Returns based on harvest
	row <- 2
	col <- 2
	par(mfrow=c(row,col))
	plot(karray,meat.dol.tot,xlab="animals killed/year",ylab="per annum meat value",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(meat.dol.tot)+(0.25*(max(meat.dol.tot))))))),type='l',col="red")
	title(main = "Scenario = ")
	lines(karray,(meat.dol.tot*exch.aud.euro),col="black")
	plot(mal6p.kill,safari.dol.tot,xlab="mature males killed/year",ylab="per annum TO safari return",xlim=range(-0.5:(max(mal6p.kill))),ylim=range(c(0,((max(safari.dol.tot)+(0.25*(max(safari.dol.tot))))))),type='l',col="red")
	title(main = scenario)
	lines(mal6p.kill,(safari.dol.tot*exch.aud.euro),col="black")
	plot(karray,(meat.dol.tot+safari.dol.tot),xlab="animals killed/year",ylab="total per annum value",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max((meat.dol.tot+safari.dol.tot))+(0.25*(max((meat.dol.tot+safari.dol.tot)))))))),type='l',col="red")
	lines(karray,((meat.dol.tot+safari.dol.tot)*exch.aud.euro),col="red")
	abline(h=(max(meat.dol.tot+safari.dol.tot)),type="l",col="black")
	mtext("Max return per annum (AU$) =", outer=TRUE, side=3, line = -20, adj=0.05, cex = 1.1, col="red")
	mtext(round(max(meat.dol.tot+safari.dol.tot)), outer=TRUE, side=3, line = -20, adj=0.5, cex = 1)
	plot(karray,min.popd.vec,xlab="animals killed/year",ylab="Min N",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(min.popd.vec)+(0.25*(max(min.popd.vec))))),type='l',col="black")
	abline(v=thresh.kill,col="red")
	par(mfrow=c(1,1))

	## Profit
	row <- 2
	col <- 2
	par(mfrow=c(row,col))
##	plot(karray,meat.prop.worth,xlab="animals killed/year",ylab="avg prop value",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(safari.prop.worth)+(0.25*(max(safari.prop.worth))))))),type='l',col="black")
##	lines(mal6p.kill,safari.prop.worth,col="red")
	plot(karray,(meat.mass.tot),xlab="animals killed/year",ylab="harvested meat (kg)",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(meat.mass.tot)+(0.25*(max(meat.mass.tot))))))),type='l',col="red")
##	abline(h=minjil.meat,col="black")
##	title(main = "Scenario = ")
##	lines(karray,((meat.dol.tot-cost.dol.tot)*exch.aud.euro),col="black")
##	mtext("Max mean safari profit per annum (AU$) =", outer=TRUE, side=3, line = -25, adj=0.05, cex = 1.1, col="red")
##	mtext(round(max((prof.dol.tot))), outer=TRUE, side=3, line = -25, adj=0.60, cex = 1)
##	mtext("Max mean meat profit per annum (AU$) =", outer=TRUE, side=3, line = -20, adj=0.05, cex = 1.1, col="red")
##	mtext(round(max((meat.dol.tot-cost.dol.tot))), outer=TRUE, side=3, line = -20, adj=0.65, cex = 1)
	plot(karray,(meat.dol.tot-cost.dol.tot),xlab="animals killed/year",ylab="total meat profit",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(meat.dol.tot-cost.dol.tot)+(0.25*(max(meat.dol.tot-cost.dol.tot))))))),type='l',col="red")
##	abline(h=minjil.meat*ppkg.meat,col="black")
##	title(main = scenario)
##	lines(karray,(prof.dol.tot*exch.aud.euro),col="black")
	plot(karray,(safari.dol.tot),xlab="animals killed/year",ylab="total safari profit",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(safari.dol.tot)+(0.25*(max(safari.dol.tot))))))),type='l',col="red")
	plot(karray,(prof.dol.tot),xlab="animals killed/year",ylab="total profit",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(prof.dol.tot)+(0.25*(max(prof.dol.tot))))))),type='l',col="red")
	plot(karray,min.popd.vec,xlab="animals killed/year",ylab="Min N",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(min.popd.vec)+(0.25*(max(min.popd.vec))))),type='l',col="black")
##	abline(v=thresh.kill,col="red")
	abline(h=min(min.popd.vec),col="red")
##	plot(karray,minmal6p.vec,xlab="animals killed/year",ylab="Min N mature males",xlim=range(-0.5:(max(karray))),ylim=range(0:(max(minmal6p.vec)+(0.25*(max(minmal6p.vec))))),type='l')
##	abline(v=thresh.kill,col="red")
##	abline(h=min(minmal6p.vec),col="red")
##	plot(karray,lambdad.vec,xlab="animals killed/year",ylab="stochastic lambda",xlim=range(-0.5:(max(karray))),ylim=c(0.95,1.05),type='l')
##	lines(karray,(rep(1,length(karray))),type='l',col="red")
##	abline(v=thresh.kill,col="red")

	par(mfrow=c(1,1))

	####################
	## Economic summary
	####################

	## Kill max
	z
	## Years over which average calculated
	tlimit
	## Costs
	max(cost.dol.tot)
	## EYCI
	eyci.now ## todays price per kg on beef market (AU$)
	## Returns (meat)
	max(meat.dol.tot)
	## Returns (safari)
	max(safari.dol.tot)
	## Returns (total)
	max(meat.dol.tot+safari.dol.tot)
	## Profit (meat)
	max(meat.dol.tot-cost.dol.tot)
	## Profit (meat) based on Australian average
	(1+(profit/100))*(max(round(cost.dol.tot))) - max(round(cost.dol.tot))
	## Profit (total)
	max(prof.dol.tot)
	## % return on investment (total)
	((round(max(meat.dol.tot+safari.dol.tot)) / round(max(cost.dol.tot)))-1)*100
	## % return on investment (meat)
	(((round(max(meat.dol.tot)) - round(max(cost.dol.tot))) / round(max(cost.dol.tot))))*100

} ## end p loop (stochastic)




## Mean & confidence intervals over p permutations
min.popd.mean <- rep(0,lkarray); min.popd.up <- rep(0,lkarray); min.popd.lo <- rep(0,lkarray)
minmal6p.mean <- rep(0,lkarray); minmal6p.up <- rep(0,lkarray); minmal6p.lo <- rep(0,lkarray)
meat.dol.mean <- rep(0,lkarray); meat.dol.up <- rep(0,lkarray); meat.dol.lo <- rep(0,lkarray)
proft.dol.mean <- rep(0,lkarray); proft.dol.up <- rep(0,lkarray); proft.dol.lo <- rep(0,lkarray)
profm.dol.mean <- rep(0,lkarray); profm.dol.up <- rep(0,lkarray); profm.dol.lo <- rep(0,lkarray)
meat.mass.mean <- rep(0,lkarray); meat.mass.up <- rep(0,lkarray); meat.mass.lo <- rep(0,lkarray)
cost.dol.mean <- rep(0,lkarray); cost.dol.up <- rep(0,lkarray); cost.dol.lo <- rep(0,lkarray)
safari.dol.mean <- rep(0,lkarray); safari.dol.up <- rep(0,lkarray); safari.dol.lo <- rep(0,lkarray)

for (d in 1:(lkarray)) {
	min.popd.mean[d] <- mean(min.popd.perm[,d])
	min.popd.up[d] <- quantile(min.popd.perm[,d],probs=0.975)
	min.popd.lo[d] <- quantile(min.popd.perm[,d],probs=0.025)

	minmal6p.mean[d] <- mean(minmal6p.perm[,d])
	minmal6p.up[d] <- quantile(minmal6p.perm[,d],probs=0.975)
	minmal6p.lo[d] <- quantile(minmal6p.perm[,d],probs=0.025)

	meat.dol.mean[d] <- mean(meat.dol.perm[,d])
	meat.dol.up[d] <- quantile(meat.dol.perm[,d],probs=0.975)
	meat.dol.lo[d] <- quantile(meat.dol.perm[,d],probs=0.025)

	proft.dol.mean[d] <- mean(proft.dol.perm[,d])
	proft.dol.up[d] <- quantile(proft.dol.perm[,d],probs=0.975)
	proft.dol.lo[d] <- quantile(proft.dol.perm[,d],probs=0.025)

	profm.dol.mean[d] <- mean(profm.dol.perm[,d])
	profm.dol.up[d] <- quantile(profm.dol.perm[,d],probs=0.975)
	profm.dol.lo[d] <- quantile(profm.dol.perm[,d],probs=0.025)

	meat.mass.mean[d] <- mean(meat.mass.perm[,d])
	meat.mass.up[d] <- quantile(meat.mass.perm[,d],probs=0.975)
	meat.mass.lo[d] <- quantile(meat.mass.perm[,d],probs=0.025)

	cost.dol.mean[d] <- mean(cost.dol.perm[,d])
	cost.dol.up[d] <- quantile(cost.dol.perm[,d],probs=0.975)
	cost.dol.lo[d] <- quantile(cost.dol.perm[,d],probs=0.025)

	safari.dol.mean[d] <- mean(safari.dol.perm[,d])
	safari.dol.up[d] <- quantile(safari.dol.perm[,d],probs=0.975)
	safari.dol.lo[d] <- quantile(safari.dol.perm[,d],probs=0.025)

} ## end d loop

if (sum(min.popd.up) == 0) min.popd.up <- min.popd.mean else min.popd.up <- min.popd.up
if (sum(min.popd.lo) == 0) min.popd.lo <- min.popd.mean else min.popd.lo <- min.popd.lo
if (sum(minmal6p.up) == 0) minmal6p.up <- minmal6p.mean else minmal6p.up <- minmal6p.up
if (sum(minmal6p.lo) == 0) minmal6p.lo <- minmal6p.mean else minmal6p.lo <- minmal6p.lo
if (sum(meat.dol.up) == 0) meat.dol.up <- meat.dol.mean else meat.dol.up <- meat.dol.up
if (sum(meat.dol.lo) == 0) meat.dol.lo <- meat.dol.mean else meat.dol.lo <- meat.dol.lo
if (sum(proft.dol.up) == 0) proft.dol.up <- proft.dol.mean else proft.dol.up <- proft.dol.up
if (sum(proft.dol.lo) == 0) proft.dol.lo <- proft.dol.mean else proft.dol.lo <- proft.dol.lo
if (sum(profm.dol.up) == 0) profm.dol.up <- profm.dol.mean else profm.dol.up <- profm.dol.up
if (sum(profm.dol.lo) == 0) profm.dol.lo <- profm.dol.mean else profm.dol.lo <- profm.dol.lo
if (sum(meat.mass.up) == 0) meat.mass.up <- meat.mass.mean else meat.mass.up <- meat.mass.up
if (sum(meat.mass.lo) == 0) meat.mass.lo <- meat.mass.mean else meat.mass.lo <- meat.mass.lo
if (sum(cost.dol.up) == 0) cost.dol.up <- cost.dol.mean else cost.dol.up <- cost.dol.up
if (sum(cost.dol.lo) == 0) cost.dol.lo <- cost.dol.mean else cost.dol.lo <- cost.dol.lo
if (sum(safari.dol.up) == 0) safari.dol.up <- safari.dol.mean else safari.dol.up <- safari.dol.up
if (sum(safari.dol.lo) == 0) safari.dol.lo <- safari.dol.mean else safari.dol.lo <- safari.dol.lo

## Make stochastic density-dependent plots
row <- 2
col <- 2
par(mfrow=c(row,col))
library(plotrix)

#plot(karray,meat.mass.mean,type="l",xlab="animals killed/year",ylab="per annum meat harvested (kg)",)
#lines(karray,meat.mass.lo,col="red")
#lines(karray,meat.mass.up,col="red")
#abline(h=minjil.meat,col="black")
#boxed.labels(3,14500,"a",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

#plot(karray,ppkg.meat*meat.mass.mean,type="l",xlab="animals killed/year",ylab="replacement value (AU$)",)
##lines(karray,meat.mass.lo,col="red")
##lines(karray,meat.mass.up,col="red")
#abline(h=minjil.meat*ppkg.meat,col="black")
#boxed.labels(3,145000,"b",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

plot(karray,profm.dol.mean,xlab="animals killed/year",ylab="per annum meat profit (AU$)",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(proft.dol.up)+(0.25*(max(proft.dol.up))))))),type='l',col="blue")
##title(main = "Scenario = ")
lines(karray,profm.dol.up,col="red")
lines(karray,profm.dol.lo,col="red")
##abline(h=minjil.meat*ppkg.meat,col="black")
boxed.labels(40,8.1e+05,"a",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

plot(karray,safari.dol.mean,xlab="animals killed/year",ylab="per annum safari revenue (AU$)",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(safari.dol.up)+(0.25*(max(safari.dol.up))))))),type='l',col="blue")
##title(main = "Scenario = ")
lines(karray,safari.dol.up,col="red")
lines(karray,safari.dol.lo,col="red")
boxed.labels(40,110000,"b",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

plot(karray,proft.dol.mean,xlab="animals killed/year",ylab="per annum total profit (AU$)",xlim=range(-0.5:(max(karray))),ylim=range(c(0,((max(proft.dol.up)+(0.25*(max(proft.dol.up))))))),type='l',col="blue")
##title(main = "Scenario = ")
lines(karray,proft.dol.up,col="red")
lines(karray,proft.dol.lo,col="red")
boxed.labels(40,8.1e+05,"c",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

plot(karray,min.popd.mean,xlab="animals killed/year",ylab="Min N",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(min.popd.up)+(0.25*(max(min.popd.up))))),type='l')
lines(karray,min.popd.up,col="red")
lines(karray,min.popd.lo,col="red")
boxed.labels(40,10500,"d",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

#plot(karray,min.popd.mean,xlab="mature males killed/year",ylab="Min N",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(min.popd.up)+(0.25*(max(min.popd.up))))),type='l')
#lines(karray,min.popd.up,col="red")
#lines(karray,min.popd.lo,col="red")
#boxed.labels(45,11000,"b",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

#plot(karray,minmal6p.mean,xlab="mature males killed/year",ylab="Min N mature male",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(minmal6p.up)+(0.25*(max(minmal6p.up))))),type='l')
#lines(karray,minmal6p.up,col="red")
#lines(karray,minmal6p.lo,col="red")
#abline(v=362)
#boxed.labels(45,3000,"c",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)

par(mfrow=c(1,1))





#############################################
## Fit logistic curve to min.popd.vec values
## yhat = a / ( 1 + (x/x0)^b)
#############################################

## Create data.frame
min.popdf.data <- data.frame(karray,min.popd.vec)
min.popd.data <- data.frame(karray,min.popd.vec)

fit.min.popdf <- nls(min.popdf.data$min.popd.vec ~ a.mp / (1 + (min.popdf.data$karray/x0.mp)^b.mp),
		data = min.popdf.data,
		start = list(a.mp = 8672, b.mp = 2, x0.mp = 1020),
		trace = TRUE)


fit.min.popd <- nls(min.popd.data$min.popd.vec ~ a.mp / (1 + (min.popd.data$karray/x0.mp)^b.mp),
		data = min.popdm.data,
		start = list(a.mp = 8672, b.mp = 2, x0.mp = 1020),
		trace = TRUE)

sum.fit.mpf <- summary(fit.min.popdf)
sum.fit.mp <- summary(fit.min.popd)

coeff.fit.mpf <- as.numeric(sum.fit.mpf$parameters)
coeff.fit.mp <- as.numeric(sum.fit.mp$parameters)

mpf.pred <- coeff.fit.mpf[1] / (1 + (karray/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mp.pred <- coeff.fit.mp[1] / (1 + (karray/coeff.fit.mp[3])^coeff.fit.mp[2])

mpf.pred.90 <- coeff.fit.mpf[1] / (1 + ((0.90*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.80 <- coeff.fit.mpf[1] / (1 + ((0.80*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.70 <- coeff.fit.mpf[1] / (1 + ((0.70*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.60 <- coeff.fit.mpf[1] / (1 + ((0.60*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.50 <- coeff.fit.mpf[1] / (1 + ((0.50*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.40 <- coeff.fit.mpf[1] / (1 + ((0.40*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.30 <- coeff.fit.mpf[1] / (1 + ((0.30*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.20 <- coeff.fit.mpf[1] / (1 + ((0.20*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])
mpf.pred.10 <- coeff.fit.mpf[1] / (1 + ((0.10*Nd)/coeff.fit.mpf[3])^coeff.fit.mpf[2])

mp.pred.90 <- coeff.fit.mp[1] / (1 + ((0.90*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.80 <- coeff.fit.mp[1] / (1 + ((0.80*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.70 <- coeff.fit.mp[1] / (1 + ((0.70*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.60 <- coeff.fit.mp[1] / (1 + ((0.60*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.50 <- coeff.fit.mp[1] / (1 + ((0.50*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.40 <- coeff.fit.mp[1] / (1 + ((0.40*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.30 <- coeff.fit.mp[1] / (1 + ((0.30*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.20 <- coeff.fit.mp[1] / (1 + ((0.20*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])
mp.pred.10 <- coeff.fit.mp[1] / (1 + ((0.10*Nd)/coeff.fit.mp[3])^coeff.fit.mp[2])

prop.Nd <- c(0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0)
mpf.pred.vec <- c(mpf.pred.90,mpf.pred.80,mpf.pred.70,mpf.pred.60,mpf.pred.50,mpf.pred.40,mpf.pred.30,mpf.pred.20,mpf.pred.10)
mp.pred.vec <- c(mp.pred.90,mp.pred.80,mp.pred.70,mp.pred.60,mp.pred.50,mp.pred.40,mp.pred.30,mp.pred.20,mp.pred.10)

## Plot function
row <- 1
col <- 1
par(mfrow=c(row,col))
plot(karray,min.popd.vec,xlab="animals killed/year",ylab="Min N",xlim=range(-0.5:(max(karray))),ylim=c(0,(max(min.popd.vec)+(0.25*(max(min.popd.vec))))),type='l',col="black")
lines(karray,mp.pred,col="red")
par(mfrow=c(1,2))

row <- 1
col <- 1
par(mfrow=c(row,col))
plot(mp.pred.vec,prop.Nd[1:9],ylab="proportion of K",xlab="number of animals killed per annum",ylim=range(c(1,0)),xlim=c(0,2500),type='l',col="black")
##lines(mp.pred.vec,prop.Nd[1:9],col="red")
par(mfrow=c(1,2))




######################
## Stochastic analysis
## no kill
######################

## Re-define number of permutations
perm <- 1000

## Re-define start population size
Nd <- 20

## Re-define tlimit
tlimit <- 150

## Re-define yr.vec
yr.vec <- rep(1,tlimit+1)
yr.vec[1] <- 0

	for (j in 1:tlimit) {
    		yr.vec[j+1] <- j
	}

## Stochastic loop
## Set storage vectors for mean and confidence intervals for population sizes
popd.sperm <- matrix(0,nrow=perm,ncol=tlimit+1)
femd.sperm <- matrix(0,nrow=perm,ncol=tlimit+1)
mald.sperm <- matrix(0,nrow=perm,ncol=tlimit+1)
male6p.sperm <- matrix(0,nrow=perm,ncol=tlimit+1)

## Start stochastic loop

for (p in 1:perm) {

	## define start survival

	        s0.d <- as.numeric(SSfpl(sum(Nd),s0.param[1],s0.param[2],s0.param[3],s0.param[4]))
	            if (s0.d < s0.max) s0.d <- s0.max else s0.d <- s0.d
		    if (s0.d > s0.min) s0.d <- s0.min else s0.d <- s0.d

	        s1p.d <- as.numeric(SSfpl(sum(Nd),s1p.param[1],s1p.param[2],s1p.param[3],s1p.param[4]))
	            if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
		    if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d
        
	## define start fecundity
	        a.md <- as.numeric(SSfpl(sum(Nd),af.param[1],af.param[2],af.param[3],af.param[4]))
	            if (a.md < af.max) a.md <- af.max else a.md <- a.md
	            if (a.md > af.min) a.md <- af.min else a.md <- a.md
	            
	        ## Re-define fertility vector
	        md.vec <- rep(0,age.max); b.fert <- coeff.fit.m.vec[3]; x0.fert <- coeff.fit.m.vec[2]
	        for (j in (1:age.max)) {
	            md<-a.md*exp(-0.5*(log(j/x0.fert)/b.fert)^2)
	            md.vec[j] <- md
	        }
	##total females in population
	f<-Nd*sr
	
	##total males in population
	mal<-Nd*(1-sr)
	
	## The normal matrix
	ad <- matrix(data<-0,nrow<-k,ncol<-k)
	
	## Fill matrix with survival & fecundity vectors
	diag(ad[2:age.max,1:(age.max-1)]) <- s1p.d ## adds s1p.d to female quadrant
	diag(ad[(age.max+2):k,(age.max+1):(k-1)]) <- s1p.d ## adds s1p.d to male quadrant
	ad[1,2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds female fecundity
	ad[(age.max+1),2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds male fecundity

	## max r
	rd <- max.r(ad)

	ssdd <- stable.stage.dist(ad)

	## pick initial vectors
	nd<-ssdd*Nd

	##set population size year step vector
	popd.vec <- rep(1,tlimit+1)
	popd.vec[1] <- sum(nd)

	##set female population size year step vector
	femd.vec<-rep(1,tlimit+1)
	femd.vec[1]<-sum(nd[1:age.max])

	femda.vec<-rep(1,tlimit+1)
	femda.vec[1]<-sum(nd[3:age.max])

	##set male population size year step vector
	mald.vec<-rep(1,tlimit+1)
	mald.vec[1]<-sum(nd[(age.max+1):k])

	malda.vec<-rep(1,tlimit+1)
	malda.vec[1]<-sum(nd[21:k])

	## set large male year step vector
	male6p.vec <- rep(0,tlimit+1)
	male6p.vec[1]<-sum(nd[23:k])

	##set year step vector
	yrd.vec <- rep(1,tlimit+1)
	yrd.vec[1] <- 0

	for (j in 1:tlimit) {
    		yrd.vec[j+1] <- j
	}

	
	##then iterate
		for (ii in 1:tlimit) {
			nd <- ad %*% nd

		## Calculate severity of a potential catastrophe
		sev.cat.input <- runif(1,0,(max(sev.cat.prob)))
		sev.cat <- ((sev.cat.coeff[1]*(sev.cat.input^sev.cat.coeff[2])))/100
		##sev.cat
		if (sev.cat > 1) sev.cat <- 0.99
		
		## calculate whether a catastrophe occurs
		cat.occur <- runif(1,0,1)
		##cat.occur
		if (cat.occur <= p.cat) nd <- nd*(1-sev.cat) else nd <- nd
		
		## Set any negative or Na values in nd to 0
		zero.sub <- which(nd<0)
		nd[zero.sub] <- 0

		na.sub <- which(is.na(nd))
		nd[na.sub] <- 0

		## Stochastic measure (variation in rainfall)
		rain.mult <- runif(1,min=rain.mult.min,max=rain.mult.max)

	        ## Set negative density feedback function for the matrix
	
	        ## Redefine survival probabilities
	        s0.d <- rain.mult*as.numeric((SSfpl(sum(nd),s0.param[1],s0.param[2],s0.param[3],s0.param[4])))
	            if (s0.d < s0.max) s0.d <- s0.max else s0.d <- s0.d
		    if (s0.d > s0.min) s0.d <- s0.min else s0.d <- s0.d
	            
	        s1p.d <- rain.mult*as.numeric((SSfpl(sum(nd),s1p.param[1],s1p.param[2],s1p.param[3],s1p.param[4])))
	            if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
		    if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d

		## re-define fecundity a coefficient
	        a.md <- as.numeric(SSfpl(sum(nd),af.param[1],af.param[2],af.param[3],af.param[4]))
	            if (a.md < af.max) a.md <- af.max else a.md <- a.md
	            if (a.md > af.min) a.md <- af.min else a.md <- a.md

	        ## Re-define fertility vector
	        md.vec <- rep(0,age.max)
		for (j in (1:age.max)) {
	            md <- a.md*exp(-0.5*(log(j/x0.fert)/b.fert)^2)
	            md.vec[j] <- md
	        }

		## Adjust for demographic stochasticity (binomial & poisson resampling)
		n.mf <- md.vec*nd[1:age.max]
		sum.n.mf <- sum(n.mf)
		if ((sum(n.mf)) == 0) sum.n.mf <- 1
		md.vec <- sum(md.vec)*(rpois(age.max,n.mf)/sum.n.mf)

	        n.1 <- round(sum(nd))
		if (n.1 == 0) n.1 <- 1
		s1p.d <- rbinom(1,n.1,s1p.d)/n.1		

		n.0 <- round(sum(md.vec*n[1:age.max]))
		s0.d <- rbinom(1,n.0,s0.d)/n.0
		if (n.0 == 0) n.0 <-1

		## Tertiary sex ratio
	        if (sum(nd[1:age.max] / n.1) > tsr.max) md.vec <- rep(0,age.max) else md.vec <- md.vec

		## Fill matrix with survival & fecundity vectors
		diag(ad[2:age.max,1:(age.max-1)]) <- s1p.d ## adds s1p.d to female quadrant
		diag(ad[(age.max+2):k,(age.max+1):(k-1)]) <- s1p.d ## adds s1p.d to male quadrant
		ad[1,2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds female fecundity
		ad[(age.max+1),2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds male fecundity
	        
		     popd.vec[ii+1]<-(sum(nd))    
		     femd.vec[ii+1]<-sum(nd[1:age.max]) ## all females
		     femda.vec[ii+1]<-sum(nd[3:age.max]) ## adult (>2 years) females
		     mald.vec[ii+1]<-sum(nd[(age.max+1):k]) ## all males
		     malda.vec[ii+1]<-sum(nd[21:k]) ## adult (>3 years) males
		     male6p.vec[ii+1]<-sum(nd[23:k]) ## males >5 years
	     
		} ## end ii loop

	log.popd.vec<-log10(popd.vec)

	## Place values in storage matrices
	popd.sperm[p,] <- popd.vec
	femd.sperm[p,] <- femd.vec
	mald.sperm[p,] <- mald.vec
	male6p.sperm[p,] <- male6p.vec

	print(p)

} ## end p loop (stochastic)

## Mean & confidence intervals over p permutations
popd.mean <- rep(0,tlimit+1); popd.up <- rep(0,tlimit+1); popd.lo <- rep(0,tlimit+1)
femd.mean <- rep(0,tlimit+1); femd.up <- rep(0,tlimit+1); femd.lo <- rep(0,tlimit+1)
mald.mean <- rep(0,tlimit+1); mald.up <- rep(0,tlimit+1); mald.lo <- rep(0,tlimit+1)
male6p.mean <- rep(0,tlimit+1); male6p.up <- rep(0,tlimit+1); male6p.lo <- rep(0,tlimit+1)


for (d in 1:(tlimit+1)) {
	popd.mean[d] <- mean(popd.sperm[,d])
	popd.up[d] <- quantile(popd.sperm[,d],probs=0.975)
	popd.lo[d] <- quantile(popd.sperm[,d],probs=0.025)

	femd.mean[d] <- mean(femd.sperm[,d])
	femd.up[d] <- quantile(femd.sperm[,d],probs=0.975)
	femd.lo[d] <- quantile(femd.sperm[,d],probs=0.025)

	mald.mean[d] <- mean(mald.sperm[,d])
	mald.up[d] <- quantile(mald.sperm[,d],probs=0.975)
	mald.lo[d] <- quantile(mald.sperm[,d],probs=0.025)

	male6p.mean[d] <- mean(male6p.sperm[,d])
	male6p.up[d] <- quantile(male6p.sperm[,d],probs=0.975)
	male6p.lo[d] <- quantile(male6p.sperm[,d],probs=0.025)

} ## end d loop

if (sum(popd.up) == 0) popd.up <- popd.mean else popd.up <- popd.up
if (sum(popd.lo) == 0) popd.lo <- popd.mean else popd.lo <- popd.lo

if (sum(femd.up) == 0) femd.up <- femd.mean else femd.up <- femd.up
if (sum(femd.lo) == 0) femd.lo <- femd.mean else femd.lo <- femd.lo

if (sum(mald.up) == 0) mald.up <- mald.mean else mald.up <- mald.up
if (sum(mald.lo) == 0) mald.lo <- mald.mean else mald.lo <- mald.lo

if (sum(male6p.up) == 0) male6p.up <- male6p.mean else male6p.up <- male6p.up
if (sum(male6p.lo) == 0) male6p.lo <- male6p.mean else male6p.lo <- male6p.lo

## Make stochastic density-dependent plots
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(yr.vec[0:tlimit],popd.mean[1:tlimit],xlab="year",ylab="N",xlim=range(-0.5:(max(tlimit))),ylim=c(0,(max(popd.up)+(0.25*(max(popd.up))))),type='l')
title(main="Total Population")
lines(yr.vec[0:tlimit],popd.up[1:tlimit],col="red")
lines(yr.vec[0:tlimit],popd.lo[1:tlimit],col="red")
plot(yr.vec[0:tlimit],femd.mean[1:tlimit],xlab="year",ylab="N",xlim=range(-0.5:(max(tlimit))),ylim=c(0,(max(popd.mean)+(0.25*(max(femd.mean))))),type='l')
lines(yr.vec[0:tlimit],femd.up[1:tlimit],col="red")
lines(yr.vec[0:tlimit],femd.lo[1:tlimit],col="red")
title(main="Total Adult Females")
plot(yr.vec[0:tlimit],mald.mean[1:tlimit],xlab="year",ylab="N",xlim=range(-0.5:(max(tlimit))),ylim=c(0,(max(popd.mean)+(0.25*(max(mald.mean))))),type='l')
lines(yr.vec[0:tlimit],mald.up[1:tlimit],col="red")
lines(yr.vec[0:tlimit],mald.lo[1:tlimit],col="red")
title(main="Total Adult Males")
plot(yr.vec[0:tlimit],male6p.mean[1:tlimit],xlab="year",ylab="N",xlim=range(-0.5:(max(tlimit))),ylim=c(0,(max(popd.mean)+(0.25*(max(male6p.mean))))),type='l')
lines(yr.vec[0:tlimit],male6p.up[1:tlimit],col="red")
lines(yr.vec[0:tlimit],male6p.lo[1:tlimit],col="red")
title(main="Total Mature Males")
par(mfrow=c(1,2))





## Parameter value checks

## logistic.dens <- dlogis(quant, location = 0.01, scale = 10, log = FALSE) ## estimates logistic density function with long left tail
logistic.dens <- 0 + ((0.997)/(1+(quant/1.8217)^3.9151))
log.dens <- (logistic.dens - min(logistic.dens)) / max(logistic.dens - min(logistic.dens))

## Set min and max s0 values
s0.min <- 0.93
s0.max <- 0.62

## Calculate logistic function between min and max s0
s0.vec <- (seq(s0.min,s0.max,-((s0.min-s0.max)/pop.step)))
range.s0.vec <- range(s0.vec)[2] - range(s0.vec)[1]
s0.vec.logis <- (range.s0.vec * log.dens) + min(s0.vec)

## Set min and max s1+ values
s1p.min <- 0.995
s1p.max <- 0.92

## Calculate logistic function between min and max s1p
s1p.vec <- (seq(s1p.min,s1p.max,-((s1p.min-s1p.max)/pop.step)))
range.s1p.vec <- range(s1p.vec)[2] - range(s1p.vec)[1]
s1p.vec.logis <- (range.s1p.vec * log.dens) + min(s1p.vec)

## Set min and max fecundity a coefficient values
af.min <- 0.65
af.max <- coeff.fit.m.vec[1]

## Calculate logistic function between min and max
af.vec <- (seq(af.min,af.max,-((af.min-af.max)/pop.step)))
range.af.vec <- range(af.vec)[2] - range(af.vec)[1]
af.vec.logis <- (range.af.vec * log.dens) + min(af.vec)

## Plot dd functions
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(num.vec,s0.vec.logis,xlab="N",ylab="s0",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
title(main="Logistic DD Functions")
plot(num.vec,s1p.vec.logis,xlab="N",ylab="s1+",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
plot(num.vec,af.vec.logis,xlab="N",ylab="a.f",xlim=range(-0.5:((max(num.vec))+1)),ylim=range(c(0,1)),type='l')
par(mfrow=c(1,2))

log.func <- data.frame(num.vec,s0.vec.logis,s1p.vec.logis,af.vec.logis) ## put new vectors into data.frame

## Estimate 4-parameter logistic model coefficients
s0.param <- getInitial(s0.vec.logis ~ SSfpl(num.vec, a.s0, b.s0, xmid.s0, scal.s0), data=log.func)
s1p.param <- getInitial(s1p.vec.logis ~ SSfpl(num.vec, a.s1p, b.s1p, xmid.s1p, scal.s1p), data=log.func)
af.param <- getInitial(af.vec.logis ~ SSfpl(num.vec, a.af, b.af, xmid.af, scal.af), data=log.func)


## Input population size to test
Ntest <- 9000

## define survival
s0.d <- SSfpl(sum(Ntest),s0.param[1],s0.param[2],s0.param[3],s0.param[4])
	if (s0.d < s0.max) s0.d <- s0.max else s0.d <- s0.d
	if (s0.d > s0.min) s0.d <- s0.min else s0.d <- s0.d

s1p.d <- SSfpl(sum(Ntest),s1p.param[1],s1p.param[2],s1p.param[3],s1p.param[4])
	if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
	if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d
        
## define fecundity
a.md <- SSfpl(sum(Ntest),af.param[1],af.param[2],af.param[3],af.param[4])
	if (a.md < af.max) a.md <- af.max else a.md <- a.md
	if (a.md > af.min) a.md <- af.min else a.md <- a.md

## Re-define fertility vector
md.vec <- rep(0,age.max); b.fert <- coeff.fit.m.vec[3]; x0.fert <- coeff.fit.m.vec[2]
	for (j in (1:age.max)) {
		md<-a.md*exp(-0.5*(log(j/x0.fert)/b.fert)^2)
		md.vec[j] <- md
	}

## Plot test values
row <- 1
col <- 1
par(mfrow=c(row,col))
plot(age.vec[2:18],md.vec,xlab="N",ylab="m",xlim=range(0:17),ylim=range(c(0,1)),type='l')
abline(h=0.65,col="red")
abline(h=(max(md.vec)))
par(mfrow=c(1,2))



#####################################
## Minimum viable population with DD
#####################################

## ***********************************************************************************************
##K population size
Nd <- 9000

## Re-define projection interval
tlimit <- 100

## Set number of iterations to estimate mean & confidence interval for projections
perm <- 1000 ## number of permutations

## Set quasi-extinction threshold (proportion of K)
Q.thresh <- 20
ex.1 <- 1 ## How many years after introduction to calculate Q?

## Define limits for start vector
init.min <- 20
init.max <- 200

## Tertiary sex ratio (0.83 = balanced in polygynous ungulates)
tsr.max <- 0.95 ## Tertiary sex ratio (breeding females:breeding adults)

## *********************************************************************************
## ***************
## End user input
## ***************
## *********************************************************************************

init.array <- seq(init.min,init.max,10) ## initial start population array array
linit.array <- length(init.array)

	min.popd.vec <- rep(0,length(init.array)) ## Minimum population size vector
	lambdad.vec <- rep(0,length(init.array)) ## mean population rate of change (r)

	## Create s vectors
	p.ext.vec <- rep(0,linit.array)
	popd.min.mean.vec <- rep(0,linit.array)
	popd.min.up.vec <- rep(0,linit.array)
	popd.min.lo.vec <- rep(0,linit.array)
	popd.end.mean.vec <- rep(0,linit.array)
	popd.end.up.vec <- rep(0,linit.array)
	popd.end.lo.vec <- rep(0,linit.array)
	init.mean.vec <- rep(0,linit.array)

	## Start simulation over init.array

	for (s in 1:linit.array) {

	N.st <- init.array[s] ## start process with new initial population size


	## Start simulation over perm permutations

	## Create p vectors
	popd.min.vec <- rep(0,perm)
	popd.end.vec <- rep(0,perm)
	qextinct.vec <- rep(0,perm)
	init.mean <- rep(0,perm)


	for (p in 1:perm) {

	## Stochastic variation in survival due to rainfall (for permutations)
	rain.multp <- runif(1,min=rain.mult.min,max=rain.mult.max)

	## pick initial vectors
	##nd<-ssdd*N.st
	af.cl <- age.max - 2 ## adult female classes
	am.cl <- age.max - 3 ## adult male classes
	a.cl <- af.cl + am.cl
	a.cl.vec <- seq(1,a.cl)
	nd <- rep(0,k)

	a.ord <- (sample(a.cl.vec))
	for (ai in 1:N.st) {
		ord.sub <- which(a.ord == ai)
		if ((sum(nd)) > N.st) next
		val <- floor((runif(1,0,0.2))*N.st)
		if (ord.sub < 3) next
		if (ord.sub == 17) next
		if (ord.sub == 18) next
		if (ord.sub == 19) next
		nd[ord.sub] <- val
	}
	nd <- ifelse((is.na(nd)) == TRUE,0,nd)
##	nd
##	sum(nd)
	nd.st <- nd

	## define start survival

	        s0.d <- rain.multp*as.numeric((SSfpl(sum(N.st),s0.param[1],s0.param[2],s0.param[3],s0.param[4])))
	            if (s0.d < s0.max) s0.d <- s0.max else s0.d <- s0.d
		    if (s0.d > s0.min) s0.d <- s0.min else s0.d <- s0.d

	        s1p.d <- rain.multp*as.numeric((SSfpl(sum(N.st),s1p.param[1],s1p.param[2],s1p.param[3],s1p.param[4])))
	            if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
		    if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d
        
	## define start fecundity
	        a.md <- rain.multp*as.numeric((SSfpl(sum(N.st),af.param[1],af.param[2],af.param[3],af.param[4])))
	            if (a.md < af.max) a.md <- af.max else a.md <- a.md
	            if (a.md > af.min) a.md <- af.min else a.md <- a.md
	            
	        ## Re-define fertility vector
	        md.vec <- rep(0,age.max); b.fert <- coeff.fit.m.vec[3]; x0.fert <- coeff.fit.m.vec[2]
	        for (j in (1:age.max)) {
	            md<-a.md*exp(-0.5*(log(j/x0.fert)/b.fert)^2)
	            md.vec[j] <- md
	        }

	## Adjust for demographic stochasticity (binomial & poisson resampling)
	n.mf <- 2*md.vec*nd.st[1:age.max]
	n.f <- sum(nd.st[1:age.max])
	if (n.f == 0) n.f <- 1
	sum.n.mf <- sum(n.mf)
	if (sum.n.mf == 0) sum.n.mf <- 1
	if ((is.na(sum.n.mf)) == TRUE) sum.n.mf <- 1

	md.vec <- sum(md.vec)*(rpois(age.max,n.f*n.mf)/(n.f*sum.n.mf))

		for (ms in 1:(length(md.vec))) {
			if (md.vec[ms] > af.min) md.vec[ms] <- af.min else md.vec[ms] <- md.vec[ms]
		}
	
	n.1 <- round(sum(nd.st))
	if (n.1 == 0) n.1 <- 1
	if ((is.na(n.1)) == TRUE) n.1 <- 1
	s1p.d <- rbinom(1,n.1,s1p.d)/n.1		
##        if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
	if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d

	n.0 <- round(sum(md.vec*n[1:age.max]))
	if (n.0 == 0) n.0 <- 1
	if ((is.na(n.0)) == TRUE) n.0 <- 1
	s0.d <- rbinom(1,n.0,s0.d)/n.0

	## The normal matrix
	ad <- matrix(data<-0,nrow<-k,ncol<-k)
	
	## Fill matrix with survival & fecundity vectors
	diag(ad[2:age.max,1:(age.max-1)]) <- s1p.d ## adds s1p.d to female quadrant
	diag(ad[(age.max+2):k,(age.max+1):(k-1)]) <- s1p.d ## adds s1p.d to male quadrant
	ad[1,2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds female fecundity
	ad[(age.max+1),2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds male fecundity

	## max r
	rd <- max.r(ad)

	##set population size year step vector
	popd.vec <- rep(1,tlimit+1)
	popd.vec[1] <- sum(nd)

	##set year step vector
	yrd.vec <- rep(1,tlimit+1)
	yrd.vec[1] <- 0

	for (j in 1:tlimit) {
    		yrd.vec[j+1] <- j
	}

	##then iterate
		for (ii in 1:tlimit) {
			nd <- ad %*% nd

		## Calculate severity of a potential catastrophe
		sev.cat.input <- runif(1,0,(max(sev.cat.prob)))
		sev.cat <- ((sev.cat.coeff[1]*(sev.cat.input^sev.cat.coeff[2])))/100
		##sev.cat
		if (sev.cat > 1) sev.cat <- 0.99

		## calculate whether a catastrophe occurs
		if (runif(1,0,1) <= p.cat) nd <- nd*(1-sev.cat) else nd <- nd

		## Set any negative values in nd to 0
		zero.sub <- which(nd<0)
		nd[zero.sub] <- 0

		na.sub <- which(is.na(nd))
		nd[na.sub] <- 0

		if (sum(nd) == 0) next

		## Stochastic measure (variation in rainfall)
		rain.mult <- runif(1,min=rain.mult.min,max=rain.mult.max)
				
	        ## Set negative density feedback function for the matrix
	
	        ## Redefine survival probabilities
	        s0.d <- rain.mult*(SSfpl(sum(nd),s0.param[1],s0.param[2],s0.param[3],s0.param[4]))
	            if (s0.d < s0.max) s0.d <- s0.max else s0.d <- s0.d
		    if (s0.d > s0.min) s0.d <- s0.min else s0.d <- s0.d

	        s1p.d <- rain.mult*(SSfpl(sum(nd),s1p.param[1],s1p.param[2],s1p.param[3],s1p.param[4]))
	            if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
		    if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d

		## re-define fecundity a coefficient
	        a.md <- SSfpl(sum(nd),af.param[1],af.param[2],af.param[3],af.param[4])
	            if (a.md < af.max) a.md <- af.max else a.md <- a.md
	            if (a.md > af.min) a.md <- af.min else a.md <- a.md

	        ## Re-define fertility vector
	        md.vec <- rep(0,age.max)
		for (j in (1:age.max)) {
	            md <- a.md*exp(-0.5*(log(j/x0.fert)/b.fert)^2)
	            md.vec[j] <- md
	        }

	        ## Tertiary sex ratio
	        	if (sum(nd[1:age.max] / (1+sum(nd))) > tsr.max) md.vec <- rep(0,age.max) else md.vec <- md.vec

		## Adjust for demographic stochasticity (binomial & poisson resampling)
		n.mf <- 2*md.vec*nd.st[1:age.max]
		n.f <- sum(nd.st[1:age.max])
		if (n.f == 0) n.f <- 1
		sum.n.mf <- sum(n.mf)
		if (sum.n.mf == 0) sum.n.mf <- 1
		if ((is.na(sum.n.mf)) == TRUE) sum.n.mf <- 1

		md.vec <- sum(md.vec)*(rpois(age.max,n.f*n.mf)/(n.f*sum.n.mf))

			for (ms in 1:(length(md.vec))) {
				if (md.vec[ms] > af.min) md.vec[ms] <- af.min else md.vec[ms] <- md.vec[ms]
			}
		
		n.1 <- round(sum(nd))
		if (n.1 == 0) n.1 <- 1
		if ((is.na(n.1)) == TRUE) n.1 <- 1
		s1p.d <- rbinom(1,n.1,s1p.d)/n.1		
	        if (s1p.d < s1p.max) s1p.d <- s1p.max else s1p.d <- s1p.d
		if (s1p.d > s1p.min) s1p.d <- s1p.min else s1p.d <- s1p.d

		n.0 <- round(sum(md.vec*n[1:age.max]))
		if (n.0 == 0) n.0 <- 1
		if ((is.na(n.0)) == TRUE) n.0 <- 1
		s0.d <- rbinom(1,n.0,s0.d)/n.0

		## Fill matrix with survival & fecundity vectors
		diag(ad[2:age.max,1:(age.max-1)]) <- s1p.d ## adds s1p.d to female quadrant
		diag(ad[(age.max+2):k,(age.max+1):(k-1)]) <- s1p.d ## adds s1p.d to male quadrant
		ad[1,2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds female fecundity
		ad[(age.max+1),2:(age.max)] <- s0.d * md.vec[1:(age.max-1)] ## adds male fecundity
	        
	        popd.vec[ii+1] <- (sum(nd))    

		} ## end ii loop

	# Plot ii loop results
	row <- 1
	col <- 1
	par(mfrow=c(row,col))
	plot(yrd.vec,popd.vec,xlab="Years",ylab="Population size",type="l")
	title(main="Population projection")
	par(mfrow=c(1,1))

	## Did the population drop below the quasi-extinction threshold (Q.thresh)?
	qext.sub <- as.numeric(which((popd.vec[ex.1:(length(popd.vec))]) < Q.thresh))
	sum.qext.sub <- sum(qext.sub)
	if (sum.qext.sub > 0) qextinct <- 1 else qextinct <- 0

	qextinct.vec[p] <- qextinct

	## Place minimum population size achieved in p vector
	popd.min.vec[p] <- min(ifelse(popd.vec>0,popd.vec,0))

	## Put end population sizes into p vector
	popd.end.vec[p] <- popd.vec[ii]

	## Put mean initial size into p vector
	init.mean[p] <- sum(nd.st)

	##print(p)

	} ## end p loop

	## Mean inital vector
	init.mean.vec[s] <- mean(init.mean)

	## Probability of extinction over p permutations
	p.ext <- sum(qextinct.vec)/perm

	## Estimate mean & confidence intervals for minimum population size
	popd.min.mean <- mean(popd.min.vec)
	popd.min.up <- as.numeric(quantile(popd.min.vec,probs=0.975))
	popd.min.lo <- as.numeric(quantile(popd.min.vec,probs=0.025))

	## Estimate mean & confidence intervals for end population size
	popd.end.mean <- mean(popd.end.vec)
	popd.end.up <- as.numeric(quantile(popd.end.vec,probs=0.975))
	popd.end.lo <- as.numeric(quantile(popd.end.vec,probs=0.025))

	print("initial population size")	
	print(init.array[s])

	## Place probability of extinction into s vector
	p.ext.vec[s] <- p.ext

	## Place minimum population sizes & confidence intervals into s vector
	popd.min.mean.vec[s] <- popd.min.mean
	popd.min.up.vec[s] <- popd.min.up
	popd.min.lo.vec[s] <- popd.min.lo

	## Place end population sizes & confidence intervals into s vector
	popd.end.mean.vec[s] <- popd.end.mean
	popd.end.up.vec[s] <- popd.end.up
	popd.end.lo.vec[s] <- popd.end.lo

	} ## end s loop


## Plot results
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(init.mean.vec,popd.end.mean.vec,xlab="Initial population size",ylab="End population size",ylim=range(c((min(popd.end.lo.vec)),(max(popd.end.up.vec)))),type="l")
lines(init.mean.vec,popd.end.up.vec,col="red")
lines(init.mean.vec,popd.end.lo.vec,col="red")
title(main="End population size")
plot(init.mean.vec,popd.min.mean.vec,xlab="Initial population size",ylab="Minimum population size",ylim=range(c((min(popd.min.lo.vec)),(max(popd.min.up.vec)))),type="l")
lines(init.mean.vec,popd.min.up.vec,col="red")
lines(init.mean.vec,popd.min.lo.vec,col="red")
title(main="Minimum population size")
plot(init.mean.vec,p.ext.vec,xlab="Initial population size",ylab="Probability",ylim=range(c(0,0.5)),type="l")
title(main="Probability of Quasi-extinction")
par(mfrow=c(1,1))

p.ext.vec1 <- p.ext.vec

row <- 1
col <- 2
par(mfrow=c(row,col),pty="s")
plot(init.mean.vec,p.ext.vec1,xlab="Initial population size",ylab="Probability",ylim=range(c(0,0.5)),type="l")
boxed.labels(30,0.48,"a",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)
plot(init.array,p.ext.vec,xlab="Initial population size",ylab="Probability",ylim=range(c(0,1.0)),type="l")
boxed.labels(25,0.95,"b",border=FALSE,xpad=0,ypad=0,font=2,cex=1.5)
par(mfrow=c(1,1))



####################################################
## Minimum viable population - no DD (stable matrix)
####################################################

## ***********************************************************************************************
##K population size
Nd <- 9000

## Re-define projection interval
tlimit <- 100

## Set number of iterations to estimate mean & confidence interval for projections
perm <- 1000 ## number of permutations

## Set quasi-extinction threshold (proportion of K)
Q.thresh <- 2
ex.1 <- 1 ## How many years after introduction to calculate Q?

## Define limits for start vector
init.min <- 10
init.max <- 200

## Tertiary sex ratio (0.83 = balanced in polygynous ungulates)
tsr.max <- 0.95 ## Tertiary sex ratio (breeding females:breeding adults)

## *********************************************************************************
## ***************
## End user input
## ***************
## *********************************************************************************

init.array <- seq(init.min,init.max,5) ## initial start population array array
linit.array <- length(init.array)

	min.popd.vec <- rep(0,length(init.array)) ## Minimum population size vector
	lambdad.vec <- rep(0,length(init.array)) ## mean population rate of change (r)

	## Create s vectors
	p.ext.vec <- rep(0,linit.array)
	popd.min.mean.vec <- rep(0,linit.array)
	popd.min.up.vec <- rep(0,linit.array)
	popd.min.lo.vec <- rep(0,linit.array)
	popd.end.mean.vec <- rep(0,linit.array)
	popd.end.up.vec <- rep(0,linit.array)
	popd.end.lo.vec <- rep(0,linit.array)

	## Start simulation over init.array

	for (s in 1:linit.array) {

	N.st <- init.array[s] ## start process with new initial population size


	## Start simulation over perm permutations

	## Create p vectors
	popd.min.vec <- rep(0,perm)
	popd.end.vec <- rep(0,perm)
	qextinct.vec <- rep(0,perm)

	for (p in 1:perm) {

	## Stochastic measure (variation in rainfall)
	rain.multp <- runif(1,min=rain.mult.min,max=rain.mult.max)

	## pick initial vectors
	nd<-ssdd*N.st

	## define start survival

	## stable survival vector after rainfall correction
	s.stable <- s.vec
	s.stable[1] <- 0.62
	s.stable <- rain.multp*s.stable
	ls.stable <- length(s.stable)

		for (q in 1:(ls.stable-1)) {
			if (s.stable[q+1] > s1p.min) s.stable[q+1] <- s1p.min else s.stable[q+1] <- s.stable[q+1]
		}

        
	## define start fecundity
	m.stable <- mf.new.vec	

	## Adjust for demographic stochasticity (binomial & poisson resampling)
	n.mf <- 2*m.stable[2:(age.max+1)]*nd[1:age.max]
	n.f <- sum(nd[1:age.max])
	sum.n.mf <- sum(n.mf)
	if (sum.n.mf == 0) sum.n.mf <- 1
	if ((is.na(sum.n.mf)) == TRUE) sum.n.mf <- 1
	mds.vec <- sum(m.stable)*(rpois(age.max,(n.mf*n.f))/(n.f*sum.n.mf))

		for (ms in 1:(length(mds.vec))) {
			if (mds.vec[ms] > af.min) mds.vec[ms] <- af.min else mds.vec[ms] <- mds.vec[ms]
		}

	n.s <- round(sum(nd))
	if (n.s == 0) n.s <- 1
	if ((is.na(n.s)) == TRUE) n.s <- 1
	s.rbinom <- rep(0,ls.stable)
		for (ss in 1:(ls.stable)) {
			s.rbinom[ss] <- rbinom(1,n.s,s.stable[ss])
		}
	sds.vec <- s.rbinom/(n.s)

		for (sf in 1:(ls.stable-1)) {
			if (sds.vec[sf+1] > s1p.min) sds.vec[sf+1] <- s1p.min else sds.vec[sf+1] <- sds.vec[sf+1]
		}
##mds.vec
##sds.vec
##plot(mds.vec,type="l",ylim=range(c(0,1)))


	## The start matrix
	as <- matrix(data<-0,nrow<-k,ncol<-k)
	
	## Add survival vectors to matrix
	diag(as[2:age.max,1:(age.max-1)]) <- sds.vec[2:age.max] ## adds s.vec to female quadrant
	diag(as[(age.max+2):k,(age.max+1):(k-1)]) <- sds.vec[2:age.max] ## adds s.vec to male quadrant
	as[1,1:(age.max)] <- sds.vec[1] * mds.vec[1:(age.max)] ## adds female fecundity
	as[(age.max+1),1:(age.max)] <- sds.vec[1] * mds.vec[1:(age.max)] ## adds male fecundity

	## max r
	rd <- max.r(ad)

	##set population size year step vector
	popd.vec <- rep(0,tlimit+1)
	popd.vec[1] <- sum(nd)

	##set year step vector
	yrd.vec <- rep(1,tlimit+1)
	yrd.vec[1] <- 0

	for (j in 1:tlimit) {
    		yrd.vec[j+1] <- j
	}

	##then iterate
		for (ii in 1:tlimit) {
			nd <- as %*% nd

		## Calculate severity of a potential catastrophe
		sev.cat.input <- runif(1,0,(max(sev.cat.prob)))
		sev.cat <- ((sev.cat.coeff[1]*(sev.cat.input^sev.cat.coeff[2])))/100
		##sev.cat
		if (sev.cat > 1) sev.cat <- 0.99

		## calculate whether a catastrophe occurs
		if (runif(1,0,1) <= p.cat) nd <- nd*(1-sev.cat) else nd <- nd

		## Set any negative values in nd to 0
		zero.sub <- which(nd<0)
		nd[zero.sub] <- 0

		na.sub <- which(is.na(nd))
		nd[na.sub] <- 0

		if (sum(nd) == 0) next

		## Stochastic measure (variation in rainfall)
		rain.mult <- runif(1,min=rain.mult.min,max=rain.mult.max)
		
		## stable survival vector after rainfall correction
		s.stable <- s.vec
		s.stable[1] <- 0.62
		s.stable <- rain.mult*s.stable

			for (q in 1:(ls.stable-1)) {
				if (s.stable[q+1] > s1p.min) s.stable[q+1] <- s1p.min else s.stable[q+1] <- s.stable[q+1]
			}
				
		## Adjust for demographic stochasticity (binomial & poisson resampling)
		n.mf <- 2*m.stable[2:(age.max+1)]*nd[1:age.max]
		n.f <- sum(nd[1:age.max])
		sum.n.mf <- sum(n.mf)
		if (sum.n.mf == 0) sum.n.mf <- 1
		if ((is.na(sum.n.mf)) == TRUE) sum.n.mf <- 1
		mds.vec <- sum(m.stable)*(rpois(age.max,(n.mf*n.f))/(sum.n.mf*n.f))

			for (mk in 1:(length(mds.vec))) {
				if (mds.vec[mk] > af.min) mds.vec[mk] <- af.min else mds.vec[mk] <- mds.vec[mk]
			}

		n.s <- round(sum(nd))
		if (n.s == 0) n.s <- 1
		if ((is.na(n.s)) == TRUE) n.s <- 1	
		ls.stable <- length(s.stable)
		s.rbinom <- rep(0,ls.stable)
			for (ss in 1:(ls.stable)) {
				s.rbinom[ss] <- rbinom(1,n.s,s.stable[ss])
			}
		sds.vec <- s.rbinom/(n.s)

			for (sf in 1:(ls.stable-1)) {
				if (sds.vec[sf+1] > s1p.min) sds.vec[sf+1] <- s1p.min else sds.vec[sf+1] <- sds.vec[sf+1]
			}
##mds.vec
##sds.vec
##plot(mds.vec,type="l",ylim=range(c(0,1)))

	        ## Tertiary sex ratio
	        	if (sum(nd[1:age.max] / (1+sum(nd))) > tsr.max) m.stable <- 0 else m.stable <- m.stable

		## The adjusted matrix
		as <- matrix(data<-0,nrow<-k,ncol<-k)
	
		## Add survival vectors to matrix
		diag(as[2:age.max,1:(age.max-1)]) <- sds.vec[2:age.max] ## adds s.vec to female quadrant
		diag(as[(age.max+2):k,(age.max+1):(k-1)]) <- sds.vec[2:age.max] ## adds s.vec to male quadrant
		as[1,1:(age.max)] <- sds.vec[1] * mds.vec[1:(age.max)] ## adds female fecundity
		as[(age.max+1),1:(age.max)] <- sds.vec[1] * mds.vec[1:(age.max)] ## adds male fecundity

	        popd.vec[ii+1] <- (sum(nd))    

		} ## end ii loop

	# Plot ii loop results
	row <- 1
	col <- 1
	par(mfrow=c(row,col))
	plot(yrd.vec,popd.vec,xlab="Years",ylab="Population size",type="l")
	title(main="Population projection")
	par(mfrow=c(1,1))

	## Did the population drop below the quasi-extinction threshold (Q.thresh)?
	qext.sub <- as.numeric(which((popd.vec[ex.1:(length(popd.vec))]) < Q.thresh))
	sum.qext.sub <- sum(qext.sub)
	if (sum.qext.sub > 0) qextinct <- 1 else qextinct <- 0

	qextinct.vec[p] <- qextinct

	## Place minimum population size achieved in p vector
	popd.min.vec[p] <- min(ifelse(popd.vec>0,popd.vec,0))

	## Put end population sizes into p vector
	popd.end.vec[p] <- popd.vec[ii]

	##print(p)

	} ## end p loop

	## Probability of extinction over p permutations
	p.ext <- sum(qextinct.vec)/perm

	## Estimate mean & confidence intervals for minimum population size
	popd.min.mean <- mean(popd.min.vec)
	popd.min.up <- as.numeric(quantile(popd.min.vec,probs=0.975))
	popd.min.lo <- as.numeric(quantile(popd.min.vec,probs=0.025))

	## Estimate mean & confidence intervals for end population size
	popd.end.mean <- mean(popd.end.vec)
	popd.end.up <- as.numeric(quantile(popd.end.vec,probs=0.975))
	popd.end.lo <- as.numeric(quantile(popd.end.vec,probs=0.025))

	print("initial population size")	
	print(init.array[s])

	## Place probability of extinction into s vector
	p.ext.vec[s] <- p.ext

	## Place minimum population sizes & confidence intervals into s vector
	popd.min.mean.vec[s] <- popd.min.mean
	popd.min.up.vec[s] <- popd.min.up
	popd.min.lo.vec[s] <- popd.min.lo

	## Place end population sizes & confidence intervals into s vector
	popd.end.mean.vec[s] <- popd.end.mean
	popd.end.up.vec[s] <- popd.end.up
	popd.end.lo.vec[s] <- popd.end.lo

	} ## end s loop


## Plot results
row <- 2
col <- 2
par(mfrow=c(row,col))
plot(init.array,popd.end.mean.vec,xlab="Initial population size",ylab="End population size",ylim=range(c((min(popd.end.lo.vec)),(max(popd.end.up.vec)))),type="l")
lines(init.array,popd.end.up.vec,col="red")
lines(init.array,popd.end.lo.vec,col="red")
title(main="End population size")
plot(init.array,popd.min.mean.vec,xlab="Initial population size",ylab="Minimum population size",ylim=range(c((min(popd.min.lo.vec)),(max(popd.min.up.vec)))),type="l")
lines(init.array,popd.min.up.vec,col="red")
lines(init.array,popd.min.lo.vec,col="red")
title(main="Minimum population size")
plot(init.array,p.ext.vec,xlab="Initial population size",ylab="Probability",ylim=range(c(0,1)),type="l")
title(main="Probability of Quasi-extinction")
par(mfrow=c(1,2))
