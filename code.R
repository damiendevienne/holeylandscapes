### April 2024
### D de Vienne/J Joseph
### Code for answer to Dochtermann et al. (2023) on holey landscapes
### 
### This code contains : 
### - source of external functions
### - Empirica data and L2/L1 computations
### - Gaussian simulations and L2/L1 computations
### - Holey simulations and L2/L1 computations
### - Plotting function to create Fig. 1 of our response.
### 


### Original functions from Dochtermann et al.
source("Wright_vs_Holey(Functions).R")
### Modified functions for the Answer to Dochtermann et al. 
source("HoleyLandscape_Newfunctions.R")

######################
### EMPIRICAL DATA ###
######################
##
fi<-readLines("empiricalmatrices.txt") #read the list of empirical matrices
M<-list()
cpt<-0
l2l1<-NULL
for (f in fi) {
	print(cpt)
	cpt<-cpt+1
	M[[cpt]]<-as.matrix(read.table(f, sep=","))
}
names(M)<-fi

M<-M[c(-19,-40,-41)] #we remove 3 problematic matrices (var=0 for some traits)
M<-M[c(-24,-25,-26,-27,-28,-39)] #we remove 6 problematic matrices (l2 abnd l1 are different than those in Dochtermann's table - their value does not correspond to ours)

#compute l2/l1 for the resulting matrices
L2L1.EMPIRICAL.Gmat<-lapply(M, function(x) eigen(x)$values[2]/eigen(x)$values[1])




###################
### SIMULATIONS ###
###################
## We keep the General run conditions ##
#### General Run Conditions ####
#N=7500 #Number of individuals
N=7500
h2=.8 #heritability
#### IVC Model (Wright LS) ####
prop=0.5 #proportion selected


#### SIMULATIONS Gaussian (as in Dochtermann et al. 2023)
#run
RES.GAUSSIAN<-list()
for (i in 1:35) {
	print(i)
	RES.GAUSSIAN[[i]]<-Gaussian_sim(100, ISSmat=NULL, nbtraits=10)
}
#compute l2/l1 (initial = ISS and final = Gmat)
L2L1.GAUSSIAN.ISS<-unlist(lapply(RES.GAUSSIAN, function(x) computel2l1(x$covSTART)))
L2L1.GAUSSIAN.Gmat<-unlist(lapply(RES.GAUSSIAN, function(x) computel2l1(x$covEND)))

#### SIMULATIONS Gaussian starting from empirical-based ISS matrices
#run
RES.GAUSSIAN_empirical<-list()
cpt<-0
for (m in M) {
	cpt<-cpt+1
	print(cpt)
	RES.GAUSSIAN_empirical[[cpt]]<-Gaussian_sim(100, ISSmat=m) #nbtraits is chosen automatically as the number of traits in m by function Gaussian_sim2
}
#compute l2/l1 (initial = ISS and final = Gmat)
L2L1.GAUSSIAN_empirical.ISS<-unlist(lapply(RES.GAUSSIAN_empirical, function(x) computel2l1(x$covSTART))) #Note that this vector is similar to L2L1.EMPIRICAL.Gmat.
L2L1.GAUSSIAN_empirical.Gmat<-unlist(lapply(RES.GAUSSIAN_empirical, function(x) computel2l1(x$covEND)))

#### SIMULATIONS Holey p = 0.5 (as in Dochtermann et al. 2023)
#run
RES.HOLEY<-Holey_sim(35, 100, nbtraits=10)
#compute l2l1 of G-matrix
L2L1.HOLEY.Gmat<-unlist(lapply(RES.HOLEY, function(x) computel2l1(x)))


#############
### PLOTS ###
#############

require(ggplot2)

### PANEL A
## prepare dataframes
DF.A<-data.frame(l2l1start=c(L2L1.GAUSSIAN_empirical.ISS, L2L1.GAUSSIAN.ISS), l2l1end=c(L2L1.GAUSSIAN_empirical.Gmat, L2L1.GAUSSIAN.Gmat), type=rep(c("Gaussian empirical-based","Gaussian Dochtermann"), each=length(L2L1.GAUSSIAN.ISS)))
plot1<-ggplot(DF.A, aes (x=l2l1start,y=l2l1end, color=type)) + geom_abline(intercept = 0, slope = 1, linetype=2)+ geom_point(cex=4, show.legend=FALSE, alpha=0.7) +theme_bw() + labs(x="\u03bb2/\u03bb1",y="\u03bb2/\u03bb1") + scale_color_manual(values = c("#91c5de","orange")) 

### PANEL B
DF.B<-data.frame(class="Gaussian empirical-based", l2l1=L2L1.GAUSSIAN_empirical.Gmat, what="Gaussian empirical-based")
DF.B<-rbind(DF.B, data.frame(class="Gaussian", l2l1=L2L1.GAUSSIAN.Gmat, what="Gaussian"))
DF.B<-rbind(DF.B, data.frame(class="Empirical (Insecta)", l2l1=L2L1.EMPIRICAL.Gmat, what="Empirical (Insecta)"))
DF.B<-rbind(DF.B, data.frame(class="Holey (p = 0.5)", l2l1=L2L1.HOLEY.Gmat, what="Holey (p = 0.5)"))

###CORRECTLY ORDER THE classes for nice plot:
DF.B$class <- factor(DF.B$class , levels=rev(c("Holey (p = 0.5)","Gaussian", "Gaussian empirical-based", "Empirical (Insecta)")))

plot2<-ggplot(DF.B,aes(x=l2l1, y=class)) + 
geom_jitter(aes(color=what), height=0.05, alpha=0.6, size=4, show.legend=FALSE) +
geom_boxplot(width=0.03, col="black", fill="black", outlier.shape = NA, lwd=0.5) + 
stat_summary(
    geom = "point",
    fun = "median",
    col = "black",
    size = 5.5,
    shape = 21,
    fill = c("#ca001e", "orange","#91c5de","#f5f5f5") 
 ) + 
scale_color_manual(values = c("#ca001e", "#91c5de","orange","grey")) + 
labs(x="\u03bb2/\u03bb1", y="") +
theme_bw()

### PLOTS side by side
grid.arrange(plot1, plot2, ncol=2, widths=c(1.5,2))
require(gridExtra)
