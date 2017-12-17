#!/usr/bin/env Rscript
# author: biotang, 2017.09.06
# function: estimate heritability based on a mode of mixture 
#           distribution and its confidence limits (sd) via
#           simulation.
# usage: Rscript estmate_QTL_heritability.R Example_conf.txt

# load R package
options(warn=-1);# close warnings
library(rootSolve);# for solving nonlinear equations
options(warn=0);# open warnings

###########################################
# command arguments and global options
###########################################
# command arguments
opt <- commandArgs(TRUE);
file_conf <- opt[1];

# global options
Times <- 10000;    # simlation times for confidence limits

########################################
# functions
########################################
# read configure file
read_conf <- function(file_conf){
  conf <- list();
  #
  FH <- file(file_conf,"r");
  while(TRUE){
    rline <- readLines(FH,n=1);
    if(length(rline) == 0) break;
    # cat(rline,"\n");
    if (length (grep ("^#", rline, perl=TRUE) ) > 0 ) next;
    if (length (grep ("^\\s", rline, perl=TRUE) ) > 0) next;
    if (nchar (rline) == 0) next;
    # sub
    rline <- sub ("#.*", "", rline, perl=TRUE);
    rline <- sub ("\\s+$", "", rline, perl=TRUE);
    rline <- sub ("\\s*?=\\s*", "=", rline, perl=TRUE);
    # cat(rline,"\n");
    arr <- strsplit (rline, "=")[[1]];
    # cat (arr[1], arr[2], "\n", sep="; ");
    conf[[ arr[1] ]] <- arr[2];
  }
  close(FH);
  return(conf);
}

########################################################
## functions: estimate heritability
########################################################
# estimate broad heritability 
# from additvie effect and dominance effect
# input parameters:
# b is genotype frequence determined by population type
# a is additive effect
# d is dominance effect
estimate_broad_heritability <- function(b, a, d){
  h2 <- (2*b*a^2 + 2*b*(1-2*b)*d^2) / (1 + 2*b*a^2 + 2*b*(1-2*b)*d^2);
  return(h2);
}

# estimate narrow heritability 
# from additvie effect and dominance effect
# narrow heritability is additive effect heritability if no obviously deviation 
estimate_narrow_heritability <- function(b, a, d){
  h2 <- (2*b*a^2) / (1 + 2*b*a^2 + 2*b*(1-2*b)*d^2);
  return(h2);
}

# convert broad heritability to narrow heritability
# rd is rate of dominance, rd = |d/a|
convert_broad_to_narrow <- function(h2B, b, rd){
  h2N <- h2B/(1+(1-2*b)*rd^2);
  return(h2N);
}

# function: convert heritability to additive effect
convert_h2B_to_effect <- function(b, rd, h2B){
  a <- sqrt( h2B / ( (1-h2B) * (2*b + 2*b*rd^2 - 4*b^2*rd^2) ) );
  return(a);
}

# estimate additive heritability by adjusted model with fQ
estimate_additive_heritability <- function(fQ, a){
  h2 <- 4*fQ*(1-fQ)*a^2/(1+4*fQ*(1-fQ)*a^2);
  return(h2);
}

#######################################################################
## Disgn A (FP+OPP): full mode 
## recommended population: filial populatin 
## e.g.: F2, F3, F4....
## Fk population: b = [1-(1/2)^(k-1)]/2, b is frequency of genotype QQ
## OPP: both low pool and high pool is required
## 
#######################################################################
estimate_by_full_mode <- function(b, fQH, fQL, pH, pL){
  # function of full mode (with four equations: Eqs.[4], [5], [6], [7])
  full_mode <- function(x){
    # x: a0, d0, xL, xH
    #F <- rep(NA, length(x));
    F <- numeric(length(x));
    F[1] <- b*pnorm(x[3]-x[1]) + (1-2*b)*pnorm(x[3]-x[2]) + b*pnorm(x[3]+x[1]) - pL;
    F[2] <- 1 - b*pnorm(x[4]-x[1]) - (1-2*b)*pnorm(x[4]-x[2]) - b*pnorm(x[4]+x[1]) - pH;
    F[3] <- 2*b*pnorm(x[3]-x[1]) + (1-2*b)*pnorm(x[3]-x[2]) - 2*fQL*pL;
    F[4] <- 2*b*(1-pnorm(x[4]-x[1])) + (1-2*b)*(1-pnorm(x[4]-x[2])) - 2*fQH*pH;
    F;
  }
  
  # solve nonlinear equations
  # a0, d0, xL, xH
  startx <- c(1, 0, -1, 1);
  #cat("estimating theta.......\n");
  # rootSolve
  result <- multiroot(f = full_mode, start = startx);
  theta  <- result$root;
  precis <- result$estim.precis;
  
  # estimate heritability
  a0 <- theta[1];
  d0 <- theta[2];
  h2B <- estimate_broad_heritability(b, a0, d0); # broad heritability
  h2N <- estimate_narrow_heritability(b, a0, d0);# narrow (additive effect) heritability
  h2D <- h2B - h2N;                              # dominance effect heritability
  #cat("Broad heritability is", sprintf("%.5f", h2B), "\n");
  #cat("Narrow heritability is", sprintf("%.5f", h2N), "\n");
  # dominance degree
  #dd <- abs(d0/a0);
  return( list(h2B=h2B, h2N=h2N, h2D=h2D, a0=a0, d0=d0, precis=precis) );
}

#######################################################################
## Design C (PP+OPP): RIL reduced mode 
## dominance effect is reduced in this mode 
## and this hypothesis is suit with permanent population
## recommended population: permanent population, e.g.: RIL, H, DH
## b = 0.5 is required
## this mode is for APP: at least one of low or high pool
## APP is based on asymmetrical pooling strategy
##
#######################################################################
estimate_by_reduced_mode <- function(b, fQH, fQL, pH, pL){
  # function of reduced mode
  # three equations: Eqs.(11), (12) and (15), 
  # Eqs.(15) is combined from Eqs.(13) and Eqs.(14)
  # b should be 0.5 in this mode
  reduced_mode <- function(x){
    # x: a0, xL, xH
    #F <- rep(NA, length(x));
    F <- numeric(length(x));
    F[1] <- b*pnorm(x[2]-x[1]) + (1-2*b)*pnorm(x[2]) + b*pnorm(x[2]+x[1]) - pL;
    F[2] <- 1 - b*pnorm(x[3]-x[1]) - (1-2*b)*pnorm(x[3]) - b*pnorm(x[3]+x[1]) - pH;
    F[3] <- ( 2*b*(1-pnorm(x[3]-x[1])) + (1-2*b)*(1-pnorm(x[3])) )/(2*pH) - ( 2*b*pnorm(x[2]-x[1]) + (1-2*b)*pnorm(x[2]) )/(2*pL) + fQL - fQH;
    F;
  }
  
  # solve nonlinear equations
  # a0, xL, xH
  startx <- c(1, -1, 1);
  #cat("estimating theta.......\n");
  # rootSolve
  result <- multiroot(f=reduced_mode, start=startx);
  theta  <- result$root;
  precis <- result$estim.precis;
  
  # estimate heritability
  a0 <- theta[1];
  d0 <- 0;
  #h2B <- estimate_broad_heritability(b, a0, d0);
  h2N <- estimate_narrow_heritability(b, a0, d0);
  #fQ <- (fQH + fQL)/2
  #h2A <- estimate_additive_heritability(fQ, a0);
  #cat("Broad heritability is", sprintf("%.5f", h2B), "\n");
  #cat("Narrow heritability is", sprintf("%.5f", h2N), "\n");
  # dominance degree
  #dd <- abs(d0/a0);
  return( list(h2A=h2N, a0=a0, d0=d0, precis=precis) );
}

#############################################################################
## Adjusted Design C (PP+OPP): RIL adjusted reduced mode 
## dominance effect is reduced in this mode
## when fQ severely deviates from 0.5, use this mode
## recommended population: permanent population, e.g.: RIL, H, DH
## b = 0.5 is required
## this mode is for APP: at least one of low or high pool
## APP is based on asymmetrical pooling strategy
##
#############################################################################
estimate_by_adjusted_reduced_mode <- function(b, fQH, fQL, fQ, pH, pL){
  # funciont of reduced mode
  # three Eqs.: Eqs.(18), Eqs.(19) and Eqs.(15)
  # Eqs.(15) is combined from Eqs.(20) and Eqs.(21)
  # real frequency fQ replace expected frequency
  adjusted_reduced_mode <- function(x){
    # x: a0, xL, xH
    #F <- rep(NA, length(x));
    F <- numeric(length(x));
    #fQ <- (fQH + fQL)/2; # if not fQ, this may represent fQ
    F[1] <- fQ*pnorm(x[2]-x[1]) + (1-fQ)*pnorm(x[2]+x[1]) - pL;
    F[2] <- 1 - fQ*pnorm(x[3]-x[1]) - (1-fQ)*pnorm(x[3]+x[1]) - pH;
    F[3] <- fQ*pnorm(x[2]-x[1])/pL - fQ*(1-pnorm(x[3]-x[1]))/pH + fQH - fQL;
    F;
  }
  
  # solve nonlinear equations
  # a0, xL, xH
  startx <- c(1, -1, 1);
  #cat("estimating theta.......\n");
  # rootSolve
  result <- multiroot(f = adjusted_reduced_mode, start = startx);
  theta  <- result$root;
  precis <- result$estim.precis;
  
  # estimate heritability
  a0 <- theta[1];
  d0 <- 0;
  #h2B <- estimate_broad_heritability(b, a0, d0);
  #h2N <- estimate_narrow_heritability(b, a0, d0);
  #fQ <- (fQH + fQL)/2
  h2A <- estimate_additive_heritability(fQ, a0);
  #cat("Broad heritability is", sprintf("%.5f", h2B), "\n");
  #cat("Narrow heritability is", sprintf("%.5f", h2N), "\n");
  # dominance degree
  #dd <- abs(d0/a0);
  return( list(h2A=h2A, a0=a0, d0=d0, precis=precis) );
  # note: only need h2A in this mode, and h2A is h2B if b=0.5
}

#######################################################################
## Design D (PP+APP): RIL asymmetrical mode (Case I: input fQL)
## if want to use Case II (input fQH), just turn fQH to fQL before estimation
## recommended population: permanent population, e.g.: RIL, H, DH
## b = 0.5 is required
## this mode is for APP: at least one of low or high pool
## APP is based on asymmetrical pooling strategy
##
#######################################################################
estimate_by_asymmetrical_mode <- function(b, fQL, pL){
  # function of asymmetrical mode (two equations: eqs.[11] and eqs.[13])
  asymmetrical_mode <- function(x){
    # x: a0, x0
    #F <- rep(NA, length(x));
    F <- numeric(length(x));
    F[1] <- b*pnorm(x[2]-x[1]) + (1-2*b)*pnorm(x[2]) + b*pnorm(x[2]+x[1]) - pL;
    F[2] <- 2*b*pnorm(x[2]-x[1]) + (1-2*b)*pnorm(x[2]) - 2*fQL*pL;
    F;
  }
  
  # solve nonlinear equations
  # a0, xL
  startx <- c(1, -1);
  #cat("estimating theta.......\n");
  # rootSolve
  result <- multiroot(f = asymmetrical_mode, start = startx);
  theta  <- result$root;
  precis <- result$estim.precis;
  
  # estimate heritability
  a0 <- theta[1];
  d0 <- 0;
  #h2B <- estimate_broad_heritability(b, a0, d0); # broad heritability
  h2N <- estimate_narrow_heritability(b, a0, d0);# narrow (additive effect) heritability
  #h2A <- estimate_additive_heritability(fQU, a0);
  #cat("Broad heritability is", sprintf("%.5f", h2B), "\n");
  #cat("Narrow heritability is", sprintf("%.5f", h2N), "\n");
  # dominance degree
  #dd <- abs(d0/a0);
  return( list(h2A=h2N, a0=a0, d0=d0, precis=precis) );
}

########################################################################
## Adjusted Design D (PP+APP): RIL adjusted asymmetrical mode 
## (Case I: input fQL and fQ)
## recommended population: permanent population, e.g.: RIL, H, DH
## b = 0.5 is required
## this mode is for APP: at least one of low or high pool
## APP is based on asymmetrical pooling strategy
##
########################################################################
estimate_by_adjusted_asymmetrical_mode <- function(b, fQL, fQ, pL){
  # function of adjusted asymmetrical mode
  # two Eqs.: Eqs.(18) and Eps.(20)
  adjusted_asymmetrical_mode <- function(x){
    # x: a0, x0
    #F <- rep(NA, length(x));
    F <- numeric(length(x));
    F[1] <- fQ*pnorm(x[2]-x[1]) + (1-fQ)*pnorm(x[2]+x[1]) - pL;
    F[2] <- fQ*pnorm(x[2]-x[1]) - fQL*pL;
    F;
  }
  
  # solve nonlinear equations
  # a0, xL
  startx <- c(1, -1);
  #cat("estimating theta.......\n");
  # rootSolve
  result <- multiroot(f = adjusted_asymmetrical_mode, start = startx);
  theta  <- result$root;
  precis <- result$estim.precis;
  
  # estimate heritability
  a0 <- theta[1];
  d0 <- 0;
  #h2B <- estimate_broad_heritability(b, a0, d0);
  #h2N <- estimate_narrow_heritability(b, a0, d0);
  h2A <- estimate_additive_heritability(fQ, a0);
  #cat("Broad heritability is", sprintf("%.5f", h2B), "\n");
  #cat("Narrow heritability is", sprintf("%.5f", h2N), "\n");
  # dominance degree
  #dd <- abs(d0/a0);
  return( list(h2A=h2A, a0=a0, d0=d0, precis=precis) );
}

########################################################
# simulation of Bulked Allele Frequency
########################################################
# esimate frequency of genotype
estimate_genotype_frequency <- function(gt){
  n <- length(gt);
  f1 <- length(gt[gt==1])/n;
  f2 <- length(gt[gt==0])/n;
  f3 <- length(gt[gt==-1])/n;
  return( list(f1=f1, f2=f2, f3=f3));
}

# simulate once
simulate_bulked_AF_once <- function(n, mL, mH, fL, fH, a, d){
  # first, genotype simulation
  # genotype value:
  # QQ, Qq, qq
  # 1,  0,  -1
  rd <- runif(n);              # sample size is n
  gt <- rep(0, n);             # default is 0 
  gt <- ifelse(rd>=fH, 1, gt); # QQ for high trait
  gt <- ifelse(rd<=fL, -1, gt);# qq for low trait

  # second, phenotype simulation
  # simulating by variance of error
  #zz <- 1-abs(gt);
  #yy <- 50 + matrix(gt, ncol=1) %*% a + matrix(zz, ncol=1) %*% d + rnorm(n, sd=sqrt(ve)); 
  #pt <- yy[,1];# convert to a vector
  # simulating by mixture distribution
  x <- 50 + rnorm(n);
  pt <- x + gt*a + (1-abs(gt))*d
  #pt <- x;
  #pt <- ifelse(gt ==  1, x+a, pt);# fH QQ, High-trait
  #pt <- ifelse(gt ==  0, x+d, pt);# fM Qq, Mid-trait
  #pt <- ifelse(gt == -1, x-a, pt);# fL qq, Low-trait
  #hist(pt);
  
  # third, generate simulated allele frequency
  # unselected pools
  U1 <- gt[1:mL];        # unselected pool 1 by random order, pool size = mL
  U2 <- gt[(n-mH+1):n];  # unselected pool 2 by random order, pool size = mH
  # sort genotype by order of phenotype value
  gt <- gt[order(pt)];
  HP <- gt[(n-mH+1):n];  # High pool with high trait-value
  LP <- gt[1:mL];        # low pool with low trait-value
  
  # estimate genotype frequency
  U1_geno <- estimate_genotype_frequency(U1);
  U2_geno <- estimate_genotype_frequency(U2);
  HP_geno <- estimate_genotype_frequency(HP);
  LP_geno <- estimate_genotype_frequency(LP);
  TT_geno <- estimate_genotype_frequency(gt);
  
  # simulated allele frequency
  U1_AF <- U1_geno[[1]] + 0.5*U1_geno[[2]];# AF in unselected pool one
  U2_AF <- U2_geno[[1]] + 0.5*U2_geno[[2]];# AF in unselected pool two
  HP_AF <- HP_geno[[1]] + 0.5*HP_geno[[2]];# AF in high pool
  LP_AF <- LP_geno[[1]] + 0.5*LP_geno[[2]];# AF in low pool
  TT_AF <- TT_geno[[1]] + 0.5*TT_geno[[2]];# total AF
  #
  return(list(HAF=HP_AF, LAF=LP_AF, U1AF=U1_AF, U2AF=U2_AF, TAF=TT_AF));
}

########################################################
# functions
########################################################
display_estimation <- function(h2){
  # broad
  avgh2 <- NA;
  avgh2 <- mean(h2, na.rm = TRUE);
  if(!is.na(avgh2) ){
    cat("mean of simulated heritability:", avgh2, "\n");
    cat("SD of simulations:", sd(h2, na.rm = TRUE), "\n");
    #cat("2.5% and 97.5% quantile:", quantile(h2, probs = 2.5/100, na.rm = TRUE), "~", quantile(h2, probs = 97.5/100, na.rm = TRUE), "\n");
    #cat("5% and 95% quantile:", quantile(h2, probs = 5/100, na.rm = TRUE), "~", quantile(h2, probs = 95/100, na.rm = TRUE), "\n");
    #cat("10% and 90% quantile:", quantile(h2, probs = 10/100, na.rm = TRUE), "~", quantile(h2, probs = 90/100, na.rm = TRUE), "\n");
    #cat("15% and 85% quantile:", quantile(h2, probs = 15/100, na.rm = TRUE), "~", quantile(h2, probs = 85/100, na.rm = TRUE), "\n");
    #cat("20% and 80% quantile:", quantile(h2, probs = 20/100, na.rm = TRUE), "~", quantile(h2, probs = 80/100, na.rm = TRUE), "\n");
    #cat("25% and 75% quantile:", quantile(h2, probs = 25/100, na.rm = TRUE), "~", quantile(h2, probs = 75/100, na.rm = TRUE), "\n");
  }
  cat("\n");
}

########################################################
## main: structure of files
########################################################
# examples
#file_conf <- "examples/Example1.rice_Chr3.conf.txt";
#file_conf <- "examples/Example2.rice_Chr6.conf.txt";
#file_conf <- "examples/Example3.yeast_Chr13.conf.txt";
#file_conf <- "examples/Example4.yeast_Chr15.conf.txt";

# prefix of output
outPrefix <- sub(pattern = "(.txt$)|(.xls$)|(.tab$)|(.tsv$)", perl = TRUE, replacement = "", file_conf);
outPrefix <- sub(pattern = "[._-]?+conf$", perl = TRUE, replacement = "", outPrefix);

# output file of simulations
file_out <- paste(outPrefix, ".simulation.txt", sep = "");

# log file of warnings
warn_log  <- paste(outPrefix, ".warnings.log", sep="");

########################################################
## main: process parameters from configure file
########################################################
# read configure file
cat("read global parameters from", file_conf, "\n");
conf <- read_conf(file_conf);
for(name in names(conf) ){
  cat("parameter ", name, "=", conf[[name]],  "\n", sep="");
}
cat("\n");

# check required parameters
#!!!not finished this part!!!
if(is.null(conf$b)) {cat("Error: parameter b is not defined\n");quit();}
if(is.null(conf$m1)) {cat("Error: parameter m1 is not defined\n");quit();}
if(is.null(conf$n)) {cat("Error: parameter n is not defined\n");quit();}
if(is.null(conf$AF1)) {cat("Error: parameter AF1 is not defined\n");quit();}

# raw global options
b   <- NA; # frequency of homozygous (AA), required
m1  <- NA; # selected indivaiduals in pool 1, required
m2  <- NA; # selected indivaiduals in pool 2, not required
n   <- NA; # population size, required
p1  <- NA; # selective press/rate of pool 1
p2  <- NA; # selective press/rate of pool 2
af1 <- NA;# allele frequency of pool 1, required
af2 <- NA;# allele frequency of pool 2, not required
afU <- NA;# allele frequency of unselected pool or total population, not required

# global options from configure file
if(!is.null(conf$b))  b   <- as.numeric(conf$b);  # frequency of homozygous (AA), required
if(!is.null(conf$m1)) m1  <- as.numeric(conf$m1); # selected indivaiduals in pool 1, required
if(!is.null(conf$m2)) m2  <- as.numeric(conf$m2); # selected indivaiduals in pool 2, not required
if(!is.null(conf$n))  n   <- as.numeric(conf$n);  # population size, required
p1  <- m1/n # selective press/rate of pool 1
p2  <- m2/n # selective press/rate of pool 2
if(!is.null(conf$AF1)) af1 <- as.numeric(conf$AF1);# allele frequency of pool 1, required
if(!is.null(conf$AF2)) af2 <- as.numeric(conf$AF2);# allele frequency of pool 2, not required
if(!is.null(conf$AFU)) afU <- as.numeric(conf$AFU);# allele frequency of unselected pool or total population, not required

# assign to low pool and high pool
fQH <- af1;
fQL <- af2;
fQU <- afU;
pH  <- p1;
pL  <- p2;
mH  <- m1;
mL  <- m2;
# convertion: exchange when both opposite pools
if(!is.na(fQL) & fQL > fQH){
  fQH <- af2;
  fQL <- af1;
  pH  <- p2;
  pL  <- p1;
  mH  <- m2;
  mL  <- m1;
}
# convertion when only asymmetrical pool
if(is.na(fQL)){
  fQL <- af1;
  fQH <- NA;
  pL  <- p1;
  pH  <- pL;
  mL  <- m1;
  mH  <- mL;
  if(fQL > 0.5){
    fQL <- 1 - fQL;
    # fQU is defined when fQU severely deviated from 0.5
    if(!is.na(fQU) ) fQU <- 1 - fQU;
  }
}

# display conditions
cat("calculation conditions:\n");
cat("frequency of homozygous genotype, b =", b, "\n");  # frequency of homozygous (AA) genotype
cat("selected indivaiduals in low pool, mL =", mL, "\n"); # selected indivaiduals in low pool
cat("selected indivaiduals in high pool, mH =", mH, "\n"); # selected indivaiduals in high pool
cat("population size, n =", n, "\n");  # population size
cat("selective press/rate of low pool, pL =", pL, "\n"); # selective press/rate of low pool
cat("selective press/rate of high pool, pH =", pH, "\n"); # selective press/rate of high pool
cat("allele frequency of low pool, fQL =", fQL, "\n");# allele frequency of low pool
cat("allele frequency of high pool, fQH =", fQH, "\n");# allele frequency of high pool
cat("allele frequency of unselected pool, fQU =", fQU, "\n");# allele frequency of unselected pool or total population
cat("\n");

########################################################
## main: estimate heritability
########################################################
# estimating by different modes
if(b<0.5)                               h2_result <- estimate_by_full_mode(b, fQH, fQL, pH, pL);
if(b==0.5 & !is.na(fQH) & is.na(fQU) )  h2_result <- estimate_by_reduced_mode(b, fQH, fQL, pH, pL);
if(b==0.5 & !is.na(fQH) & !is.na(fQU) ) h2_result <- estimate_by_adjusted_reduced_mode(b, fQH, fQL, fQU, pH, pL);
if(b==0.5 & is.na(fQH) & is.na(fQU) )   h2_result <- estimate_by_asymmetrical_mode(b, fQL, pL);
if(b==0.5 & is.na(fQH) & !is.na(fQU) )  h2_result <- estimate_by_adjusted_asymmetrical_mode(b, fQL, fQU, pL);

# list of explanation
h2_explain <- list(
  h2B = "Broad (total) heratibility",
  h2N = "Narrow (additive effect) heratibility",
  h2A = "Additive effect heratibility",
  h2D = "Dominance effect heratibility",
  a0  = "Additive effect",
  d0  = "Dominance effect",
  precis = "precision of equations solving"
);

# display the result of estimation
cat("Estimation of heritabilty:\n")
for(name in names(h2_result) ){
  cat(h2_explain[[name]], ": ", h2_result[[name]],  "\n", sep="");
}
cat("\n")

########################################################
## main: estimate confidence limits of heritability
## based on simulation
########################################################
# frequency of genotype
# breakpoint of segregation
fL <- b;      # genotype should be -1 (qq) ~ low trait
#fM <- 1-2*b; # genotype should be 0  (Qq) ~ middle trait
fH <- 1-b;    # genotype should be 1  (QQ) ~ high trait

# for adjusted mode
if(b == 0.5 & !is.na(fQU) ){
  fL <- fQU;
  fH <- 1 - fQU;
}

# estimated effect (real value)
a <- h2_result$a0; # additive effect
d <- h2_result$d0; # dominance effect

# simulation conditions
cat("Simulation conditions:\n");
cat("genotype frequency (b):", b, "\n");
cat("population size:", n, "\n");
cat("selected individuals:", mL, "\n");
cat("selection rate:", pL, "\n");
cat("selected individuals:", mH, "\n");
cat("selection rate:", pH, "\n");
cat("a0:", a, "\n");
cat("d0:", d, "\n");
cat("simulation times:", Times, "\n");
cat("\n");

# simulating bulked allele frequency
haf  <- c();
laf  <- c();
uaf  <- c();
u2af <- c();
taf  <- c();
for(i in 1:Times){
  # simulating
  ss  <- simulate_bulked_AF_once(n, mL, mH, fL, fH, a, d);
  # simulation data
  haf <- c(haf, ss$HAF);
  laf <- c(laf, ss$LAF);
  uaf <- c(uaf, ss$U1AF);
  u2af<- c(u2af,ss$U2AF);
  taf <- c(taf, ss$TAF);
}

# overview statistics of simulated bulked allele frequency 
cat("overview of simulated bulked allele frequency\n")
cat("simulated data (display average value): \n");
cat("Total AF =", mean(taf), "\n");
#cat("2.5% and 97.5% quantile:", quantile(taf, probs = 2.5/100), "~", quantile(taf, probs = 97.5/100), "\n");
cat("fQH =", mean(haf), "\n");
#cat("25% and 75% quantile:", quantile(haf, probs = 25/100), "~", quantile(haf, probs = 75/100), "\n");
cat("fQL =", mean(laf), "\n");
#cat("25% and 75% quantile:", quantile(laf, probs = 25/100), "~", quantile(laf, probs = 75/100), "\n");
#afd  <- (haf - laf);
#cat("AFD = fQH - fQL =", mean(afd), "\n");
#cat("25% and 75% quantile:", quantile(afd, probs = 25/100), "~", quantile(afd, probs = 75/100), "\n");
cat("fQu1 =", mean(uaf), "\n");
#cat("2.5% and 97.5% quantile:", quantile(uaf, probs = 2.5/100), "~", quantile(uaf, probs = 97.5/100), "\n");
cat("fQu2 =", mean(u2af), "\n");
#cat("2.5% and 97.5% quantile:", quantile(u2af, probs = 2.5/100), "~", quantile(u2af, probs = 97.5/100), "\n");
#afd  <- (uaf - laf);
#cat("AFD = fQU1 - fQL =", mean(afd), "\n");
#cat("25% and 75% quantile:", quantile(afd, probs = 25/100), "~", quantile(afd, probs = 75/100), "\n");
#afd  <- (u2af - laf);
#cat("AFD = fQU2 - fQL =", mean(afd), "\n");
#cat("25% and 75% quantile:", quantile(afd, probs = 25/100), "~", quantile(afd, probs = 75/100), "\n");
#afd  <- (haf - uaf);
#cat("AFD = fQH - fQU1 =", mean(afd), "\n");
#cat("25% and 75% quantile:", quantile(afd, probs = 25/100), "~", quantile(afd, probs = 75/100), "\n");
#afd  <- (haf - u2af);
#cat("AFD = fQH - fQU2 =", mean(afd), "\n");
#cat("25% and 75% quantile:", quantile(afd, probs = 25/100), "~", quantile(afd, probs = 75/100), "\n");
#afd  <- abs(uaf - u2af);
#cat("AFD = abs(fQU1 - fQU2) =", mean(afd), "\n");
#cat("95% quantile:", quantile(afd, probs = 95/100), "\n");
cat("\n");

# when solving, make it sink to log file.
sink(file=warn_log);
# close warning temporarily
options(warn=-1);

# estimating heritability for simulations
# estimations
fh2B <- c();
fh2N <- c();
fh2A <- c();
fh2D <- c();

# marker
precis_mk <- "NO";
# estimating
for(i in 1:Times){
  # check deviation
  if(haf[i] <= 0.5 | laf[i] >= 0.5){ 
    fh2B <- c(fh2B, NA);
    fh2N <- c(fh2N, NA);
    fh2A <- c(fh2A, NA);
    fh2D <- c(fh2D, NA);    
    next
  }
  
  # estimating by different modes
  if(b<0.5)                 fhh <- estimate_by_full_mode(b, haf[i], laf[i], pH, pL);
  if(b==0.5 & !is.na(fQH) ) fhh <- estimate_by_adjusted_reduced_mode(b, haf[i], laf[i], taf[i], pH, pL);
  if(b==0.5 &  is.na(fQH) ) fhh <- estimate_by_adjusted_asymmetrical_mode(b, 1-haf[i], 1-taf[i], pL);
  # save
  if(fhh$precis > 0.001 ){
    precis_mk <- "YES";
    fhh$h2B <- NA;
    fhh$h2N <- NA;
    fhh$h2A <- NA;
    fhh$h2D <- NA;
  }
  if(is.null(fhh$h2B) ) fhh$h2B <- NA;
  if(is.null(fhh$h2N) ) fhh$h2N <- NA;
  if(is.null(fhh$h2A) ) fhh$h2A <- NA;
  if(is.null(fhh$h2D) ) fhh$h2D <- NA;
  fh2B <- c(fh2B, fhh$h2B);
  fh2N <- c(fh2N, fhh$h2N);
  fh2A <- c(fh2A, fhh$h2A);
  fh2D <- c(fh2D, fhh$h2D);
}

# sink to stdout
sink();
options(warn=0);

# display results
cat("overview of heritability estimation from simulated data:\n");
if(length( which( is.na(fh2B) ) ) < Times ){ 
  cat("Broad (total) heritability.\n")
  display_estimation(fh2B);
}
if(length( which( is.na(fh2N) ) ) < Times ){ 
  cat("Narrow (additive effect) heritability.\n")
  display_estimation(fh2N);
}
if(length( which( is.na(fh2A) ) ) < Times ){ 
  cat("Additive effect heritability.\n")
  display_estimation(fh2A);
}
if(length( which( is.na(fh2D) ) ) < Times ){ 
  cat("Domiance effect heritability.\n")
  display_estimation(fh2D);
}

# write out
#out <- data.frame(taf, haf, laf, sprintf("%.4f", haf - laf), uaf, u2af, sprintf("%.4f", abs(uaf - u2af)), fh2B, fh2N, fh2A, fh2D);
#colnames(out) <- c("Total AF", "fQH", "fQL", "fQH - fQL", "fQU1", "fQU2", "abs fQU1 - fQU2", "h2B", "h2N", "h2A", "h2D");
out <- data.frame(taf, haf, laf, uaf, u2af, fh2B, fh2N, fh2A, fh2D);
colnames(out) <- c("Total_AF", "fQH", "fQL", "fQU1", "fQU2", "h2B", "h2N", "h2A", "h2D");
write.table(out, file_out, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t", append=FALSE)
