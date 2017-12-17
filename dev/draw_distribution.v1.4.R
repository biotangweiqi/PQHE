library(ggplot2)
library(RColorBrewer)
library(reshape2)

# global
a0 <- 1
d0 <- 0.6
b  <- 0.25
pL <- 0.1
pH <- 1 - 0.1

# global options for calculating cut point
STEP  <- 1;
THETA <- 1E-10;

########################################################
## functions: bi-search cut point
########################################################
#
cal_selection_proportion <- function(b, a, d, x0){
  p <- b*pnorm(x0-a) + (1-2*b)*pnorm(x0-d) + b*pnorm(x0+a);
  return(p);
}

#
bisearch_cut_point <- function(b, a, d, p0){
  L <- -2;
  R <- 2;
  #
  pl <- cal_selection_proportion(b, a, d, L);
  while( (pl-p0) > 0 ){
    L <- L-STEP;
    pl<- cal_selection_proportion(b, a, d, L);
  }
  
  #
  pr <- cal_selection_proportion(b, a, d, R);
  while( (pr-p0) < 0 ){
    R <- R+STEP;
    pr<- cal_selection_proportion(b, a, d, R);
  }
  #
  x0<- (L+R)/2;
  p <- cal_selection_proportion(b, a, d, x0);
  #print(paste(L,pl,R,pr,x0,p));
  #
  while( abs(p-p0)>THETA*p0 ){
    if( (p-p0)<0 ) L<-x0;
    if( (p-p0)>0 ) R<-x0;
    x0 <- (L+R)/2;
    p <- cal_selection_proportion(b, a, d, x0);
    #print(paste(L,pl,R,pr,x0,p));
  }
  return(x0);
}

cal_y1 <- function(x, a0, d0){
  b/(sqrt(2*pi)*exp((x+a0)^2/2) )
}

cal_y2 <- function(x, a0, d0){
  (1-2*b)/(sqrt(2*pi)*exp((x-d0)^2/2) )
}

cal_y3 <- function(x, a0, d0){
  b/(sqrt(2*pi)*exp((x-a0)^2/2) )
}

cal_fx <- function(x, a0, d0){
  cal_y1(x, a0, d0) + cal_y2(x, a0, d0) + cal_y3(x, a0, d0)
}
########################################################
## main
########################################################
#
x <- seq(-4, 4, by = 0.001)

#
#y1 <- b/(sqrt(2*pi)*exp((x+a0)^2/2) )       # qq
#y2 <- (1-2*b)/(sqrt(2*pi)*exp((x-d0)^2/2) ) # Qq
#y3 <- b/(sqrt(2*pi)*exp((x-a0)^2/2) )       # QQ
y1 <- cal_y1(x, a0, d0)
y2 <- cal_y2(x, a0, d0)
y3 <- cal_y3(x, a0, d0)

#
dat <- data.frame(y1, y2, y3)
dat <- cbind(x, dat, fx = rowSums(dat) )

xL <- bisearch_cut_point(b, a0, d0, pL)
xH <- bisearch_cut_point(b, a0, d0, pH)
sum(dat$fx[dat$x<xL])/sum(dat$fx)
sum(dat$fx[dat$x>xH])/sum(dat$fx)

yL <- cal_y1(xL,a0,d0) + cal_y2(xL,a0,d0) + cal_y3(xL,a0,d0)
yH <- cal_y1(xH,a0,d0) + cal_y2(xH,a0,d0) + cal_y3(xH,a0,d0)

#
#plot(x = dat$x, y = dat$fx,
#     type = "l", lwd = 2, col = "red",
#     xlab = "x",
#     ylab = "f(x)",
#     panel.first = grid()
#)
#points(x = dat$x, y = dat$y1, type = "l", lwd = 2, lty = 2, col = "blue")
#points(x = dat$x, y = dat$y2, type = "l", lwd = 2, lty = 2, col = "green")
#points(x = dat$x, y = dat$y3, type = "l", lwd = 2, lty = 2, col = "orange")
#abline(v = xL)
#abline(v = xH)

#
#CLR <- c("firebrick1", "forestgreen", "darkviolet", "dodgerblue4")
#p <- ggplot(data = dat, aes(x = x, y = fx) ) + geom_line(color = "firebrick1", lwd = 1.5)
#p <- p + geom_line(data = dat, aes(x = x, y = y1), color = "forestgreen", lwd = 1.2, lty = 2)
#p <- p + geom_line(data = dat, aes(x = x, y = y2), color = "darkviolet", lwd = 1.2, lty = 2)
#p <- p + geom_line(data = dat, aes(x = x, y = y3), color = "dodgerblue4", lwd = 1.2, lty = 2)
#p

# reshape data
melt_data <- melt(
  data = dat, 
  id.vars = "x", # "Standardized trait (x)"
  variable.name = "distribution", # , 
  value.name = "y", # "Probability density (y)", 
  measure.vars = c("y1", "y2", "y3", "fx")
)


p <- ggplot(data = melt_data, aes(x, y) ) + geom_line(size = 1.2, aes(linetype = distribution) )
#
#p <- p + theme_minimal()

# title and labels
#p <- p + ggtitle("Mixture distribution") + theme(plot.title = element_text(hjust = 0.5, size = rel(1.8)) )
p <- p + xlab(expression (paste("Standardized trait (", italic("x"),")", sep = "") ) )
p <- p + ylab(expression (paste("Probability density ",italic(f), "(", italic(x), ")", sep = "") ) )
p <- p + theme(axis.title = element_text(size = rel(1.5) ), axis.text = element_text(size = rel(1.2)) )

# text
#p <- p + annotate("text", x = min(dat$x), y = max(dat$fx), hjust = 0, vjust = 1, label = "Set:")
#eqn_b  <- as.character( as.expression(substitute(italic("b")==b, list(b = b) ) ) )
#p <- p + annotate("text", x = max(dat$x), y = max(dat$fx)*0.92, hjust = 1, vjust = 1, label = eqn_b, parse = TRUE )
#eqn_a0 <- as.character( as.expression(substitute(italic("a")[0]==a0, list(a0 = a0) ) ) )
#p <- p + annotate("text", x = max(dat$x), y = max(dat$fx)*0.84, hjust = 1, vjust = 1, label = eqn_a0, parse = TRUE )
#eqn_d0 <- as.character( as.expression(substitute(italic("d")[0]==d0, list(d0 = d0) ) ) )
#p <- p + annotate("text", x = max(dat$x), y = max(dat$fx)*0.76, hjust = 1, vjust = 1, label = eqn_d0, parse = TRUE )

# legend
lab.y3 <- expression(paste(pi[QQ], italic(f)[QQ], "(",italic(x),")=", italic(b), italic(phi), "(", italic(x), "-", italic(a)[0], ")", sep = "") )
lab.y2 <- expression(paste(pi[Qq], italic(f)[Qq], "(",italic(x),")=(", "1-2", italic(b), ")", italic(phi), "(", italic(x), "-", italic(d)[0], ")", sep = "") )
lab.y1 <- expression(paste(pi[qq], italic(f)[qq], "(",italic(x),")=", italic(b), italic(phi), "(", italic(x), "+", italic(a)[0], ")", sep = "") )
lab.fx <- expression(paste(italic(f), "(",italic(x),")=", pi[QQ], italic(f)[QQ], "(",italic(x),")+", pi[Qq], italic(f)[Qq], "(",italic(x),")+", pi[qq], italic(f)[qq], "(",italic(x),")", sep = "") )
#p <- p + scale_colour_discrete(name = "Distribution curve:", breaks = c("fx", "y3", "y2", "y1"), labels = c(lab.fx, lab.y3, lab.y2, lab.y1) )
#p <- p + scale_colour_discrete(
#  name = "Distribution curve:",  
#  breaks = c("fx", "y3", "y2", "y1"), 
#  labels = c(lab.fx, lab.y3, lab.y2, lab.y1) 
#)
p <- p + scale_linetype_discrete(
  name = NULL,
  limits = c("fx", "y3", "y2", "y1"),
  breaks = c("fx", "y3", "y2", "y1"), 
  labels = c(lab.fx, lab.y3, lab.y2, lab.y1)
)
p <- p + geom_point(data = data.frame(xx=c(0, 0, 0), parameters=c("bb", "aa", "dd"), yy=c(0, 0, 0)), aes(x=xx, y=yy, fill = parameters), alpha = 0.001 ) 
eqn_b  <- substitute(paste(italic("b")=="", b, sep = ""), list(b=b) )
eqn_a0 <- substitute(paste(italic("a")[0]=="", a0, sep = ""), list(a0=a0) )
eqn_d0 <- substitute(paste(italic("d")[0]=="", d0, sep = ""), list(d0=d0) )
p <- p + scale_fill_discrete(name = NULL, breaks = c("bb", "aa", "dd"), labels = c(eqn_b, eqn_a0, eqn_d0) )
# theme of legend
#p <- p + theme(legend.position = c(0.01,0.98), legend.justification = c(0,1) )
#p <- p + theme(legend.background = element_rect(fill = "white", color = "black") )
p <- p + theme(legend.text = element_text(size = rel(1.36) ), legend.title = element_text(size = rel(1.5) ) )
p <- p + theme(legend.key = element_rect(fill = NA) )
p <- p + theme(legend.key.height = unit(1.2, "cm"), legend.text.align = 0)

#
eqn_xL <- as.character( as.expression(substitute(italic("x")["L"]==xL, list(xL = sprintf("%.3f",xL) ) ) ) )
p <- p + annotate("text", x = xL, y = yL, hjust = 1.1, vjust = 0, label = eqn_xL, parse = TRUE, size = 6 )
eqn_xH <- as.character( as.expression(substitute(italic("x")["H"]==xH, list(xH = sprintf("%.3f",xH) ) ) ) )
p <- p + annotate("text", x = xH, y = yH, hjust = -0.1, vjust = 0, label = eqn_xH, parse = TRUE, size = 6 )
# draw arrow
p <- p + annotate("segment", x = xL, xend = xL, y = yL, yend = 0, arrow = arrow(length = unit(0.3, "cm") ), color = "black", size = 0.5)
p <- p + annotate("segment", x = xH, xend = xH, y = yH, yend = 0, arrow = arrow(length = unit(0.3, "cm") ), color = "black", size = 0.5)

# draw area
data_xL <- data.frame(x = seq(min(dat$x), xL, by = 0.001))
data_xL$y <- cal_fx(data_xL$x, a0, d0)
#p <- p + geom_area(data = data_xL, aes(x, y), fill = "firebrick1", alpha = 0.3)
p <- p + geom_area(data = data_xL, aes(x, y), fill = "grey", alpha = 0.6)

data_xH <- data.frame(x = seq(xH, max(dat$x), by = 0.001))
data_xH$y <- cal_fx(data_xH$x, a0, d0)
#p <- p + geom_area(data = data_xH, aes(x, y), fill = "steelblue2", alpha = 0.3)
p <- p + geom_area(data = data_xH, aes(x, y), fill = "grey", alpha = 0.6)
#
#fx_area <- function(min, max, a0, d0){
#  function(x){
#    y <- cal_fx(xx, a0, d0)
#    y[x<min | x>max] <- NA
#    return(y)
#  }
#}
#p <- p + stat_function(data = ,fun = fx_area(min(dat$x), xL, a0, d0), geom = "area", fill = "firebrick", alpha = 0.25 )
#p <- p + stat_function(fun = fx_area(xH, max(dat$x), a0, d0), geom = "area", fill = "steelblue", alpha = 0.25 )
#p <- p + stat_function(fun = fx_area(xL, xH, a0, d0), geom = "area", fill = "forestgreen", alpha = 0.1)

#
xpL <- xL - (xL-min(dat$x))/4
ypL <- cal_y1(xpL, a0, d0) + cal_y2(xpL, a0, d0) + cal_y3(xpL, a0, d0)
xpH <- xH + (max(dat$x)-xH)/4
ypH <- cal_y1(xpH, a0, d0) + cal_y2(xpH, a0, d0) + cal_y3(xpH, a0, d0)
#
eqn_pL <- as.character( as.expression(substitute(italic("p")["L"]==pL, list(pL = pL ) ) ) )
eqn_pH <- as.character( as.expression(substitute(italic("p")["H"]==pH, list(pH = 1-pH ) ) ) )
p <- p + annotate("text", x = xpL-0.05, y = ypL+0.008, hjust = 1, vjust = 0, label = eqn_pL, parse = TRUE, size = 6 )
p <- p + annotate("text", x = xpH+0.05, y = ypH+0.008, hjust = 0, vjust = 0, label = eqn_pH, parse = TRUE, size = 6 )
p <- p + annotate("segment", x = xpL-0.05, xend = xpL+0.1, y = ypL+0.005, yend = ypL-0.01, arrow = arrow(length = unit(0.3, "cm") ), color = "black", size = 1, alpha = 0.6)
p <- p + annotate("segment", x = xpH+0.05, xend = xpH-0.1, y = ypH+0.005, yend = ypH-0.01, arrow = arrow(length = unit(0.3, "cm") ), color = "black", size = 1, alpha = 0.6)
#p <- p + annotate("text", x = min(dat$x), y = max(dat$fx)*0.68, hjust = 0, vjust = 1, label = eqn_pL, parse = TRUE )
#p <- p + annotate("text", x = min(dat$x), y = max(dat$fx)*0.60, hjust = 0, vjust = 1, label = eqn_pH, parse = TRUE )

#
#p <- p + annotate("text", x = min(dat$x)+1, y = max(dat$fx), hjust = 0, vjust = 1, label = "Solve:" )
#p <- p + annotate("text", x = min(dat$x)+1, y = max(dat$fx)*0.92, hjust = 0, vjust = 1, label = eqn_xL, parse = TRUE )
#p <- p + annotate("text", x = min(dat$x)+1, y = max(dat$fx)*0.84, hjust = 0, vjust = 1, label = eqn_xH, parse = TRUE )

# save
file_tiff <- paste("mixture_distribution_b_is_", b, "_a0_is_", a0, "_d0_is_", d0, ".tiff", sep = "")
ggsave(file_tiff, width = 12, height = 5, device = "tiff", dpi = 500, compression = "lzw")



