
data1 <- readRDS("sigmaSquare2_N16383_log2Gen21_alpha1.6.rds")
data2 <- readRDS("sigmaSquare2_N32767_log2Gen21_alpha1.6.rds")
data3 <- readRDS("sigmaSquare2_N65535_log2Gen21_alpha1.6.rds")
data4 <- readRDS("sigmaSquare2_N131071_log2Gen21_alpha1.6.rds")
data5 <- readRDS("sigmaSquare2_N262143_log2Gen21_alpha1.6.rds")
data6 <- readRDS("sigmaSquare2_N524287_log2Gen21_alpha1.6.rds")


#write.table(data1, file = "sigmaSquare2_N16383_log2Gen21_alpha1.6.txt", col.names = FALSE, row.names = FALSE)
#write.table(data2, file = "sigmaSquare2_N32767_log2Gen21_alpha1.6.txt", col.names = FALSE, row.names = FALSE)
#write.table(data3, file = "sigmaSquare2_N65535_log2Gen21_alpha1.6.txt", col.names = FALSE, row.names = FALSE)
#write.table(data4, file = "sigmaSquare2_N131071_log2Gen21_alpha1.6.txt", col.names = FALSE, row.names = FALSE)
#write.table(data5, file = "sigmaSquare2_N262143_log2Gen21_alpha1.6.txt", col.names = FALSE, row.names = FALSE)
#write.table(data6, file = "sigmaSquare2_N524287_log2Gen21_alpha1.6.txt", col.names = FALSE, row.names = FALSE)

sigmaSquare <- list(data1,data2,data3,data4,data5,data6)

mlty <- length(sigmaSquare):1
mcol <- rainbow(21+1)
g1 <- 0
nbSite = 250 #100000

## plot of the variance sigma2 vs the sites (index)
## do not show all of the sites
## plot one line per N
## FIGURE 2 of the preprint
pdf("variance_vs_site_zoom.pdf", width = 14, height = 14)

for(i in 1:length(sigmaSquare)) {

    N <- dim(sigmaSquare[[i]])[1] -1
    log2Generations <- dim(sigmaSquare[[i]])[2]-1

    odds <- seq(2,N,by=2)  # odd sites because origin is at index = 1

    # First, the qualitative shape (no rescaling)
    for(log2gen in g1:log2Generations) {
       if(i == 1 & log2gen == g1) plot(odds,sigmaSquare[[i]][odds, 1+log2gen] , xlab="site_index", ylab="variance", xlim=c(0,nbSite), ylim=c(0,1), type="l", col=mcol[log2gen+1], lty = mlty[i])
       points(odds, sigmaSquare[[i]][odds,1+log2gen], type="l", col=mcol[log2gen+1], lty = mlty[i])
    }
}

legend("topright", legend = c(paste0("g = 2^", 0:21),"",paste("N =", c(16383,32767,65535,131071,262143,524287))), col = c(mcol, "white", rep("black",6)), lty = c(rep(1,23),mlty), ncol = 3)

dev.off()




#### ----------------------------------------------
## fit

fit <- vector("list", length(sigmaSquare))
mcol <- c("green4","cyan","blue","purple","red","black")

## fit the sigma2 with a function (actually 1-sigma2)
## we found that a function in A/x + B/(log2(x)*x^1.5 fitted best the data
## with generation = 2^x
## check the fit for all N
## and look at the best values for A and B
## lastly, plot the different sigma2 and the fit on the same plot
## FIGURE 3 of the preprint
pdf("fit.pdf")

for(i in 1:length(sigmaSquare)) {
    x <- 6:17
    y <- 1 - unlist(sigmaSquare[[i]][1,x+1])
    fit[[i]] <- nls(y ~ A/x + B/(log2(x)*x^1.5), start = list(A = 1.4, B = -2))

    x2 <- seq(from = 6, to = 17, by = 0.01)
    y2 <- predict(fit[[i]],list(x = x2))
    x3 <- seq(from = 0, to = 21, by = 0.01)
    y3 <- predict(fit[[i]],list(x = x3))

    x <- 0:21
    y <- 1 - unlist(sigmaSquare[[i]][1,x+1])

    if(i == 1) plot(x,y*x, ylim =c(0,2), col = mcol[i])
    points(x,y*x, ylim =c(0,2),col = mcol[i])
    lines(x2,y2*x2, col = mcol[i])
    lines(x3,y3*x3, col = mcol[i], lty = 2)
}

legend("topleft", legend = paste("N =", c(16383,32767,65535,131071,262143,524287)), col = mcol, lty = 1)

dev.off()

## extrapolation in N for a given g

N <- c(16383,32767,65535,131071,262143,524287)

data <- sapply(sigmaSquare, function(x) x[1,19+1])
epsilon <- 9

plot(N,round(data,epsilon), type = "b", ylim = c(min(data, na.rm = TRUE)-10^-epsilon, max(data, na.rm = TRUE)+10^-epsilon))
abline(h = max(round(data,epsilon)), col = "red")


               
#### ----------------------------------------------
## plot sigma2 vs g for the different N

N <- c(16383,32767,65535,131071,262143,524287)
mcol <- c("green4","cyan","blue","purple","red","black")

data <- sapply(sigmaSquare, function(x) unlist(x[1,]))

pdf("VarVsGeneration.pdf")
for(i in 1:length(data)) {
  x <- 0:(length(data[[i]])-1)
  if(i == 1) plot(x, data[[i]], col = mcol[i], type = "b", xlim = c(0,21), xlab = "generation (power of two)", ylab = "1-variance")
  points(x, data[[i]], col = mcol[i], type = "b")
}
legend("topleft", legend = paste("N =", c(16383,32767,65535,131071,262143,524287)), col = mcol, lty = 1)
dev.off()



#### ----------------------------------------------
#### ----------------------------------------------
## evolution of the mean mu (selection: alpha = 1.6)


data <- read.table("mu2_N16383_log2Gen21_alpha1.6.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data1 <- data
data <- read.table("mu2_N32767_log2Gen21_alpha1.6.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data2 <- data
data <- read.table("mu2_N65535_log2Gen21_alpha1.6.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data3 <- data
data <- read.table("mu2_N131071_log2Gen21_alpha1.6.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data4 <- data
data <- read.table("mu2_N262143_log2Gen21_alpha1.6.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data5 <- data
data <- read.table("mu2_N524287_log2Gen21_alpha1.6.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data6 <- data

mu <- list(data1,data2,data3,data4,data5,data6)
N <- c(16383,32767,65535,131071,262143,524287)
mcol <- c("green4","cyan","blue","purple","red","black")

for(i in 1:length(mu)) {
  if(i == 1) plot(mu[[i]]$generation, mu[[i]]$mu, col = mcol[1], type = "b")
  points(mu[[i]]$generation, mu[[i]]$mu, col = mcol[i], type = "b")
}

legend("topleft", legend = paste("N =", N), col = mcol, lty = 1)


#############################
## with different mu

data <- read.table("mu2_N16383_log2Gen21_alpha100.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data1 <- data
data <- read.table("mu2_N16383_log2Gen21_alpha0.1.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data2 <- data
data <- read.table("mu2_N16383_log2Gen21_alpha1.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data3 <- data
data <- read.table("mu2_N16383_log2Gen21_alpha10.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data4 <- data
data <- read.table("mu2_N16383_log2Gen21_alpha1.6.txt", fill = TRUE, header = FALSE, col.names = c("row", "generation", "mu"))
data <- data[seq(2, nrow(data), by = 2), -1]
data$generation <- as.numeric(data$generation)
data5 <- data

mu <- list(data1,data2,data3,data4,data5,data6)
alpha <- c(1.6, 5, 10, 15, 20, 25)
mcol <- c("green4","cyan","blue","purple","red","black")

for(i in 1:length(mu)) {
  if(i == 1) plot(mu[[i]]$generation, mu[[i]]$mu, col = mcol[1], type = "b")
  points(mu[[i]]$generation, mu[[i]]$mu, col = mcol[i], type = "b")
}

legend("topleft", legend = paste("alpha =", alpha), col = mcol, lty = 1)
