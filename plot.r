## plots of the infinitesimal model


# Save these data
# saveRDS(sigmaSquare, file="sigmaSquare_N4000_log2Gen15.rds")
# to restore, do sigmaSquare <- readRDS("sigmaSquare_N4000_log2Gen15.rds")
sigmaSquare <- readRDS("sigmaSquare_N32767_log2Gen20.rds")

N <- dim(sigmaSquare)[1] -1
log2Generations <- dim(sigmaSquare)[2]-1

odds <- seq(2,N,by=2)  # odd sites because origin is at index = 1
mcol <- rainbow(log2Generations+1)

# First, the qualitative shape (no rescaling)
plot(odds,sigmaSquare[odds, 1] , xlab="site_index", ylab="variance", xlim=c(0,250), ylim=c(0,1), type="l", col=mcol[1], main = paste("N = ", N), lty = 2)
for(log2gen in 1:log2Generations) {
  points(odds, sigmaSquare[odds,1+log2gen], type="l", col=mcol[log2gen+1], lty = 2)
}


## fit
x <- 6:17
y <- 1 - unlist(sigmaSquare[1,x+1])
fit <- nls(y ~ A/x + B/(log2(x)*x^1.5), start = list(A = 1.4, B = -2))
x2 <- 0:20
y2 <- predict(fit,list(x = x2))

x <- 0:20
y <- 1 - unlist(sigmaSquare[1,x+1])

plot(x,y*x, ylim =c(0,1.4),col = "red")
lines(x2,y2*x2)

dev.off()
