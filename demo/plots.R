library(makehams) 

# survival function
plot( function(t) {
        tpx(t,x=20,s=0)
      }, 0, 10, type='l', ylab="tp[x]", xlab="Time", 
      main="Survival function", col='blue')

# force of mortality
plot(function(t) uxt(t,x=20), 0, 5, type='l', col='blue', 
     ylab="u[x]+t", xlab="Time", main="Force of mortality")
abline(v=2, col='red')
legend("bottomright", leg=c("force of mortality", "select period"),
       lty=c(1,1), col=c('blue', 'red'))

# discount factor
plot(function(t) v(0.05,n=t), 0, 20, ylab="Discount factor", xlab="Time",
     main="Comparison of discount factors", col='blue', ylim=c(0.2,1.0))
plot(function(t) tEx(t, x=60, s = 0), 0, 20, add=TRUE, col='red')
legend("bottomleft", leg=c("Present value factor", "Actuarial discount factor"),
       lty=c(1,1), col=c('blue','red'))

