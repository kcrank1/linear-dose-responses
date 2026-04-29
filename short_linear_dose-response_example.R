#Tutorial for linearization
library(pracma) #for fderiv

#define function- make sure it's in log-log format

approx_bp <- function(x){
  x <- 10^x #we're giving the function a log-transformed dose
  a <- 0.104
  b <- 32.3
  log10(1 - (1 + (x/b))^-a) #the function returns a log-transformed response
}

#define a sequence of doses to use for numerical differentiation:lower than -13 hits float problems. 
xs <- seq(-13, 6, length.out=2000)

#run the function over the defined sequence
ys <- approx_bp(xs)

#take the first derivative using backwards finite differences
y1 <- fderiv(approx_bp, xs, n=1, method="backward" )
#take the second derivative using backwards finite differences
y2 <- fderiv(approx_bp, xs, n=2, method="backward")

#identify the infelction point - this is iterative, check visually to adjust the range to look for the inflection point
#identify the index value of the inflection point
inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
#find the x and y values for that index
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

# Plot
plot(xs, ys, type="l", col="gray", lwd=2, main="Approximate Beta-Poisson",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

##check here if the range for inflection point is correct, visually the point should match where the inflection point is. There can be inprecision, "wiggles", so narrow the range of xs if needed in line 27


#Find the linear fit
#pull out the the x and y values that are in the linear range
linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]

#use linear regression on the values
lin_model <- lm(y_linear ~ x_linear)



