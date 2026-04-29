Dose-Response Linear Ranges
================
2026-04-29

For full code and figures, download the .rmd file! For a shorter non-specific tutorial, see the short_linear_dose-response_example.R file.

# Introduction

This is the supporting code for *Simplified Linear Dose-Response Models
for Low-Dose QMRAs.* 

Authors: Katherine Crank, Emily Clements, Daniel
Gerrity

Here, we analyze **dose-response functions** for common enteric
pathogens and define the linear range of the function and the equation
of the line. We also report the maximum doses that can be utilized with the
linear models.

For each function, we:

1.  Transform the dose-response curve to log-log space (log<sub>10</sub>
    of dose and probability of infection) .

2.  Compute numerical derivatives using **finite differences**: the
    first derivative provides the slope at each point, and the second
    derivative identifies inflection points.

3.  Identify the inflection point, defined as the dose at which the
    second derivative is minimized. This point corresponds to the
    transition where the function begins to have linear behavior.

4.  Fit a linear regression to the portion of the curve to the left of
    the inflection point, capturing the linear region.

5.  Overlay the linear fit and its equation directly on the curve for
    visual reference.

This approach allows us to quantify the low-dose, linear region of each
dose-response function and compare slope and intercept parameters across
different models.

``` r
library(pracma) #numerical differentiation
library(kableExtra) #nice tables
library(gsl)#hypergeometric function
library(ggplot2)#graphing
library(dplyr)#data manipulation
library(tidyr)#data manipulation
library(DT)#data manipulation
library(Cairo)#graphics rendering device - not required

# Initialize results table
results <- data.frame(
  results=character(),
  Function=character(),
  Max_dose=numeric(),
 
  Linear_Intercept=numeric(),
  Linear_Slope=numeric(),
  Name=character(),
  Pathogen=character(),
  stringsAsFactors=FALSE
)
fit_list <- list()
## sequence for numerical second differentiation. lower than -13 hits float problems. 
xs <- seq(-13, 6, length.out=2000)
```

# Norovirus Dose-Response Functions

``` r
# -------------------------------
# Dose-response model formulas, parameters, and metadata
# -------------------------------
nov_info <- data.frame(
  Name = c(
    "Approximate Beta-Poisson",
    "Hypergeometric (Teunis 2008)",
    "Hypergeometric (Teunis 2020)",
    "Fractional Poisson (Aggregated)",
    "Fractional Poisson (No Aggregation)"
  ),
  Function = c("approx_bp", "hyper", "new_hyper", "fp_agg", "fp"),
  Formula = c(
    "y = 1 - (1 + x/b)^-a",
    "y = 1 - 1F1(alpha, alpha+beta, -x)",
    "y = 1 - 1F1(alpha_n, alpha_n+beta_n, -x)",
    "y = p * (1 - exp(-x/alpha1))",
    "y = p * (1 - exp(-x/alpha2))"
  ),
  Parameters = c(
    "a=0.104, b=32.3",
    "alpha=0.04, beta=0.055",
    "alpha_n=0.393, beta_n=0.767",
    "p=0.72, alpha1=1106",
    "p=0.72, alpha2=1"
  ),
  Citation = c(
    "Van Abel et al., 2017",
    "Teunis et al., 2008",
    "Teunis et al., 2020",
    "Messner et al., 2014",
    "Messner et al., 2014"
  ),
  stringsAsFactors = FALSE
)

# Display 
kable(nov_info, caption="Norovirus dose-response models") %>%
  kable_styling(full_width=F, position="center")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Norovirus dose-response models
</caption>
<thead>
<tr>
<th style="text-align:left;">
Name
</th>
<th style="text-align:left;">
Function
</th>
<th style="text-align:left;">
Formula
</th>
<th style="text-align:left;">
Parameters
</th>
<th style="text-align:left;">
Citation
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Approximate Beta-Poisson
</td>
<td style="text-align:left;">
approx_bp
</td>
<td style="text-align:left;">
y = 1 - (1 + x/b)^-a
</td>
<td style="text-align:left;">
a=0.104, b=32.3
</td>
<td style="text-align:left;">
Van Abel et al., 2017
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypergeometric (Teunis 2008)
</td>
<td style="text-align:left;">
hyper
</td>
<td style="text-align:left;">
y = 1 - 1F1(alpha, alpha+beta, -x)
</td>
<td style="text-align:left;">
alpha=0.04, beta=0.055
</td>
<td style="text-align:left;">
Teunis et al., 2008
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypergeometric (Teunis 2020)
</td>
<td style="text-align:left;">
new_hyper
</td>
<td style="text-align:left;">
y = 1 - 1F1(alpha_n, alpha_n+beta_n, -x)
</td>
<td style="text-align:left;">
alpha_n=0.393, beta_n=0.767
</td>
<td style="text-align:left;">
Teunis et al., 2020
</td>
</tr>
<tr>
<td style="text-align:left;">
Fractional Poisson (Aggregated)
</td>
<td style="text-align:left;">
fp_agg
</td>
<td style="text-align:left;">
y = p \* (1 - exp(-x/alpha1))
</td>
<td style="text-align:left;">
p=0.72, alpha1=1106
</td>
<td style="text-align:left;">
Messner et al., 2014
</td>
</tr>
<tr>
<td style="text-align:left;">
Fractional Poisson (No Aggregation)
</td>
<td style="text-align:left;">
fp
</td>
<td style="text-align:left;">
y = p \* (1 - exp(-x/alpha2))
</td>
<td style="text-align:left;">
p=0.72, alpha2=1
</td>
<td style="text-align:left;">
Messner et al., 2014
</td>
</tr>
</tbody>
</table>

## 1. Approximate Beta-Poisson

Van Abel et al 2017, alpha = 0.104, beta=32.3

``` r
# -------------------------------
# 1. Approximate Beta-Poisson (Van Abel 2017)
# -------------------------------
approx_bp <- function(x){
  x <- 10^x
  a <- 0.104
  b <- 32.3
  log10(1 - (1 + (x/b))^-a)
}


ys <- approx_bp(xs)
y1 <- fderiv(approx_bp, xs, n=1, method="backward" )
y2 <- fderiv(approx_bp, xs, n=2, method="backward")

# Inflection
inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

# Plot
plot(xs, ys, type="l", col="gray", lwd=2, main="Approximate Beta-Poisson",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

# Linear fit
linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

# Add formula
slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)


# Legend
legend("topleft", legend=c("Function","1st derivative","2nd derivative",formula_text),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Save results
results <- rbind(results, data.frame(
  Function="approx_bp",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
     Name= "Approximate Beta-Poisson",
  Pathogen = "Norovirus"
))
```

## 2. Hypergeometric

Teunis et al., 2008; alpha=0.04, beta=0.055, no aggregation

``` r
# -------------------------------
# 2. Hypergeometric (Teunis 2008)
# -------------------------------
hyper <- function(x){
  x <- 10^x
  alpha <- 0.04
  beta <- 0.055
  log10(1 - hyperg_1F1(alpha, (alpha+beta), -x))
}


ys <- hyper(xs)
y1 <- fderiv(hyper, xs, n=1, method="backward")
y2 <- fderiv(hyper, xs, n=2, method="backward")

inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

plot(xs, ys, type="l", col="gray", lwd=2, main="Hypergeometric (Teunis 2008)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)

legend("topleft", legend=c("Function","1st derivative","2nd derivative",formula_text),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="hyper",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
  Name="Hypergeometric (Teunis 2008)",
  Pathogen = "Norovirus"
))
```

## 3. New Hypergeometric (Teunis 2020)

Teunis et al., 2020; alpha=0.04, beta=0.05

``` r
# -------------------------------
# 3. New Hypergeometric (Teunis 2020)
# -------------------------------
new_hyper <- function(x){
  x <- 10^x
  alpha_n <- 0.393
  beta_n <- 0.767
  log10(1 - hyperg_1F1(alpha_n, (alpha_n + beta_n), -x))
}


ys <- new_hyper(xs)
y1 <- fderiv(new_hyper, xs, n=1, method="backward")
y2 <- fderiv(new_hyper, xs, n=2, method="backward")

inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

plot(xs, ys, type="l", col="gray", lwd=2, main="Hypergeometric (Teunis 2020)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)

legend("topleft", legend=c("Function","1st derivative","2nd derivative",formula_text),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="new_hyper",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
   Name="Hypergeometric (Teunis 2020)",
  Pathogen = "Norovirus"
))
```

## 4. Fractional Poisson (Aggregated)

Messner et al., 2014; Atmar et al., 2014; Atmar et al., 2008; P=0.72,
alpha=1106

``` r
# -------------------------------
# 4. Fractional Poisson Aggregated
# -------------------------------
fp_agg <- function(x){
  x <- 10^x
  p <- 0.72
  alpha1 <- 1106
  log10(p * (1 - exp(-x / alpha1)))
}


ys <- fp_agg(xs)
y1 <- fderiv(fp_agg, xs, n=1, method="backward")
y2 <- fderiv(fp_agg, xs, n=2, method="backward")

inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

plot(xs, ys, type="l", col="gray", lwd=2, main="Fractional Poisson (Aggregated)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)

legend("topleft", legend=c("Function","1st derivative","2nd derivative",formula_text),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="fp_agg",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
  Name= "Fractional Poisson (Aggregated)",
  Pathogen = "Norovirus"
))
```

## 5. Fractional Poisson (No Aggregation)

Messner et al., 2014; Atmar et al., 2014; Atmar et al., 2008; P=0.72,
alpha=1

``` r
# -------------------------------
# 5. Fractional Poisson (No Aggregation)
# -------------------------------
fp <- function(x){
  x <- 10^x
  p1 <- 0.72
  alpha2 <- 1
  log10(p1 * (1 - exp(-x / alpha2)))
}


ys <- fp(xs)
y1 <- fderiv(fp, xs, n=1, method="backward")
y2 <- fderiv(fp, xs, n=2, method="backward")

inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

plot(xs, ys, type="l", col="gray", lwd=2, main="Fractional Poisson (No Aggregation)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)

legend("topleft", legend=c("Function","1st derivative","2nd derivative","formula_text"),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="fp",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
  Name="Fractional Poisson (No Aggregation)",
  Pathogen = "Norovirus"
))
```

# Enterovirus Dose-Response Functions

``` r
# -------------------------------
# Dose-response model formulas, parameters, and metadata
# -------------------------------

env_info <- data.frame(
  Name = c(
    " Beta-Poisson"
  ),
  Function = c("ent_bp"),
  Formula = c(
    "y = 1 - ((1 + (dose / beta))^(-1 * alpha))"),
  Parameters = c(
  
    "alpha=0.253, beta=0.426"

  ),
  Citation = c(
    "Ward et al., 1986"
  ),
  stringsAsFactors = FALSE
)

# Display nicely in R Markdown
kable(env_info, caption="Enterovirus dose-response models") %>%
  kable_styling(full_width=F, position="center")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Enterovirus dose-response models
</caption>
<thead>
<tr>
<th style="text-align:left;">
Name
</th>
<th style="text-align:left;">
Function
</th>
<th style="text-align:left;">
Formula
</th>
<th style="text-align:left;">
Parameters
</th>
<th style="text-align:left;">
Citation
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Beta-Poisson
</td>
<td style="text-align:left;">
ent_bp
</td>
<td style="text-align:left;">
y = 1 - ((1 + (dose / beta))^(-1 \* alpha))
</td>
<td style="text-align:left;">
alpha=0.253, beta=0.426
</td>
<td style="text-align:left;">
Ward et al., 1986
</td>
</tr>
</tbody>
</table>

## 1. Beta-Poisson

``` r
# -------------------------------
# 1.  Beta-Poisson (Ward 1986)
# -------------------------------
ent_bp <- function(x){
  x <- 10^x
  a <- 0.253
  b <- 0.426
  log10(1 - ((1 + (x / b))^(-1 * a)))
}


ys <- ent_bp(xs)
y1 <- fderiv(ent_bp, xs, n=1, method="backward")
y2 <- fderiv(ent_bp, xs, n=2, method="backward")

# Inflection
inflect_idx <- which(xs > -1)[which.min(y2[xs > -1])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

# Plot
plot(xs, ys, type="l", col="gray", lwd=2, main="Beta-Poisson",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

# Linear fit
linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

# Add formula
slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)


# Legend
legend("topleft", legend=c("Function","1st derivative","2nd derivative",formula_text),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# Save results
results <- rbind(results, data.frame(
  Function="ent_bp",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
  Name="Beta-Poisson",
  Pathogen = "Enterovirus"
))
```

# Adenovirus Dose-Response Functions

``` r
# -------------------------------
# Adenovirus dose-response model formulas, parameters, and metadata
# -------------------------------

adeno_info <- data.frame(
  Name = c(
    "Exponential (Crabtree et al., 1997)",
    "Exact Beta-Poisson (Teunis et al., 2016)"
  ),
  Function = c("adeno_exp_crabtree", "adeno_bp_exact"),
  Formula = c(
    "y = 1 - exp(-1 * r * dose)",
    "y = 1 - hyperg_1F1(alpha, alpha+beta, -1 * dose)"
  ),
  Parameters = c(
    "r = 0.4172",
    "alpha = 5.11, beta = 2.8"
  ),
  Citation = c(
    "Crabtree et al., 1997",
    "Teunis et al., 2016"
  ),
  stringsAsFactors = FALSE
)

kable(adeno_info, caption="Adenovirus dose-response models") %>%
  kable_styling(full_width=F, position="center")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Adenovirus dose-response models
</caption>
<thead>
<tr>
<th style="text-align:left;">
Name
</th>
<th style="text-align:left;">
Function
</th>
<th style="text-align:left;">
Formula
</th>
<th style="text-align:left;">
Parameters
</th>
<th style="text-align:left;">
Citation
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Exponential (Crabtree et al., 1997)
</td>
<td style="text-align:left;">
adeno_exp_crabtree
</td>
<td style="text-align:left;">
y = 1 - exp(-1 \* r \* dose)
</td>
<td style="text-align:left;">
r = 0.4172
</td>
<td style="text-align:left;">
Crabtree et al., 1997
</td>
</tr>
<tr>
<td style="text-align:left;">
Exact Beta-Poisson (Teunis et al., 2016)
</td>
<td style="text-align:left;">
adeno_bp_exact
</td>
<td style="text-align:left;">
y = 1 - hyperg_1F1(alpha, alpha+beta, -1 \* dose)
</td>
<td style="text-align:left;">
alpha = 5.11, beta = 2.8
</td>
<td style="text-align:left;">
Teunis et al., 2016
</td>
</tr>
</tbody>
</table>

inp## 1. Exponential (Crabtree et al., 1997)

``` r
adeno_exp_crabtree <- function(x){
  x <- 10^x; r <- 0.4172
  log10(1 - exp(-1 * r * x))
}

ys <- adeno_exp_crabtree(xs)
y1 <- fderiv(adeno_exp_crabtree, xs, n=1, method="backward")
y2 <- fderiv(adeno_exp_crabtree, xs, n=2, method="backward")
#inflect_idx <- which.min(y2)
inflect_idx <- which(xs > -1)[which.min(y2[xs > -1])] #this narrows the range where the minimum index value is found, to those where the value of x (not the index) are greater than in this case -5. Can be specified using visual examination of the curve
inflect_x <- xs[inflect_idx]; inflect_y <- ys[inflect_idx]
plot(xs, ys, type="l", col="gray", lwd=2, main="Exponential (Crabtree et al. 1997)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2); lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
linear_idx <- which(xs <= inflect_x)
lin_model <- lm(ys[linear_idx] ~ xs[linear_idx])
lines(xs[linear_idx], predict(lin_model), col="blue", lwd=2)
slope <- coef(lin_model)[2]; intercept <- coef(lin_model)[1]
legend("topleft", legend=c("Function","1st derivative","2nd derivative",
       bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="adeno_exp_crabtree", Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3), Linear_Slope=round(slope,3),
  Name="Exponential (Crabtree et al. 1997)", Pathogen="Adenovirus"))
```

## 2. Exact Beta-Poisson (Teunis et al., 2016)

``` r
adeno_bp_exact <- function(x){
  x <- 10^x
  alpha <- 5.11
  beta <- 2.8
  log10(1 - hyperg_1F1(alpha, (alpha + beta), -1 * x))
}


ys <- adeno_bp_exact(xs)
y1 <- fderiv(adeno_bp_exact, xs, n=1, method="backward")
y2 <- fderiv(adeno_bp_exact, xs, n=2, method="backward")

# Inflection
inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

# Plot
plot(xs, ys, type="l", col="gray", lwd=2, main="Exact Beta-Poisson (Teunis et al., 2016)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

# Linear fit
linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

# Add formula
slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)

# Legend
legend("topleft", legend=c("Function","1st derivative","2nd derivative",formula_text),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# Save results
results <- rbind(results, data.frame(
  Function="adeno_bp_exact",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
  Name="Exact Beta-Poisson (Teunis et al., 2016)",
  Pathogen="Adenovirus"
))
```

# Giardia Dose-Response Functions

``` r
# -------------------------------
# Dose-response model formulas, parameters, and metadata
# -------------------------------

gia_info <- data.frame(
  Name = c(
    "Exponential"
  ),
  Function = c("gia_exp"),
  Formula = c(
    "y = 1-exp(-1 * r * dose)"),
  Parameters = c(
  
    "r = 0.0199"

  ),
  Citation = c(
    "Teunis et al., 1997"
  ),
  stringsAsFactors = FALSE
)

# Display nicely in R Markdown
kable(gia_info, caption="Giardia dose-response models") %>%
  kable_styling(full_width=F, position="left")
```

<table class="table" style="width: auto !important; ">
<caption>
Giardia dose-response models
</caption>
<thead>
<tr>
<th style="text-align:left;">
Name
</th>
<th style="text-align:left;">
Function
</th>
<th style="text-align:left;">
Formula
</th>
<th style="text-align:left;">
Parameters
</th>
<th style="text-align:left;">
Citation
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Exponential
</td>
<td style="text-align:left;">
gia_exp
</td>
<td style="text-align:left;">
y = 1-exp(-1 \* r \* dose)
</td>
<td style="text-align:left;">
r = 0.0199
</td>
<td style="text-align:left;">
Teunis et al., 1997
</td>
</tr>
</tbody>
</table>

## 1. Exponential

``` r
# -------------------------------
# 1.  Exponential (Teunis et al. 1997)
# -------------------------------
gia_exp <- function(x){
  x <- 10^x
  r <- 0.0199
  log10(1-exp(-1 * r * x))
}


ys <- gia_exp(xs)
y1 <- fderiv(gia_exp, xs, n=1, method="backward")
y2 <- fderiv(gia_exp, xs, n=2, method="backward")

# Inflection
inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]
inflect_y <- ys[inflect_idx]

# Plot
plot(xs, ys, type="l", col="gray", lwd=2, main="Exponential",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2)
lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
text(inflect_x, inflect_y, labels=paste0("Inflection\n(", round(inflect_x,3), ",", round(inflect_y,3),")"), pos=4, col="darkred")

# Linear fit
linear_idx <- which(xs <= inflect_x)
x_linear <- xs[linear_idx]
y_linear <- ys[linear_idx]
lin_model <- lm(y_linear ~ x_linear)
y_fit <- predict(lin_model)
lines(x_linear, y_fit, col="blue", lwd=2)

# Add formula
slope <- coef(lin_model)[2]
intercept <- coef(lin_model)[1]
formula_text <- bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)


# Legend
legend("topleft", legend=c("Function","1st derivative","2nd derivative",formula_text),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# Save results
results <- rbind(results, data.frame(
  Function="gia_exp",
  Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3),
  Linear_Slope=round(slope,3),
  Name="Exponential",
  Pathogen = "Giardia"
))
```

# Cryptosporidium Dose-Response Functions

``` r
# -------------------------------
# Cryptosporidium dose-response model formulas, parameters, and metadata
# -------------------------------

crypto_info <- data.frame(
  Name = c(
    "Exponential (EPA 2006)",
    "Exponential (Haas et al., 1999; Barbeau et al., 2000; Zhang et al., 2012)",
    "Fractional Poisson (Messner & Berger, 2016)",
    "Beta-Poisson (Messner & Berger, 2016)",
    "Exponential with Immunity (Messner & Berger, 2016)"
  ),
  Function = c(
    "crypto_exp_epa",
    "crypto_exp_haas",
    "crypto_fracpois",
    "crypto_bp",
    "crypto_exp_immunity"
  ),
  Formula = c(
    "y = 1 - exp(-1 * r * dose)",
    "y = 1 - exp(-1 * r * dose)",
    "y = p * (1 - exp((-1 * dose) / alpha))",
    "y = 1 - (1 + (dose / beta))^(-alpha)",
    "y = p * (1 - exp(-1 * r * dose))"
  ),
  Parameters = c(
    "r = 0.09",
    "r = 0.00419",
    "p = 0.737, alpha = 1",
    "alpha = 0.116, beta = 0.121",
    "p = 0.737, r = 0.608"
  ),
  Citation = c(
    "EPA, 2006",
    "Haas et al., 1999; Barbeau et al., 2000; Zhang et al., 2012",
    "Messner & Berger, 2016",
    "Messner & Berger, 2016",
    "Messner & Berger, 2016"
  ),
  stringsAsFactors = FALSE
)

kable(crypto_info, caption="Cryptosporidium dose-response models") %>%
  kable_styling(full_width=F, position="center")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Cryptosporidium dose-response models
</caption>
<thead>
<tr>
<th style="text-align:left;">
Name
</th>
<th style="text-align:left;">
Function
</th>
<th style="text-align:left;">
Formula
</th>
<th style="text-align:left;">
Parameters
</th>
<th style="text-align:left;">
Citation
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Exponential (EPA 2006)
</td>
<td style="text-align:left;">
crypto_exp_epa
</td>
<td style="text-align:left;">
y = 1 - exp(-1 \* r \* dose)
</td>
<td style="text-align:left;">
r = 0.09
</td>
<td style="text-align:left;">
EPA, 2006
</td>
</tr>
<tr>
<td style="text-align:left;">
Exponential (Haas et al., 1999; Barbeau et al., 2000; Zhang et al.,
2012)
</td>
<td style="text-align:left;">
crypto_exp_haas
</td>
<td style="text-align:left;">
y = 1 - exp(-1 \* r \* dose)
</td>
<td style="text-align:left;">
r = 0.00419
</td>
<td style="text-align:left;">
Haas et al., 1999; Barbeau et al., 2000; Zhang et al., 2012
</td>
</tr>
<tr>
<td style="text-align:left;">
Fractional Poisson (Messner & Berger, 2016)
</td>
<td style="text-align:left;">
crypto_fracpois
</td>
<td style="text-align:left;">
y = p \* (1 - exp((-1 \* dose) / alpha))
</td>
<td style="text-align:left;">
p = 0.737, alpha = 1
</td>
<td style="text-align:left;">
Messner & Berger, 2016
</td>
</tr>
<tr>
<td style="text-align:left;">
Beta-Poisson (Messner & Berger, 2016)
</td>
<td style="text-align:left;">
crypto_bp
</td>
<td style="text-align:left;">
y = 1 - (1 + (dose / beta))^(-alpha)
</td>
<td style="text-align:left;">
alpha = 0.116, beta = 0.121
</td>
<td style="text-align:left;">
Messner & Berger, 2016
</td>
</tr>
<tr>
<td style="text-align:left;">
Exponential with Immunity (Messner & Berger, 2016)
</td>
<td style="text-align:left;">
crypto_exp_immunity
</td>
<td style="text-align:left;">
y = p \* (1 - exp(-1 \* r \* dose))
</td>
<td style="text-align:left;">
p = 0.737, r = 0.608
</td>
<td style="text-align:left;">
Messner & Berger, 2016
</td>
</tr>
</tbody>
</table>

## 1. Exponential (EPA 2006)

``` r
crypto_exp_epa <- function(x){
  x <- 10^x
  r <- 0.09
  log10(1 - exp(-1 * r * x))
}

ys <- crypto_exp_epa(xs)
y1 <- fderiv(crypto_exp_epa, xs, n=1, method="backward")
y2 <- fderiv(crypto_exp_epa, xs, n=2, method="backward")
inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]; inflect_y <- ys[inflect_idx]
plot(xs, ys, type="l", col="gray", lwd=2, main="Exponential (EPA 2006)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2); lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
linear_idx <- which(xs <= inflect_x)
lin_model <- lm(ys[linear_idx] ~ xs[linear_idx])
lines(xs[linear_idx], predict(lin_model), col="blue", lwd=2)
slope <- coef(lin_model)[2]; intercept <- coef(lin_model)[1]
legend("topleft", legend=c("Function","1st derivative","2nd derivative",
       bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="crypto_exp_epa", Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3), Linear_Slope=round(slope,3),
  Name="Exponential (EPA 2006)", Pathogen="Cryptosporidium"))
```

## 2. Exponential (Haas et al., 1999; Barbeau et al., 2000; Zhang et al., 2012)

``` r
crypto_exp_haas <- function(x){
  x <- 10^x; r <- 0.00419
  log10(1 - exp(-1 * r * x))
}

ys <- crypto_exp_haas(xs)
y1 <- fderiv(crypto_exp_haas, xs, n=1, method="backward")
y2 <- fderiv(crypto_exp_haas, xs, n=2, method="backward")
inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]; inflect_y <- ys[inflect_idx]
plot(xs, ys, type="l", col="gray", lwd=2, main="Exponential (Haas et al. 1999 etc.)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2); lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
linear_idx <- which(xs <= inflect_x)
lin_model <- lm(ys[linear_idx] ~ xs[linear_idx])
lines(xs[linear_idx], predict(lin_model), col="blue", lwd=2)
slope <- coef(lin_model)[2]; intercept <- coef(lin_model)[1]
legend("topleft", legend=c("Function","1st derivative","2nd derivative",
       bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="crypto_exp_haas", Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3), Linear_Slope=round(slope,3),
  Name="Exponential (Haas/Barbeau/Zhang)", Pathogen="Cryptosporidium"))
```

## 3. Fractional Poisson (Messner & Berger, 2016)

``` r
crypto_fracpois <- function(x){
  x <- 10^x; p <- 0.737; alpha <- 1
  log10(p * (1 - exp((-1 * x) / alpha)))
}

ys <- crypto_fracpois(xs)
y1 <- fderiv(crypto_fracpois, xs, n=1, method="backward")
y2 <- fderiv(crypto_fracpois, xs, n=2, method="backward")
inflect_idx <- which(xs > -1)[which.min(y2[xs > -1])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]; inflect_y <- ys[inflect_idx]
plot(xs, ys, type="l", col="gray", lwd=2, main="Fractional Poisson (M&B 2016)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2); lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
linear_idx <- which(xs <= inflect_x)
lin_model <- lm(ys[linear_idx] ~ xs[linear_idx])
lines(xs[linear_idx], predict(lin_model), col="blue", lwd=2)
slope <- coef(lin_model)[2]; intercept <- coef(lin_model)[1]
legend("topleft", legend=c("Function","1st derivative","2nd derivative",
       bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="crypto_fracpois", Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3), Linear_Slope=round(slope,3),
  Name="Fractional Poisson (M&B 2016)", Pathogen="Cryptosporidium"))
```

## 4. Beta-Poisson (Messner & Berger, 2016)

``` r
crypto_bp <- function(x){
  x <- 10^x; alpha <- 0.116; beta <- 0.121
  log10(1 - (1 + (x / beta))^(-alpha))
}

ys <- crypto_bp(xs)
y1 <- fderiv(crypto_bp, xs, n=1, method="backward")
y2 <- fderiv(crypto_bp, xs, n=2, method="backward")
inflect_idx <- which(xs > -1)[which.min(y2[xs > -1])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]; inflect_y <- ys[inflect_idx]
plot(xs, ys, type="l", col="gray", lwd=2, main="Beta-Poisson (M&B 2016)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2); lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
linear_idx <- which(xs <= inflect_x)
lin_model <- lm(ys[linear_idx] ~ xs[linear_idx])
lines(xs[linear_idx], predict(lin_model), col="blue", lwd=2)
slope <- coef(lin_model)[2]; intercept <- coef(lin_model)[1]
legend("topleft", legend=c("Function","1st derivative","2nd derivative",
       bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="crypto_bp", Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3), Linear_Slope=round(slope,3),
  Name="Beta-Poisson (M&B 2016)", Pathogen="Cryptosporidium"))
```

## 5. Exponential with Immunity (Messner & Berger, 2016)

``` r
crypto_exp_immunity <- function(x){
  x <- 10^x; p <- 0.737; r <- 0.608
  log10(p * (1 - exp(-1 * r * x)))
}
ys <- crypto_exp_immunity(xs)
y1 <- fderiv(crypto_exp_immunity, xs, n=1, method="backward")
y2 <- fderiv(crypto_exp_immunity, xs, n=2, method="backward")
inflect_idx <- which(xs > 0)[which.min(y2[xs > 0])] #visually identified range where inflection point is, bounded the search for the minimum in that area to avoid wiggles. which.min looks for rank not value.
inflect_x <- xs[inflect_idx]; inflect_y <- ys[inflect_idx]
plot(xs, ys, type="l", col="gray", lwd=2, main="Exponential w/ Immunity (M&B 2016)",
     xlab="log(x-dose)", ylab="log(y-prob)")
lines(xs, y1, col="black", lty=2); lines(xs, y2, col="red", lty=3)
points(inflect_x, inflect_y, col="darkred", pch=19, cex=1.2)
linear_idx <- which(xs <= inflect_x)
lin_model <- lm(ys[linear_idx] ~ xs[linear_idx])
lines(xs[linear_idx], predict(lin_model), col="blue", lwd=2)
slope <- coef(lin_model)[2]; intercept <- coef(lin_model)[1]
legend("topleft", legend=c("Function","1st derivative","2nd derivative",
       bquote(y == .(round(intercept,3)) + .(round(slope,3))*x)),
       col=c("gray","black","red","blue"), lty=c(1,2,3,1), lwd=c(2,1,1,2), bty="n")
```

![](Dose-Response-Linear-Approximation_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
results <- rbind(results, data.frame(
  Function="crypto_exp_immunity", Max_dose=round(inflect_x,3),
  Linear_Intercept=round(intercept,3), Linear_Slope=round(slope,3),
  Name="Exponential with Immunity (M&B 2016)", Pathogen="Cryptosporidium"))
```

