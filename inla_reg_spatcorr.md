---
title: "Ground-Motion Model Regression with Spatial Correlations using INLA"
author:
  - Nicolas Kuehn^[Unversity of California, Los Angeles, kuehn@ucla.edu]
date: "04 April, 2023"
output:
  html_document:
  #pdf_document:
    keep_md: true
    toc: true
    toc_depth: 2
    number_sections: true
    highlight: tango
link-citations: yes
linkcolor: blue
citecolor: blue
urlcolor: blue
bibliography: /Users/nico/BIBLIOGRAPHY/BIBTEX/references.bib
---
  
  


# Introduction

This documet contains code to estimtate empirical ground-motion models (GMMs), while taking spatial correlations (of within-event/within-site residuals) into account.
The models are estimated with the R-INLA package (<https://www.r-inla.org/>) [@Rue2017].
Such regressions are similar to the ones proposed by @Jayaram2010a and @Ming2019, but in a Bayesian context.
The data is the Italian dat of the ITA18 model [@Lanzano2019], used also in @Caramenti2022.
The underlying functional form is the same as for the ITA18 model.

We fit a model of the form
$$
Y = f(\vec{c},\vec{x}) + \delta B + \delta S + \delta C(\vec{t}_s) + \delta W_0
$$
where $\delta B$ and $\delta S$ are the typical event and site terms, and $\delta C(\vec{t}_s)$ is a spatially correlated within-event/within-site residual term, dependent on the spatial coordinate $\vec{t}_s$ of the stations.
$f(\vec{c},\vec{x})$ is the base functional form, wth coefficients $\vec{c}$ and predictors $\vec{x}$.

Regressions are carried out using INLA @Rue2009 and the SPDE approximation @Lindgren2011.

# Getting Started

First, we load the packages required for the analysis.


```r
# load required packages
library(ggplot2)
library(INLA)
library(inlabru)
library(tidyverse)
library(cowplot)
```


```r
lw <- 1.5
sp <- 4
wid <- 8
asp <- 0.8

theme_set(theme_bw() + theme(#panel.grid.minor = element_blank(),
  axis.title = element_text(size = 30),
  axis.text = element_text(size = 20),
  plot.title = element_text(size = 30),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20),
  legend.key.width = unit(1, "cm"),
  legend.box.background = element_rect(colour = "black"),
  panel.grid = element_line(color = "gray",linewidth = 0.75)
))

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

cols <- c("darkblue", "dodgerblue1", "cadetblue2", "white")
cols2 <- rev(cols)
```

# Data

Load data, and define linear predictors and data frame for regression.


```r
### INLA
data <- read.csv(file.path('DATA', 'italian_data_pga_id_utm_stat.csv'))
data_dm <- rstan::read_rdump(file.path('DATA','dm_25x25.Rdata'))
data$R_epi <- sqrt(rowSums((data[,c('X_ev','Y_ev')] - data[,c('X_stat','Y_stat')])^2))

# Set linear predictors
mh = 5.5
mref = 5.324
h = 6.924
attach(data)
b1 = (mag-mh)*(mag<=mh)
b2 = (mag-mh)*(mag>mh)
c1 = (mag-mref)*log10(sqrt(JB_complete^2+h^2))
c2 = log10(sqrt(JB_complete^2+h^2))
c3 = sqrt(JB_complete^2+h^2)
f1 = as.numeric(fm_type_code == "SS")
f2 = as.numeric(fm_type_code == "TF")
k = log10(vs30/800)*(vs30<=1500)+log10(1500/800)*(vs30>1500)
y = log10(rotD50_pga)
detach(data)

n_rec <- length(b1)
eq <- data$EQID
stat <- data$STATID
n_eq <- max(eq)
n_stat <- max(stat)
n_cell <- data_dm$NCELL

data_reg <- data.frame(Y = y,
                       intercept = 1,
                       M1 = b1,
                       M2 = b2,
                       MlogR = c1,
                       logR = c2,
                       R = c3,
                       Fss = f1,
                       Frv = f2,
                       lnVS = k,
                       eq = eq,
                       stat = stat
)

# divide cell-specific distances by 100, to avoid small values of cell-specific
# coefficients and standard deviation
dm_sparse <- as(data_dm$RC/100,"dgCMatrix")
data_reg$idx_cell <- 1:nrow(data_reg)
```

# Model Estimation

Set priors for hyperparameters.
Defne both a Gamma prior fo the precision, as well as a penalized complexity prior.


```r
prior_prec_tau_lg    <- list(prec = list(prior = "loggamma", param = c(1.74, 0.0153)))
prior_prec_phiS2S_lg <- list(prec = list(prior = "loggamma", param = c(1.74, 0.0153)))
prior_prec_phiSS_lg  <- list(prec = list(prior = "loggamma", param = c(1.74, 0.0153)))

prior_prec_tau  <- list(prec = list(prior = 'pc.prec', param = c(0.4, 0.01)))
prior_prec_phiS2S  <- list(prec = list(prior = 'pc.prec', param = c(0.4, 0.01)))
prior_prec_phi0  <- list(prec = list(prior = 'pc.prec', param = c(0.4, 0.01)))
prior_prec_cell <- list(prec = list(prior = 'pc.prec', param = c(1, 0.01))) 
```

First, fit base models without spatial correlations (one with event and site terms,
one with only event terms).


```r
form <- Y ~ M1 + M2 + logR + MlogR + R + Fss + Frv + lnVS +
  f(eq, model = "iid", hyper = prior_prec_tau_lg) + 
  f(stat, model = "iid",hyper = prior_prec_phiS2S_lg)

fit_inla <- inla(form, 
                 data = data_reg,
                 family="gaussian",
                 control.family = list(hyper = prior_prec_phi0),
                 control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE)
)

form_eq <- Y ~ M1 + M2 + logR + MlogR + R + Fss + Frv + lnVS +
  f(eq, model = "iid", hyper = prior_prec_tau)

fit_inla_eq <- inla(form_eq, 
                    data = data_reg,
                    family="gaussian",
                    control.family = list(hyper = prior_prec_phi0),
                    control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE)
)

# add residals to data frame
data_reg$deltaW <- data_reg$Y - fit_inla_eq$summary.fitted.values$mean
data_reg$deltaWS <- data_reg$Y - fit_inla$summary.fitted.values$mean
data_reg$deltaW2 <- data_reg$Y - 
  (fit_inla$summary.fitted.values$mean - fit_inla$summary.random$stat$mean[data_reg$stat])
```

Next, estimate the model with spatial correlations.
First, we define the mesh, and plot the mesh togethe with the dta, colorcoded by event.


```r
### make mesh
co_eq <- data[,c("X_ev","Y_ev")]
co_stat <- data[,c("X_stat","Y_stat")]

# define relative coordinates
co_stat2 <- co_stat - co_eq
data_reg$X_stat <- co_stat2[,1]
data_reg$Y_stat <- co_stat2[,2]

max.edge2    <- 5 #0.04
bound.outer2 <- 20 #0.3
mesh = inla.mesh.2d(loc=as.matrix(co_stat2),
                    max.edge = c(1,5)*max.edge2,
                    # - use 5 times max.edge in the outer extension/offset/boundary
                    cutoff = max.edge2,
                    offset = c(5 * max.edge2, bound.outer2))
print(mesh$n) # number of mesh nodes
```

```
## [1] 4929
```

```r
print(length(mesh$graph$tv[,1])) # number of triangles
```

```
## [1] 9728
```

```r
ggplot() + theme_bw() + gg(mesh) +
  geom_point(data = data.frame(co_stat2), aes(x=X_stat,y=Y_stat, color = as.factor(eq))) +
  labs(x="X (km)", y="Y (km)") +
  theme(axis.title = element_text(size=30), axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        legend.position = 'none')
```

<img src="pictures_inla/make-mesh-1.png" width="50%" />



```r
# spde prior and define model
spde_rec <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = c(50, 0.5),
  # P(sigma > 1) = 0.01
  prior.sigma = c(.4, 0.01))

# group model with stations
A_rec_gr     <- inla.spde.make.A(mesh, loc = as.matrix(co_stat2),
                                 group = eq, n.group = n_eq)
idx_rec_gr   <- inla.spde.make.index("idx_rec_gr",spde_rec$n.spde,
                                     n.group = n_eq)

# build stack
stk_spatial_rec <- inla.stack(
  data = list(y = data_reg$Y),
  A = list(A_rec_gr, 1),
  effects = list(idx_rec_gr = idx_rec_gr,
                 data_reg))

form_rec <- y ~ 0 + intercept + M1 + M2 + logR + MlogR + R + Fss + Frv + lnVS +
  f(eq, model = "iid", hyper = prior_prec_tau_lg) + 
  f(stat, model = "iid",hyper = prior_prec_phiS2S_lg) +
  f(idx_rec_gr, model = spde_rec, group = idx_rec_gr.group,
                                  control.group = list(model = "iid"))

fit_inla_cor <- inla(form_rec,
                     data = inla.stack.data(stk_spatial_rec),
                     family="gaussian",
                     control.family = list(hyper = prior_prec_phi0),
                     control.predictor = list(A = inla.stack.A(stk_spatial_rec)),
                     control.compute = list(waic = TRUE),
                     control.inla = list(int.strategy = 'eb')
)
print(summary(fit_inla_cor))
```

```
## 
## Call:
##    c("inla.core(formula = formula, family = family, contrasts = contrasts, 
##    ", " data = data, quantiles = quantiles, E = E, offset = offset, ", " 
##    scale = scale, weights = weights, Ntrials = Ntrials, strata = strata, 
##    ", " lp.scale = lp.scale, link.covariates = link.covariates, verbose = 
##    verbose, ", " lincomb = lincomb, selection = selection, control.compute 
##    = control.compute, ", " control.predictor = control.predictor, 
##    control.family = control.family, ", " control.inla = control.inla, 
##    control.fixed = control.fixed, ", " control.mode = control.mode, 
##    control.expert = control.expert, ", " control.hazard = control.hazard, 
##    control.lincomb = control.lincomb, ", " control.update = 
##    control.update, control.lp.scale = control.lp.scale, ", " 
##    control.pardiso = control.pardiso, only.hyperparam = only.hyperparam, 
##    ", " inla.call = inla.call, inla.arg = inla.arg, num.threads = 
##    num.threads, ", " blas.num.threads = blas.num.threads, keep = keep, 
##    working.directory = working.directory, ", " silent = silent, inla.mode 
##    = inla.mode, safe = FALSE, debug = debug, ", " .parent.frame = 
##    .parent.frame)") 
## Time used:
##     Pre = 10.9, Running = 2870, Post = 21.6, Total = 2903 
## Fixed effects:
##             mean    sd 0.025quant 0.5quant 0.975quant   mode kld
## intercept  3.410 0.066      3.280    3.410      3.540  3.410   0
## M1         0.249 0.066      0.119    0.249      0.378  0.249   0
## M2         0.072 0.094     -0.113    0.072      0.256  0.072   0
## logR      -1.428 0.038     -1.502   -1.428     -1.355 -1.428   0
## MlogR      0.264 0.030      0.206    0.264      0.322  0.264   0
## R         -0.003 0.000     -0.003   -0.003     -0.002 -0.003   0
## Fss        0.152 0.041      0.071    0.152      0.233  0.152   0
## Frv        0.063 0.038     -0.011    0.063      0.137  0.063   0
## lnVS      -0.458 0.040     -0.537   -0.458     -0.378 -0.458   0
## 
## Random effects:
##   Name	  Model
##     eq IID model
##    stat IID model
##    idx_rec_gr SPDE2 model
## 
## Model hyperparameters:
##                                           mean     sd 0.025quant 0.5quant
## Precision for the Gaussian observations  73.67  2.908     68.101    73.62
## Precision for eq                         92.01 30.644     48.998    86.31
## Precision for stat                       25.02  1.559     22.087    24.97
## Range for idx_rec_gr                    131.38  9.608    113.668   130.95
## Stdev for idx_rec_gr                      0.22  0.008      0.204     0.22
##                                         0.975quant    mode
## Precision for the Gaussian observations     79.557  73.525
## Precision for eq                           167.832  75.696
## Precision for stat                          28.228  24.880
## Range for idx_rec_gr                       151.502 129.968
## Stdev for idx_rec_gr                         0.237   0.219
## 
## Watanabe-Akaike information criterion (WAIC) ...: -4872.50
## Effective number of parameters .................: 1635.80
## 
## Marginal log-Likelihood:  786.38 
##  is computed 
## Posterior summaries for the linear predictor and the fitted values are computed
## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```

## Cell-Specific Attenuation

Now also include cell-specific attenuation.


```r
form_rec_cell <- y ~ 0 + intercept + M1 + M2 + logR + MlogR + R + Fss + Frv + lnVS +
  f(eq, model = "iid", hyper = prior_prec_tau_lg) + 
  f(stat, model = "iid",hyper = prior_prec_phiS2S_lg) +
  f(idx_rec_gr, model = spde_rec, group = idx_rec_gr.group,
    control.group = list(model = "iid")) +
  f(idx_cell, model = "z", Z = dm_sparse, hyper = prior_prec_cell)

fit_inla_cell_cor <- inla(form_rec_cell,
                     data = inla.stack.data(stk_spatial_rec),
                     family="gaussian",
                     control.family = list(hyper = prior_prec_phi0),
                     control.predictor = list(A = inla.stack.A(stk_spatial_rec)),
                     control.compute = list(waic = TRUE),
                     control.inla = list(int.strategy = 'eb')
)
```

```
## Warning in inla.model.properties.generic(inla.trim.family(model), mm[names(mm) == : Model 'z' in section 'latent' is marked as 'experimental'; changes may appear at any time.
##   Use this model with extra care!!! Further warnings are disabled.
```

```r
form_cell <- y ~ 0 + intercept + M1 + M2 + logR + MlogR + R + Fss + Frv + lnVS +
  f(eq, model = "iid", hyper = prior_prec_tau_lg) + 
  f(stat, model = "iid",hyper = prior_prec_phiS2S_lg) +
  f(idx_cell, model = "z", Z = dm_sparse, hyper = prior_prec_cell)

fit_inla_cell <- inla(form_cell,
                      data = data_reg,
                      family="gaussian",
                      control.family = list(hyper = prior_prec_phi0),
                      control.compute = list(waic = TRUE)
)
```

## Non-Stationary Model

Now a non-stationary model, based on ideas from @Kuehn2020 and @Liu2022a.
Here, the spatial range $\ell$ is dependent on distance to the source
$$
\ell = \exp \left[ a + b \frac{\min(R_{EPI}, 80)}{80}\right]
$$
For implementation details, see Chapter 5 of @Krainski2019 (<https://becarioprecario.bitbucket.io/spde-gitbook/index.html>).


```r
# set mesh distances
ref_d <- 80
D_b_mesh <- sqrt(mesh$loc[,1]^2 + mesh$loc[,2]^2)
D_b_mesh[D_b_mesh > ref_d] <- ref_d
D_b_mesh <- D_b_mesh/ref_d

# correlation depends on distance to basin
nu <- 1 
alpha <- nu + 2 / 2
# log(kappa)
logkappa0 <- log(8 * nu) / 2
# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0

# SPDE model
spde_rec_ns <- inla.spde2.matern(mesh, 
                                 B.tau = cbind(logtau0, -1, nu, nu * D_b_mesh), 
                                 B.kappa = cbind(logkappa0, 0, -1, -1 * D_b_mesh),
                                 theta.prior.mean = rep(0, 3), 
                                 #theta.prior.mean = c(-0.7, 2.3, 0.5),
                                 theta.prior.prec = rep(1, 3))

# build stack
stk_spatial_rec <- inla.stack(
  data = list(y = data_reg$Y),
  A = list(A_rec_gr, 1),
  effects = list(idx_rec_gr = idx_rec_gr,
                 data_reg))

form_rec_ns <- y ~ 0 + intercept + M1 + M2 + logR + MlogR + R + Fss + Frv + lnVS +
  f(eq, model = "iid", hyper = prior_prec_tau_lg) + 
  f(stat, model = "iid",hyper = prior_prec_phiS2S_lg) +
  f(idx_rec_gr, model = spde_rec_ns, group = idx_rec_gr.group,
    control.group = list(model = "iid"))

fit_inla_cor_ns <- inla(form_rec_ns,
                        data = inla.stack.data(stk_spatial_rec),
                        family="gaussian",
                        control.family = list(hyper = prior_prec_phi0),
                        control.predictor = list(A = inla.stack.A(stk_spatial_rec)),
                        control.compute = list(waic = TRUE)
)
print(summary(fit_inla_cor_ns))
```

```
## 
## Call:
##    c("inla.core(formula = formula, family = family, contrasts = contrasts, 
##    ", " data = data, quantiles = quantiles, E = E, offset = offset, ", " 
##    scale = scale, weights = weights, Ntrials = Ntrials, strata = strata, 
##    ", " lp.scale = lp.scale, link.covariates = link.covariates, verbose = 
##    verbose, ", " lincomb = lincomb, selection = selection, control.compute 
##    = control.compute, ", " control.predictor = control.predictor, 
##    control.family = control.family, ", " control.inla = control.inla, 
##    control.fixed = control.fixed, ", " control.mode = control.mode, 
##    control.expert = control.expert, ", " control.hazard = control.hazard, 
##    control.lincomb = control.lincomb, ", " control.update = 
##    control.update, control.lp.scale = control.lp.scale, ", " 
##    control.pardiso = control.pardiso, only.hyperparam = only.hyperparam, 
##    ", " inla.call = inla.call, inla.arg = inla.arg, num.threads = 
##    num.threads, ", " blas.num.threads = blas.num.threads, keep = keep, 
##    working.directory = working.directory, ", " silent = silent, inla.mode 
##    = inla.mode, safe = FALSE, debug = debug, ", " .parent.frame = 
##    .parent.frame)") 
## Time used:
##     Pre = 9.52, Running = 48942, Post = 73.2, Total = 49025 
## Fixed effects:
##             mean    sd 0.025quant 0.5quant 0.975quant   mode kld
## intercept  3.427 0.077      3.275    3.427      3.577  3.427   0
## M1         0.266 0.071      0.126    0.266      0.404  0.266   0
## M2         0.096 0.096     -0.093    0.096      0.285  0.096   0
## logR      -1.446 0.046     -1.536   -1.446     -1.356 -1.446   0
## MlogR      0.245 0.034      0.180    0.245      0.311  0.245   0
## R         -0.003 0.000     -0.003   -0.003     -0.002 -0.003   0
## Fss        0.133 0.042      0.052    0.133      0.215  0.133   0
## Frv        0.071 0.039     -0.005    0.071      0.147  0.071   0
## lnVS      -0.463 0.040     -0.540   -0.463     -0.385 -0.463   0
## 
## Random effects:
##   Name	  Model
##     eq IID model
##    stat IID model
##    idx_rec_gr SPDE2 model
## 
## Model hyperparameters:
##                                           mean     sd 0.025quant 0.5quant
## Precision for the Gaussian observations  92.48  3.779      85.15    92.45
## Precision for eq                        109.43 38.544      56.92   101.79
## Precision for stat                       25.74  1.564      22.80    25.68
## Theta1 for idx_rec_gr                    -1.56  0.027      -1.61    -1.56
## Theta2 for idx_rec_gr                     3.22  0.092       3.05     3.22
## Theta3 for idx_rec_gr                     2.17  0.112       1.94     2.17
##                                         0.975quant  mode
## Precision for the Gaussian observations     100.03 92.45
## Precision for eq                            205.87 87.90
## Precision for stat                           28.96 25.58
## Theta1 for idx_rec_gr                        -1.50 -1.56
## Theta2 for idx_rec_gr                         3.41  3.21
## Theta3 for idx_rec_gr                         2.38  2.18
## 
## Watanabe-Akaike information criterion (WAIC) ...: -5729.54
## Effective number of parameters .................: 1779.54
## 
## Marginal log-Likelihood:  951.38 
##  is computed 
## Posterior summaries for the linear predictor and the fitted values are computed
## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```

# Model Comparison

## Posterior Plots

First, plot the posterior of the spatial range of the models which include spatial correlation.
For thenon-stationary model, we plot the value of the spatial range corresponding o $R_{EPI} = 0$km.


```r
xlab <- 'spatial range (km)'
ylab <- 'posterior density'
rbind(data.frame(fit_inla_cor$marginals.hyperpar$`Range for idx_rec_gr`, mod = "Stationary"),
      data.frame(fit_inla_cell_cor$marginals.hyperpar$`Range for idx_rec_gr`, mod = "Cell"),
      data.frame(inla.tmarginal(function(x) exp(x), 
                                      fit_inla_cor_ns$internal.marginals.hyperpar$`Theta2 for idx_rec_gr`),
                 mod = "Non-Stationary")) %>%
  ggplot() +
  geom_line(aes(x = x, y = y, color = mod), linewidth = lw) +
  labs(x = xlab, y = ylab) +
  guides(color = guide_legend(title=NULL)) +
  scale_y_continuous(expand = c(0, 0))
```

<div class="figure">
<img src="pictures_inla/plot-range-1.png" alt="Posterior distribution of spatial range for different models." width="50%" />
<p class="caption">Posterior distribution of spatial range for different models.</p>
</div>

Now, plot the posterior distribution of the standard deviation $\phi_c$, which is associated with the spatial correlation structure.


```r
xlab <- expression(atop(paste(phi[c])))
ylab <- 'posterior density'
rbind(data.frame(fit_inla_cor$marginals.hyperpar$`Stdev for idx_rec_gr`, mod = "Stationary"),
      data.frame(fit_inla_cell_cor$marginals.hyperpar$`Stdev for idx_rec_gr`, mod = "Cell"),
      data.frame(inla.tmarginal(function(x) exp(x), 
                                      fit_inla_cor_ns$internal.marginals.hyperpar$`Theta1 for idx_rec_gr`),
                 mod = "Non-Stationary")) %>%
  ggplot() +
  geom_line(aes(x = x, y = y, color = mod), linewidth = lw) +
  labs(x = xlab, y = ylab) +
  guides(color = guide_legend(title=NULL)) +
  scale_y_continuous(expand = c(0, 0))
```

<div class="figure">
<img src="pictures_inla/plot-phi-c-1.png" alt="Posterior distribution of phi_c for different models." width="50%" />
<p class="caption">Posterior distribution of phi_c for different models.</p>
</div>

## Plot of Correlations

Now we plot the correlation kernel for the stationary and non-stationary correlation functions.
The non-stationary correlation depends on the distance to the source, so is different for sies with different coordinates.

We use the median of the posterior distribution as the parameters of the spatial correlation structure.
First, we calculate the precision matrix using `inla.spde2.precision`, which is then used to calculate the correlation values at the mesh nodes.
The function `book.spatial.correlation` is from @Krainsk2019 (<https://www.r-inla.org/learnmore/books>).


```r
theta_ns <- fit_inla_cor_ns$summary.hyperpar[c(4,5,6),"0.5quant"]
Q_ns <- inla.spde2.precision(spde_rec_ns, theta = theta_ns)

theta_stat <- fit_inla_cor$internal.summary.hyperpar[c(4,5),"0.5quant"]
Q_stat <- inla.spde2.precision(spde_rec, theta = theta_stat)

corr_ns1 <- book.spatial.correlation(Q_ns, c(10, 0), mesh)
corr_ns2 <- book.spatial.correlation(Q_ns, c(100, 0), mesh)
corr_stat <- book.spatial.correlation(Q_stat, c(100, 0), mesh)

df <- data.frame(label = "Non-stationary - Site 1", x = -80, y= 240)
pl_ns1 <- ggplot()+gg(mesh, color=corr_ns1, nx = 300, ny =300) +
  scale_fill_gradientn(colours = cols2, name = "", limits = c(0.,1)) +
  geom_point(data.frame(x = 0, y = 0), mapping = aes(x= x, y = y), shape = 8, size = 2) +
  geom_point(data.frame(x = 10, y = 0), mapping = aes(x= x, y = y), shape = 17, color = 'red', size = 2) +
  theme(panel.grid = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(limits = c(-260,260),expand = c(0,0)) +
  scale_y_continuous(limits = c(-260,260),expand = c(0,0)) +
  geom_label(df, mapping = aes(x = x, y = y, label = label), size = 30/2.8)

df <- data.frame(label = "Non-stationary - Site 2", x = -80, y= 240)
pl_ns2 <- ggplot()+gg(mesh, color=corr_ns2, nx = 300, ny =300) +
  scale_fill_gradientn(colours = cols2, name = "", limits = c(0.,1)) +
  geom_point(data.frame(x = 0, y = 0), mapping = aes(x= x, y = y), shape = 8, size = 2) +
  geom_point(data.frame(x = 100, y = 0), mapping = aes(x= x, y = y), shape = 17, color = 'red', size = 2) +
  theme(panel.grid = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(limits = c(-260,260),expand = c(0,0)) +
  scale_y_continuous(limits = c(-260,260),expand = c(0,0)) +
  geom_label(df, mapping = aes(x = x, y = y, label = label), size = 30/2.8)

df <- data.frame(label = "Stationary", x = -100, y= 240)
pl_stat <- ggplot()+gg(mesh, color=corr_stat, nx = 500, ny =500) +
  scale_fill_gradientn(colours = cols2, name = "", limits = c(0.,1)) +
  geom_point(data.frame(x = 0, y = 0), mapping = aes(x= x, y = y), shape = 8, size = 2) +
  geom_point(data.frame(x = 100, y = 0), mapping = aes(x= x, y = y), shape = 17, color = 'red', size = 2) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.15,0.2)) +
  scale_x_continuous(limits = c(-260,260),expand = c(0,0)) +
  scale_y_continuous(limits = c(-260,260),expand = c(0,0)) +
  geom_label(df, mapping = aes(x = x, y = y, label = label), size = 30/2.8)

plot_grid(pl_stat, pl_ns1, pl_ns2, ncol = 2, nrow = 2)
```

<img src="pictures_inla/plot-correlations-1.png" width="100%" />


# References