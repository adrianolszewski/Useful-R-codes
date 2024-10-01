This code allows one to pool over a set of emmeans objects.

**Scenario**

Assume you work with missing  data imputed via mice.
Then, you want to apply a statistical model from a package which is partially supported by the emmeans package.
"Partially" means that the emmean object can be succesfully created, but not all features of the package are easily accessible from the emmeans interface.

Let me show you an example where I use GEE estimation method to fit a simple linear general linear model for repeated observations.
Let's say it's a GEE-fitted MMRM.

/ Side note: MMRM stands for Mixed Model for Repeated Measurements. Despite its name it's a fixed-effect only model accounting for within-subject correlations and potential heteroscedasticity.
Normally it's fitted using the Generalized Least Squares (GLS) estimation, but since the normality of normalized residuals is compromised, I want to employ GEE, which does not need this assumption for valid asymptotic inference.
So instead of using the `nlme::gls()` or `mmrm::mmrm()` routines, I employ the best (to me) and most flexible package for GEE estimation: `glmtoolbox::glmgee`
/

To create an emmean object, I want the robust, bias-corrected estimator of empirical covariance (clustered sandwich with enhancements).
I can provide the `vcov = vcov(m, type = "bias-corrected")` parameter, but how to manage the fact that the analyses are pooled?

Let me show you an example.

-----------------

First, the data:
```{r}
> d <- structure(list(ID = structure(c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 
3L, 4L, 4L, 4L, 5L, 5L, 5L, 6L, 6L, 6L, 7L, 7L, 7L, 8L, 8L, 8L, 
9L, 9L, 9L, 10L, 10L, 10L, 11L, 11L, 11L, 12L, 12L, 12L, 13L, 
13L, 13L, 14L, 14L, 14L), levels = c("1", "2", "3", "4", "5", 
"6", "7", "8", "9", "10", "11", "12", "13", "14"), class = "factor"), 
Arm = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L), levels = c("A", "B"), class = "factor"), Visit = structure(c(1L, 
2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 
2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 
2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L), levels = c("V1", 
"V2", "V3"), class = c("ordered", "factor")), Visit_non_ord = structure(c(1L, 
2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 
2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 
2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L), levels = c("V1", 
"V2", "V3"), class = "factor"), Painscore = c(1L, 2L, 4L, 
5L, 6L, NA, 4L, 2L, 3L, 4L, 3L, 4L, 5L, NA, 5L, 6L, 7L, 6L, 
4L, 3L, 5L, 1L, NA, 4L, 5L, 4L, 5L, 4L, 6L, NA, 0L, 4L, 3L, 
4L, NA, 4L, 3L, 4L, NA, 6L, 7L, 8L)), row.names = c(NA, -42L
), class = "data.frame")
```

which looks like that:
```
  ID Arm Visit Visit_non_ord Painscore
1   1   A    V1            V1         1
2   1   A    V2            V2         2
3   1   A    V3            V3         4
4   2   A    V1            V1         5
5   2   A    V2            V2         6
6   2   A    V3            V3        NA
7   3   A    V1            V1         4
8   3   A    V2            V2         2
9   3   A    V3            V3         3
10  4   A    V1            V1         4
11  4   A    V2            V2         3
12  4   A    V3            V3         4
13  5   A    V1            V1         5
14  5   A    V2            V2        NA
... etc
```

It's data from a virtual longitudinal (repeated-data) study comparing some outcome called PainScore across two treatment arms.
Each patient has three visits at which some PainScore is assessed.
Within-subject (over-time) trends of means are assessed.
Let's skip the discussion about the nture of this variable and the justification of validity of this analysis - it's just a working example.

# 1) Imputation
Don't take this much seriously, it's just an illustration.

```{r}
> library(mice)

> imp <- mice(d, m=5)
> imp$predictorMatrix[1:4, 1:5] <- 0
> imp$predictorMatrix[5, 1] <- 0
> imp <- mice(d, m=5, predictorMatrix = imp$predictorMatrix)
> imp_list <- complete(imp, "all")
```

# 2) Fitting GEE model and creating the emmeans objects:

**Notice the** `vcov =  vcov(m, type = "bias-corrected"))` **line**.
The body of the lapply's function can be very simple or very complex - whatever you need to do the anallysis. What matters is that it produces a list of valid emmeans objects.

```
> trend_emms <- lapply(imp_list, 
                    function(dat) {
                        m <- glmgee(Painscore ~ Visit * Arm, 
                                    family = gaussian(link = "identity"), 
                                    id = ID,  
                                    data=dat, 
                                    corstr = "unstructured")
                        
                        emmeans(m, 
                                specs = ~ Visit | Arm,
                                adjust="none", 
                                vcov =  vcov(m, type = "bias-corrected"))   # this is not accessible from within emmeans
                      })
```

resulting in:
```
> trend_emms
$`1`
Arm = A:
 Visit emmean    SE df lower.CL upper.CL
 V1      4.14 0.642 36     2.84     5.45
 V2      3.57 0.845 36     1.86     5.29
 V3      4.57 0.398 36     3.76     5.38

Arm = B:
 Visit emmean    SE df lower.CL upper.CL
 V1      3.29 0.873 36     1.52     5.06
 V2      4.14 0.797 36     2.53     5.76
 V3      4.86 0.642 36     3.55     6.16

Covariance estimate used: user-supplied 
Confidence level used: 0.95 

$`2`
Arm = A:
 Visit emmean    SE df lower.CL upper.CL
 V1      4.14 0.642 36     2.84     5.45
 V2      4.29 0.934 36     2.39     6.18
 V3      4.86 0.549 36     3.74     5.97

Arm = B:
 Visit emmean    SE df lower.CL upper.CL
 V1      3.29 0.873 36     1.52     5.06
 V2      5.43 0.661 36     4.09     6.77
 V3      4.71 0.735 36     3.22     6.20

Covariance estimate used: user-supplied 
Confidence level used: 0.95 
...
```

OK, we have the list of emmeans objects, now we want to pool them according to Rubin's rules and using the Barnard's small-sample method for determining the degrees of freedom, implemented in `emmeans:::emm_basis.mira()`.

# 3) The pool_emmeans()
We will use the code from this function and adjust it slightly:
```{r}
pool_emmeans <- function(emmeans_list) {
  bas = emmeans_list[[1]]
  k = length(emmeans_list)
  V = 1/k * bas@V
  allb = cbind(bas@bhat, matrix(0, nrow = length(bas@bhat), ncol = k - 1))
  
  for (i in 1 + seq_len(k - 1)) {
    basi = emmeans_list[[i]]
    V = V + 1/k * basi@V
    allb[, i] = basi@bhat
  }
  bas@bhat = apply(allb, 1, mean)
  notna = which(!is.na(bas@bhat))
  bas@dfargs$m = k
  bas@dfargs$df1 = bas@dffun
  bas@dfargs$B = cov(t(allb[notna, , drop = FALSE]))
  bas@V = bas@dfargs$T = V + (k + 1)/k * bas@dfargs$B
  bas@dffun = function(a, dfargs) {
    dfcom = dfargs$df1(a, dfargs)
    with(dfargs, {
      b = sum(a * (B %*% a))
      t = sum(a * (T %*% a))
      lambda = (1 + 1/m) * b/t
      dfold = (m - 1)/lambda^2
      dfobs = (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
      ifelse(is.infinite(dfcom), dfold, dfold * dfobs/(dfold + dfobs))
    })
  }
   return(bas)
}
```

# 4) Let's pool the emmeans objects and calculate some contrasts (Wald's inference)

```{r}
> pooled_ems <- pool_emmeans(trend_emms)
> contrast(pooled_ems, "poly", max.degree=2, by = "Arm")
Arm = A:
 contrast  estimate    SE   df t.ratio p.value
 linear       0.457 0.583 26.9   0.784  0.4399
 quadratic    1.029 1.482 15.7   0.694  0.4979

Arm = B:
 contrast  estimate    SE   df t.ratio p.value
 linear       1.514 0.680 19.7   2.227  0.0378
 quadratic   -1.571 1.748 10.8  -0.899  0.3883
```

As you can see - it works and the degrees of freedom (df) are pooled properly as well.

----

Just to be sure, let's recalculate it with different estimator of covariance:
```{r}
> trend_emms <- lapply(imp_list, 
                    function(dat) {
                        m <- glmgee(Painscore ~ Visit * Arm, 
                                    family = gaussian(link = "identity"), 
                                    id = ID,  
                                    data=dat, 
                                    corstr = "unstructured")
                        
                        emmeans(m, 
                                specs = ~ Visit | Arm,
                                adjust="none", 
                                vcov =  vcov(m, type = "model"))   # different estimator
                      })

> pooled_ems <- pool_emmeans(trend_emms)
> contrast(pooled_ems, "poly", max.degree=2, by = "Arm")
Arm = A:
 contrast  estimate    SE    df t.ratio p.value
 linear       0.457 0.389 16.88   1.174  0.2565
 quadratic    1.029 1.090  7.09   0.944  0.3762

Arm = B:
 contrast  estimate    SE    df t.ratio p.value
 linear       1.514 0.464  8.53   3.265  0.0105
 quadratic   -1.571 1.340  4.51  -1.173  0.2991
```

OK, the `vcov` parameter was handled properly.

--------

**Last, but not least...**

I would like to sincerely thank the creators of both packages: 
Professor **Russell V. Lenth** (_Department of Statistics and Actuarial Science, The University of Iowa Iowa City_), author of the [emmeans](https://github.com/rvlenth/emmeans),
and Professor **Luis Hernando Vanegas** (_Department of Statistics, The National University of Colombia_) [glmtoolbox](https://cran.r-project.org/web/packages/glmtoolbox/index.html) 
for their hard work, for sharing it with the community, and for their active support in resolving challenging issues. A wonderful tandem of tools has been created!
