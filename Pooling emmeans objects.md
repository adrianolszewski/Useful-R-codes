This code allows one to pool over a set of emmeans objects.

**Scenario**

Assume you work with missing  data imputed via mice.
Then, you want to apply a statistical model from a package which is partially supported by the emmeans package.
"Partially" means that the emmean object can be succesfully created, but not all features of the package are easily accessible from the emmeans interface if operating directly on imputed datasets.

Let me show you an example where I use GEE estimation method to fit a simple linear general linear model for repeated observations.
Let's say it's a GEE-fitted MMRM.

/ Side note: MMRM stands for Mixed Model for Repeated Measurements. Despite its name it's a fixed-effect only model accounting for within-subject correlations and potential heteroscedasticity.
Normally it's fitted using the Generalized Least Squares (GLS) estimation, but since the normality of normalized residuals is compromised, I want to employ GEE, which does not need this assumption for valid asymptotic inference.
So instead of using the `nlme::gls()` or `mmrm::mmrm()` routines, I employ the best (to me) and most flexible package for GEE estimation: `glmtoolbox::glmgee`
/

To create an emmean object, I want the robust, bias-corrected estimator of empirical covariance (clustered sandwich with enhancements).
In emmeans I can provide the `vcov = vcov(m, type = "bias-corrected")` parameter, but how to manage the fact that the analyses are pooled?

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
    2L, 2L), levels = c("A", "B"), class = "factor"), PainScore = c(1L, 
    2L, 4L, 5L, 6L, NA, 4L, 2L, 3L, 4L, 3L, 4L, 5L, NA, 5L, 6L, 
    7L, 6L, 4L, 3L, 5L, 1L, NA, 4L, 5L, 4L, 5L, 4L, 6L, NA, 0L, 
    4L, 3L, 4L, NA, 4L, 3L, 4L, NA, 6L, 7L, 8L), Visit = structure(c(1L, 
    2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 
    2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 
    2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L), levels = c("V1", 
    "V2", "V3"), class = "factor")), row.names = c(NA, -42L), class = "data.frame")
```

which looks like that:
```
> head(d, 15)
   ID Arm PainScore Visit
1   1   A         1    V1
2   1   A         2    V2
3   1   A         4    V3
4   2   A         5    V1
5   2   A         6    V2
6   2   A        NA    V3
7   3   A         4    V1
8   3   A         2    V2
9   3   A         3    V3
10  4   A         4    V1
11  4   A         3    V2
12  4   A         4    V3
13  5   A         5    V1
14  5   A        NA    V2
15  5   A         5    V3
... etc
```

It's fake data from a virtual longitudinal (repeated-data) study, which aim is to assess the impact of some treatment on the patient-reported outcome called PainScore. 
Each patient has three visits at which the PainScore is assessed. **Let's perform per-visit comparisons of mean PainScore between the treatment arms**.

/ _(let's assume, just for the sake of simplicity, that PainScore is a numerical endpoint, meaningfully sumamrized with arithmetic means)._ /


# 1) Imputation
Don't take this much seriously, it's just an illustration.

```{r}
> library(mice)

> imp <- mice(d, m=5)
> imp$predictorMatrix[1:3, 1:4] <- 0
> imp$predictorMatrix[4, 1] <- 0
> imp <- mice(d, m=5, predictorMatrix = imp$predictorMatrix)

> imp_list <- complete(imp, "all") # we will need this later
```

# 2) Fitting GEE model and creating the emmeans objects

Let's start with the simplest approach to better illustrate what we want to achieve.

Let's fit the models on the imputed datasets.
```{r}
> my_models_per_arm <- with(imp, glmgee(PainScore ~ Visit * Arm, 
                                      family = gaussian(link = "identity"), 
                                      id = ID,  
                                      corstr = "exchangeable"))
# Did it work?
> my_models_per_arm$analyses[1:2]
[[1]]

        Variance function:  gaussian
                     Link:  identity
    Correlation structure:  Exchangeable 

[[2]]

        Variance function:  gaussian
                     Link:  identity
    Correlation structure:  Exchangeable

# Let's have a glance...
> summary(my_models_per_arm$analyses[[1]])

Sample size
   Number of observations:  42
       Number of clusters:  14 
             Cluster size:  3 
*************************************************************
Model
        Variance function:  gaussian
            Link function:  identity
    Correlation structure:  Exchangeable
*************************************************************
Coefficients
             Estimate Std.Error  z-value   Pr(>|z|)
(Intercept)   4.14286   0.55064  7.52368 5.3255e-14
VisitV2      -0.42857   0.48894 -0.87652   0.380746
VisitV3       0.14286   0.55064  0.25944   0.795298
ArmB         -0.85714   0.92896 -0.92269   0.356170
VisitV2:ArmB  1.85714   0.77497  2.39640   0.016557
VisitV3:ArmB  1.00000   0.77873  1.28415   0.199090
                                                   
Dispersion    2.78571                              
*************************************************************
Working correlation
    [1]  [2]  [3] 
[1] 1.00 0.54 0.54
[2] 0.54 1.00 0.54
```

Good! Now, let's create pooled emmeans over all those models and test some contrast using the Wald's approach...
```
> (my_emmeans_per_arm <- emmeans(my_models_per_arm, 
                                specs = ~ Arm * Visit,
                                adjust="none"))

 Arm Visit emmean    SE   df lower.CL upper.CL
 A   V1      4.14 0.551  NaN     3.06     5.22
 B   V1      3.29 0.748 34.1     1.77     4.81
 A   V2      3.89 0.717 31.4     2.42     5.35
 B   V2      4.83 0.501 27.3     3.80     5.86
 A   V3      4.37 0.503 14.6     3.30     5.45
 B   V3      4.71 0.742 18.6     3.16     6.27

Covariance estimate used: robust  # <---- notice this information
Confidence level used: 0.95 

> set.seed(1000)  # for the MVT adjustment using Monte Carlo integration
> (update(contrast(my_emmeans_per_arm,
                 list(
                     "V1: B vs. A" = c(-1,1,   0,0,   0,0), # direction of comparisons is up to you
                     "V2: B vs. A" = c( 0,0,  -1,1,   0,0),
                     "V3: B vs. A" = c( 0,0,   0,0,  -1,1)
                 )),
        adjust="mvt", level = 0.95, infer = c(TRUE, TRUE)))

 contrast    estimate    SE   df lower.CL upper.CL t.ratio p.value
 V1: B vs. A   -0.857 0.929 34.1    -3.12     1.41  -0.923  0.6654
 V2: B vs. A    0.943 0.852 33.0    -1.14     3.03   1.106  0.5403
 V3: B vs. A    0.343 0.803 29.4    -1.64     2.32   0.427  0.9455

Confidence level used: 0.95 
Conf-level adjustment: mvt method for 3 estimates 
P value adjustment: mvt method for 3 tests 
```

OK! We were able to fit models pear each individual imputed dataset and create a pooled emmeans object.
**But now we want a different covariance estimator, namely the bias-corrected one.**

Seems a pretty simple task! 
The `emmeans()` allows us to specify the `vcov` parameter.
The `glmtoolbox` package offers a way to calculate different estimators by calling the `vcov()` function, for example: `vcov(model, type = "bias-corrected")`.
So what stops us from combining the two?

Well, there is no just a single specific model to refer to but a rather **a list of them**.
And the _emmeans_ object **already contains pooled estimates**.

Let's try another approach. We will:
1) iterate over the list of the imputed datasets,
2) in each iteration - fit a model to each dataset,
3) create a corresponding emmeans object with the necessary covariance estimator
4) pool the estimated model coefficients and covariance matrices stored in every emmeans object. We will refer to it as "_pooling the emmeans objects_".

Let's do it!
The body of the lapply's function can be very simple or very complex - whatever you need to do the anallysis. What matters is that it produces a list of valid emmeans objects.

First let's recreate the previous results choosing the "_robust_" estimator. Just to validate the new approach.

```
> my_emmeans_per_arm <- lapply(imp_list, 
                    function(dat) {
                        m <- glmgee(PainScore ~ Visit * Arm, 
                                    family = gaussian(link = "identity"), 
                                    id = ID,  
                                    data=dat, 
                                    corstr = "exchangeable")
                        
                        emmeans(m, 
                                specs = ~ Arm * Visit,
                                adjust="none", 
                                vcov =  vcov(m, type = "robust"))   # now the model "m" is accessible!
                     })
```

...resulting in:
```
> my_emmeans_per_arm[1:2]
$`1`
 Arm Visit emmean    SE df lower.CL upper.CL
 A   V1      4.14 0.551 36     3.03     5.26
 B   V1      3.29 0.748 36     1.77     4.80
 A   V2      3.71 0.691 36     2.31     5.12
 B   V2      4.71 0.439 36     3.82     5.60
 A   V3      4.29 0.389 36     3.50     5.08
 B   V3      4.43 0.601 36     3.21     5.65

Covariance estimate used: user-supplied    # <--- A-ha! emmeans noticed our request!
Confidence level used: 0.95 

$`2`
 Arm Visit emmean    SE df lower.CL upper.CL
 A   V1      4.14 0.551 36     3.03     5.26
 B   V1      3.29 0.748 36     1.77     4.80
 A   V2      3.86 0.683 36     2.47     5.24
 B   V2      4.71 0.439 36     3.82     5.60
 A   V3      4.00 0.571 36     2.84     5.16
 B   V3      4.57 0.829 36     2.89     6.25

Covariance estimate used: user-supplied 
Confidence level used: 0.95
```

Now, having the list of emmeans objects, we want to pool them according to Rubin's rules and the Barnard's small-sample method for pooling the degrees of freedom. We will slightly adjust the code implemented in the `emmeans:::emm_basis.mira()` function:

# 3) The pool_emmeans() function
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

# 4) Let's pool the emmeans objects and calculate the same contrasts as previously
```{r}
> (pooled_ems <- pool_emmeans(my_emmeans_per_arm))

 Arm Visit emmean    SE   df lower.CL upper.CL
 A   V1      4.14 0.551  NaN     3.06     5.22
 B   V1      3.29 0.748 34.1     1.77     4.81
 A   V2      3.89 0.717 31.4     2.42     5.35
 B   V2      4.83 0.501 27.3     3.80     5.86
 A   V3      4.37 0.503 14.6     3.30     5.45
 B   V3      4.71 0.742 18.6     3.16     6.27

Covariance estimate used: user-supplied 
Confidence level used: 0.95 

> set.seed(1000)  # for the MVT adjustment using Monte Carlo integration
> (update(contrast(pooled_ems,
                 list(
                     "V1: B vs. A" = c(-1,1,   0,0,   0,0), # direction of comparisons is up to you
                     "V2: B vs. A" = c( 0,0,  -1,1,   0,0),
                     "V3: B vs. A" = c( 0,0,   0,0,  -1,1)
                 )),
        adjust="mvt", level = 0.95, infer = c(TRUE, TRUE)))

 contrast    estimate    SE   df lower.CL upper.CL t.ratio p.value
 V1: B vs. A   -0.857 0.929 34.1    -3.12     1.41  -0.923  0.6654
 V2: B vs. A    0.943 0.852 33.0    -1.14     3.03   1.106  0.5403
 V3: B vs. A    0.343 0.803 29.4    -1.64     2.32   0.427  0.9455

Confidence level used: 0.95 
Conf-level adjustment: mvt method for 3 estimates 
P value adjustment: mvt method for 3 tests 
```

The results perfectly agree with the previous approach.
**Now, knowing that we can control things this way, let's finally request the bias-corrected estimtor of covariance!**

```{r}
 my_emmeans_per_arm <- lapply(imp_list, 
                               function(dat) {
                                   m <- glmgee(PainScore ~ Visit * Arm, 
                                               family = gaussian(link = "identity"), 
                                               id = ID,  
                                               data=dat, 
                                               corstr = "exchangeable")
                                   
                                   emmeans(m, 
                                           specs = ~ Arm * Visit,
                                           adjust="none", 
                                           vcov =  vcov(m, type = "bias-corrected")) # Finally!
                               })

> (pooled_ems <- pool_emmeans(my_emmeans_per_arm))
 Arm Visit emmean    SE   df lower.CL upper.CL
 A   V1      4.14 0.642  NaN     2.88     5.40   #  <--- Notice larger standard errors (SE)!
 B   V1      3.29 0.873 34.1     1.51     5.06
 A   V2      3.89 0.830 32.2     2.19     5.58
 B   V2      4.83 0.575 29.3     3.65     6.00
 A   V3      4.37 0.562 17.9     3.19     5.55
 B   V3      4.71 0.837 22.0     2.98     6.45

Covariance estimate used: user-supplied 
Confidence level used: 0.95

> set.seed(1000)  # for the MVT adjustment using Monte Carlo integration
> (update(contrast(pooled_ems,
                 list(
                     "V1: B vs. A" = c(-1,1,   0,0,   0,0), # direction of comparisons is up to you
                     "V2: B vs. A" = c( 0,0,  -1,1,   0,0),
                     "V3: B vs. A" = c( 0,0,   0,0,  -1,1)
                 )),
        adjust="mvt", level = 0.95, infer = c(TRUE, TRUE)))

 contrast    estimate    SE   df lower.CL upper.CL t.ratio p.value
 V1: B vs. A   -0.857 1.084 34.1    -3.50     1.79  -0.791  0.7535
 V2: B vs. A    0.943 0.991 33.4    -1.48     3.37   0.952  0.6444
 V3: B vs. A    0.343 0.926 30.9    -1.93     2.62   0.370  0.9628

Confidence level used: 0.95 
Conf-level adjustment: mvt method for 3 estimates 
P value adjustment: mvt method for 3 tests 
```

**As we could see - the `vcov` parameter was handled properly.**

Remember this trick with lapply + pool_emmeans(), because this way you will be able to do more than using defaults.

--------

**Last, but not least...**

I would like to sincerely thank the creators of both packages: 
Professor **Russell V. Lenth** (_Department of Statistics and Actuarial Science, The University of Iowa Iowa City_), author of the [emmeans](https://github.com/rvlenth/emmeans),
and Professor **Luis Hernando Vanegas** (_Department of Statistics, The National University of Colombia_) [glmtoolbox](https://cran.r-project.org/web/packages/glmtoolbox/index.html) 
for their hard work, for sharing it with the community, and for their active support in resolving challenging issues. A wonderful tandem of tools has been created!
