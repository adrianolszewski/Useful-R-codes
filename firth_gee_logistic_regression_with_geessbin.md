# Firth-like penalised GEE-fit logistic regression (for repeated observations e.g. in longitudinal studies) in R

This code allows one to test contrasts in longitudinal studies with binary endpoint under small trial size and small prevalence (low fraction of the recorded events).
Normally you would fit a logistic regression using **GEE** (_Generalized Estimating Equations_) to obtain population-average estimates (like in clinical trials), 
or **GLMM** (_Generalized Linear Mixed Model_) to obtain conditional, subject-specific estimates. **But this is likely to fail under very small number of events.**

Therefore we will focus on the GEE estimation combined with the **Firth-like penalisation**.

Used packages: 
- [geessbin](https://github.com/rtishii/geessbin): Modified Generalized Estimating Equations for Binary Outcome. Analyze small-sample clustered or longitudinal data with binary outcome using modified generalized estimating equations (GEE) with bias-adjusted covariance estimator.
- [emmeans](https://rvlenth.github.io/emmeans/index.html): Estimated Marginal Means, aka Least-Squares Means. Obtain estimated marginal means (EMMs) for many linear, generalized linear, and mixed models. Compute contrasts or linear functions of EMMs, trends, and comparisons of slopes. Plots and other displays.
- [glmtoolbox](https://cran.r-project.org/web/packages/glmtoolbox/index.html): glmtoolbox: Set of Tools to Data Analysis using Generalized Linear Models. Set of tools for the statistical analysis of data using: (1) normal linear models; (2) generalized linear models; (3) negative binomial regression models as alternative to the Poisson regression models under the presence of overdispersion; (4) beta-binomial and random-clumped binomial regression models as alternative to the binomial regression models under the presence of overdispersion; (5) Zero-inflated and zero-altered regression models to deal with zero-excess in count data; (6) generalized nonlinear models; (7) generalized estimating equations for cluster correlated data.

# Introduction to the problem

Let's assume that you analyse a longitudinal clinical trial with a binary endpoint. _Longitudinal_ means that observations were made for each study subject (patient) multiple times, at subsequent timepoints (e.g. study visits).
Your goal is to formally compare the fraction (%) of some event between the treatment arms (groups) through a hypothesis test. For example, you may want to compare the % of (somehow defined) **clinical successes** between the new investigated treatment and some standard of care.
You may also want to adjust these estimates for some numerical covariates, test interactions, etc.

**But the problem is that the study is small **(say N=20 patients per arm) and also **the % of events is very small at certain timepoints** (like 1-5) **or even drops to zero** over time.

This creates the problem called a **quasi separation** (one value, 0 or 1 dominates in the response variable in certain group(s)) or **full separation** (all 0 or all 1 in the response variable in certain group(s)). 
In other words:
- quasi separation = some combination of predictors perfectly predicts the outcome in parts of the data
- full separation = a predictor completely determines the outcome (e.g., group A has only 0s, group B has only 1s).

In either case, the classic logistic regression may either be biased or not converge and you may notice warning or error messages thrown by the estimation procedure like: `iterations limit exceeded` or `fitted probabilities numerically 0 or 1 occurred`.

That's because for full separation the Maximum Likelihood Estimation (MLE) fails: the model tries to assign infinite-heading coefficients (log-odds -> ±∞) to match the perfect separation, which leads to non-convergence. For the quasi separation the situation ends up with unstable estimates for the predictors - when you look at the coefficients table, some estimates are infinite (or close to it) or even missing from the printout (depending on the implementation). The corresponding standard errors may be huge. You will notice it for sure!

There are various approaches to this problem, e.g. via: 
- exact logistic regression
- Firth-penalised logistic regression
- Bayesian logistic regression
- independent exact tests (Fisher, Boschloo, Barnard - depending on your case!); they don't allow for interactions and covariate adjustments.

We will focus on the Firth approach, but in an unusual manner - applied to GEE.
Why unusual? Because the penalisation is applied to the likelihood function. There is no such thing in GEE, which a semi-parametric method.
This is why we call it "Firth-like".

Let's have a look at an exemplary text in the statistical report:

> The statistical analysis plan (SAP) specified the use of GEE estimation for the logistic regression to account for the longitudinal nature of the study and provide population-averaged estimates.
> The observed data, however, presented challenges due to the low event rates and 0 events at later visits leading to partially and then complete separation issue, which caused convergence difficulties
>  and instability in parameter estimates under the planned approach.
> To address this, a Firth-like penalized GEE method was applied. It extends the ordinary GEE framework with penalties akin to Firth’s correction to achieve first-order bias reduction via
> penalized maximum likelihood using Jeffreys prior (Firth, 1993). To further improve the small-sample properties, the bias-corrected estimator of variance was employed, as proposed in (Mancl & DeRouen, 2001), and is also robust to the misspecification of the working correlation matrix.
>  This approach retains the intended marginal (population-level) interpretation while improving estimation stability in the presence of sparse binary outcomes,
>  offering a methodologically sound and robust alternative consistent with the original analytic intent.

Now, let me show you an example.

-----------------

First, the data:
```{r}
> dput(data)
data <- structure(list(PatientId = structure(c(1L, 1L, 1L, 1L, 1L, 2L, 
2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 
5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 
8L, 8L, 9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 10L, 11L, 11L, 
11L, 11L, 11L, 12L, 12L, 12L, 12L, 12L, 13L, 13L, 13L, 13L, 13L, 
14L, 14L, 14L, 14L, 14L, 15L, 15L, 15L, 15L, 15L, 16L, 16L, 16L, 
16L, 16L, 17L, 17L, 17L, 17L, 17L, 18L, 18L, 18L, 18L, 18L, 19L, 
19L, 19L, 19L, 19L, 20L, 20L, 20L, 20L, 20L, 21L, 21L, 21L, 21L, 
21L, 22L, 22L, 22L, 22L, 22L, 23L, 23L, 23L, 23L, 23L, 24L, 24L, 
24L, 24L, 24L, 25L, 25L, 25L, 25L, 25L, 26L, 26L, 26L, 26L, 26L, 
27L, 27L, 27L, 27L, 27L, 28L, 28L, 28L, 28L, 28L, 29L, 29L, 29L, 
29L, 29L, 30L, 30L, 30L, 30L, 30L, 31L, 31L, 31L, 31L, 31L, 32L, 
32L, 33L, 33L, 33L, 33L, 33L, 34L, 34L, 34L, 34L, 34L, 35L, 35L, 
35L, 35L, 35L, 36L, 36L, 36L, 36L, 36L), levels = c("05018801", 
"05018802", "05018803", "05018804", "05018805", "05018806", "05018807", 
"05018808", "05018809", "05018810", "05018811", "05018812", "05018813", 
"05018814", "05018815", "05018816", "05018817", "05018818", "05018819", 
"05018820", "05018821", "05018822", "05018823", "05018824", "05018826", 
"05018827", "05018828", "05018829", "05018830", "05018831", "05018832", 
"05018833", "05018834", "05018835", "05018836", "05018837"), class = "factor"), 
    Arm = structure(c(2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 
    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
    2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
    2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
    2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
    1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 
    2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
    2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 
    1L, 1L), levels = c("Active", "Control"), class = "factor"), 
    Visit = structure(c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
    1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 1L, 2L, 3L, 
    4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 
    4L, 5L), levels = c("Baseline", "T1", "T2", "T3", "T4"), class = "factor"), 
    Response = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 
    0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 0, 0, 0, 0, 1, 1, 0, 0, 0), Age = c(42L, 42L, 42L, 42L, 
    42L, 41L, 41L, 41L, 41L, 41L, 60L, 60L, 60L, 60L, 60L, 46L, 
    46L, 46L, 46L, 46L, 68L, 68L, 68L, 68L, 68L, 42L, 42L, 42L, 
    42L, 42L, 61L, 61L, 61L, 61L, 61L, 49L, 49L, 49L, 49L, 49L, 
    47L, 47L, 47L, 47L, 47L, 50L, 50L, 50L, 50L, 50L, 47L, 47L, 
    47L, 47L, 47L, 36L, 36L, 36L, 36L, 36L, 49L, 49L, 49L, 49L, 
    49L, 52L, 52L, 52L, 52L, 52L, 66L, 66L, 66L, 66L, 66L, 41L, 
    41L, 41L, 41L, 41L, 57L, 57L, 57L, 57L, 57L, 51L, 51L, 51L, 
    51L, 51L, 43L, 43L, 43L, 43L, 43L, 46L, 46L, 46L, 46L, 46L, 
    53L, 53L, 53L, 53L, 53L, 62L, 62L, 62L, 62L, 62L, 44L, 44L, 
    44L, 44L, 44L, 41L, 41L, 41L, 41L, 41L, 58L, 58L, 58L, 58L, 
    58L, 61L, 61L, 61L, 61L, 61L, 63L, 63L, 63L, 63L, 63L, 57L, 
    57L, 57L, 57L, 57L, 55L, 55L, 55L, 55L, 55L, 61L, 61L, 61L, 
    61L, 61L, 51L, 51L, 51L, 51L, 51L, 44L, 44L, 36L, 36L, 36L, 
    36L, 36L, 41L, 41L, 41L, 41L, 41L, 54L, 54L, 54L, 54L, 54L, 
    60L, 60L, 60L, 60L, 60L), Wave = c(1, 2, 3, 4, 5, 1, 2, 3, 
    4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 
    3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 
    2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 
    1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 
    5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 
    4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 
    3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 
    2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 1, 2, 3, 
    4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5)), class = "data.frame", row.names = c(NA, 
-177L), Version = 1L, Date = structure(1745839693.267, tzone = "UTC", class = c("POSIXct", 
"POSIXt")))
```

which looks like that:
``` {r}
> head(data)
  PatientId     Arm    Visit Response Age Wave
1  05018801 Control Baseline        1  42    1
2  05018801 Control       T1        1  42    2
3  05018801 Control       T2        0  42    3
4  05018801 Control       T3        0  42    4
5  05018801 Control       T4        0  42    5
6  05018802  Active Baseline        0  41    1
> 
> tail(data)
    PatientId     Arm    Visit Response Age Wave
172  05018836 Control       T4        0  54    5
173  05018837  Active Baseline        1  60    1
174  05018837  Active       T1        1  60    2
175  05018837  Active       T2        0  60    3
176  05018837  Active       T3        0  60    4
177  05018837  Active       T4        0  60    5
```

# 1) Summary of the event counts

Let's sumamrize the number of events:

```{r}
> data %>% 
     group_by(Arm, Visit, Response) %>% 
     count() %>% 
     ungroup() %>% 
     complete(Arm, Visit, Response, fill = list(n=0)) %>% 
     mutate(Ntot = sum(n), .by = c("Arm", "Visit")) %>% 
     mutate(n = sprintf("%d (%.1f%%)", n, 100*n/Ntot)) %>% 
     filter(Response == 1) %>% 
     select(-Ntot) %>% 
     pivot_wider(names_from = Arm, values_from = n) %>% 
     select(-Response)

# A tibble: 5 × 3
  Visit    Active     Control  
  <fct>    <chr>      <chr>    
1 Baseline 10 (55.6%) 6 (33.3%)
2 T1       10 (55.6%) 4 (22.2%)
3 T2       0 (0.0%)   1 (5.9%) 
4 T3       0 (0.0%)   0 (0.0%) 
5 T4       0 (0.0%)   0 (0.0%) 
```
Not only the number of events is low, but also drops to zero since the T2 visit.

# 2) Fitting the classic GEE logistic regression

Let's fit the classic logistic regression.
Notice, that for "unstructured" covariance it won't even complete, saying  `Error in solve.default(Vi) : system is computationally singular: reciprocal condition number = 6.25574e-17`.
So we need to simplify the structure to AR(1).

```{r}
> library(glmtoolbox)
> library(emmeans)

# This may take a few seconds - 500 iterations

> gee_model <- glmgee(formula = Response ~ Arm * Visit + Age + Age : Visit, 
                    id = PatientId,
                    data = data,
                    waves = data$Wave,
                    corstr = "AR-M-dependent(1)",
                    family = binomial(link = "logit"),
                    maxit = 500) # default = 50

Warning message:
Iteration limit exceeded!!
 
> summary(gee_model)

Sample size
   Number of observations:  177
       Number of clusters:  36 
                            Min  25%  50%  75%  Max
            Cluster sizes:    2    5    5    5    5
*************************************************************
Model
        Variance function:  binomial
            Link function:  logit
    Correlation structure:  AR-M-dependent(1)
*************************************************************
Coefficients
                    Estimate Std.Error   z-value Pr(>|z|)
(Intercept)        -7.78e+15  3.57e+15 -2.18e+00 0.029314
ArmControl         -9.50e+14  9.85e+14 -9.65e-01 0.334473
VisitT1             1.17e+16  5.68e+15  2.06e+00 0.039254
VisitT2             7.63e+15  3.59e+15  2.12e+00 0.033694
VisitT3            -2.79e+15  3.57e+15 -7.81e-01 0.434915
VisitT4            -4.67e+15  3.57e+15 -1.31e+00 0.190737
Age                 1.55e+14  6.38e+13  2.42e+00 0.015334
ArmControl:VisitT1 -1.97e+14  1.61e+15 -1.23e-01 0.902351
ArmControl:VisitT2  3.20e+15  1.02e+15  3.14e+00 0.001665
ArmControl:VisitT3  1.13e+15  9.85e+14  1.15e+00 0.251736
ArmControl:VisitT4  1.01e+15  9.85e+14  1.02e+00 0.306202
VisitT1:Age        -2.30e+14  9.81e+13 -2.35e+00 0.018934
VisitT2:Age        -2.18e+14  6.42e+13 -3.39e+00 0.000701
VisitT3:Age        -8.94e+12  6.38e+13 -1.40e-01 0.888544
VisitT4:Age         2.25e+13  6.38e+13  3.53e-01 0.724126
                                                         
Dispersion          8.06e+14                             
*************************************************************
Working correlation
      [1]    [2]    [3]    [4]    [5] 
[1]  1.000 -0.089  0.008 -0.001  0.000
[2] -0.089  1.000 -0.089  0.008 -0.001
[3]  0.008 -0.089  1.000 -0.089  0.008
[4] -0.001  0.008 -0.089  1.000 -0.089
[5]  0.000 -0.001  0.008 -0.089  1.000
```

Ouch! Look at the estimated coefficients (log-odds scale) and their standard errors.
**No surprise - at 0 events MLE just could NOT do any better.**

OK, let's switch to the probability scale
_Note: model-based/predicted conditional means for the Bernoulli distribution = per-group covariate-adjusted probability._

```{r}
> emmeans(gee_model, specs = ~Arm * Visit + Age + Age : Visit, regrid="response")
 Arm     Visit     Age prob     SE  df asymp.LCL asymp.UCL
 Active  Baseline 51.1    1 0.1610 Inf     0.685     1.315
 Control Baseline 51.1    0 0.1400 Inf    -0.274     0.274
 Active  T1       51.1    1 0.1490 Inf     0.708     1.292
 Control T1       51.1    0 0.1540 Inf    -0.302     0.302
 Active  T2       51.1    0 0.0015 Inf    -0.003     0.003
 Control T2       51.1    0 0.0579 Inf    -0.114     0.114
 Active  T3       51.1    0 0.0000 Inf     0.000     0.000
 Control T3       51.1    0 0.0001 Inf     0.000     0.000
 Active  T4       51.1    0 0.0000 Inf     0.000     0.000
 Control T4       51.1    0 0.0000 Inf     0.000     0.000

Covariance estimate used: robust 
Confidence level used: 0.95
```

Notice also the negative bounds of the confidence intervals (CIs).

Again - no surprise. These are Wald's CIs, whih are symmetric on the current scale.
The predicted probability was 0, the SE is non-zero, so the CIs must take some "space", and since it's symmetric - it "hooked" the negative side.

We can fix this by computing the CIs on the linear-predictor scale and then back-tranform them on the response scale by replacing `regrid="response"` into `"type="response"`.

# 3) Fitting the Firth-like penalised GEE logistic regression

Let's first fit a model without age adjustment, just to compare the obtained % vs. the crude ones calculated below from the raw data.
```{r}
> library(geessbin)

gee_firth_model <- geessbin(formula = Response ~ Arm * Visit,
                        id = PatientId,
                        data = data,
                        corstr = "unstructured",
                        repeated = Visit, 
                        beta.method = "PGEE",  # bias reduction by adding a Firth-type penalty
                        SE.method = "MD", maxitr = 500)      # Mancl and DeRouen bias-corrected estimator

> summary(gee_firth_model)
Call:
geessbin(formula = Response ~ Arm * Visit, data = data, id = PatientId, 
    corstr = "unstructured", repeated = Visit, beta.method = "PGEE", 
    SE.method = "MD", maxitr = 500)

Correlation Structure:  unstructured 
Estimation Method for Regression Coefficients:  PGEE 
Estimation Method for Standard Errors:  MD 

Coefficients:
                   Estimate Std.err      Z Pr(>|Z|)
(Intercept)           0.156   0.499  0.313 7.55e-01
ArmControl           -0.845   0.727 -1.162 2.45e-01
VisitT1               0.188   0.589  0.319 7.50e-01
VisitT2              -4.374   0.567 -7.709 1.27e-14
VisitT3              -4.345   0.567 -7.658 1.90e-14
VisitT4              -4.345   0.567 -7.658 1.90e-14
ArmControl:VisitT1   -0.635   1.001 -0.635 5.25e-01
ArmControl:VisitT2    2.543   1.216  2.092 3.65e-02
ArmControl:VisitT3    0.886   0.823  1.076 2.82e-01
ArmControl:VisitT4    0.886   0.823  1.076 2.82e-01

Odds Ratios with 95% Confidence Intervals :
                   Odds Ratio Lower Limit Upper Limit
ArmControl             0.4296     0.10329      1.7867
VisitT1                1.2064     0.38056      3.8245
VisitT2                0.0126     0.00415      0.0383
VisitT3                0.0130     0.00426      0.0394
VisitT4                0.0130     0.00426      0.0394
ArmControl:VisitT1     0.5297     0.07452      3.7653
ArmControl:VisitT2    12.7210     1.17362    137.8846
ArmControl:VisitT3     2.4261     0.48302     12.1856
ArmControl:VisitT4     2.4261     0.48302     12.1856

Estimated Scale Parameter:  0.517
Number of Iterations:  8 

Working Correlation:
         Baseline     T1      T2      T3      T4
Baseline   1.0000 0.3127 -0.2294 -0.0123 -0.0123
T1         0.3127 1.0000  0.5373  0.0132  0.0132
T2        -0.2294 0.5373  1.0000  0.0307  0.0307
T3        -0.0123 0.0132  0.0307  1.0000  0.0419
T4        -0.0123 0.0132  0.0307  0.0419  1.0000
```

Looks much better!
Let's switch to the probability scale and check if we obtained % comaprable to the crude ones computed before.

**Unfortunately, emmeans does not support package geessbin.**
We will have to use qdrg() to make it work.

```{r}
> refgrid <- qdrg(formula = Response ~ Arm * Visit,
                data = data,
                coef = coef(gee_firth_model),
                vcov = vcov(gee_firth_model),
                regrid="response",
                link = "logit")

> (emmean <- emmeans(refgrid, specs = ~ Arm * Visit)
 Arm     Visit    response     SE  df asymp.LCL asymp.UCL
 Active  Baseline    0.539 0.1240 Inf     0.296     0.782
 Control Baseline    0.334 0.1180 Inf     0.104     0.565
 Active  T1          0.585 0.1240 Inf     0.342     0.829
 Control T1          0.243 0.1040 Inf     0.039     0.447
 Active  T2          0.015 0.0036 Inf     0.007     0.022
 Control T2          0.075 0.0593 Inf    -0.042     0.191
 Active  T3          0.015 0.0037 Inf     0.008     0.022
 Control T3          0.016 0.0040 Inf     0.008     0.023
 Active  T4          0.015 0.0037 Inf     0.008     0.022
 Control T4          0.016 0.0040 Inf     0.008     0.023

Confidence level used: 0.95 
```

Good! Finally we can get some reasonly looking numbers!

Let's compare:
| Timepoint | Arm     | Crude (%) | Firth GEE (%) |
|-----------|---------|-----------|---------------|
| Baseline  | Active  | 55.6      | 53.9          |
| Baseline  | Control | 33.3      | 33.4          |
| T1        | Active  | 55.6      | 58.5          |
| T1        | Control | 22.2      | 24.3          |
| T2        | Active  | 0.0       | 1.5 (artifact)|
| T2        | Control | 5.9       | 7.5           |
(all others are just artifacts).

**I like these estimates!** considering so poorly conditioned dataset.

PS: Although no events were observed at late time points in the sample (i.e., complete separation),
the penalized GEE method used here applies a first-order bias correction that shrinks extreme estimates toward finite values.
As a result, even when the observed event count is zero, the model does not force the predicted probability to be exactly zero.
Instead, it returns a small, non-zero estimate reflecting both the finite-sample uncertainty and the underlying structure of the model (e.g., covariates, correlation).
This behavior is expected and desirable: it avoids the infinite log-odds estimates and degenerate inference that would arise from unpenalized models
under separation and provides more stable estimates and confidence intervals.

Thus, the reported probabilities (e.g., 0.015 instead of 0) and comparisons (0.001 instead of 0) represent a model-based shrinkage toward plausible values in the population,
not a contradiction of the observed data. **In other words, the estimate is not just a restatement of sample proportions - penalized regression is model-based, not descriptive.**

OK, now let's refit the model with age adjustment and do some between-arm comparison:

```{r}
> gee_firth_model <- geessbin(formula = Response ~ Arm * Visit + Age + Age : Visit,
                        id = PatientId,
                        data = data,
                        corstr = "unstructured",
                        repeated = Visit, 
                        beta.method = "PGEE",  # bias reduction by adding a Firth-type penalty
                        SE.method = "MD", maxitr = 500)      # Mancl and DeRouen bias-corrected estimator

> refgrid <- qdrg(formula = Response ~ Arm * Visit + Age + Age : Visit,
                data = data,
                coef = coef(gee_firth_model),
                vcov = vcov(gee_firth_model),
                regrid="response",
                link = "logit")

> emmean <- emmeans(refgrid, specs = ~ Arm * Visit + Age + Age : Visit)

> update(contrast(emmean,
                 list(                     # Bas   T1    T2    T3    T4
                   "Bas: ACT vs. CTL"  = c(-1,1,  0,0,  0,0,  0,0,  0,0),
                   "T1:  ACT vs. CTL"  = c( 0,0, -1,1,  0,0,  0,0,  0,0),
                   "T2:  ACT vs. CTL"  = c( 0,0,  0,0, -1,1,  0,0,  0,0),
                   "T3:  ACT vs. CTL"  = c( 0,0,  0,0,  0,0, -1,1,  0,0),
                   "T4:  ACT vs. CTL"  = c( 0,0,  0,0,  0,0,  0,0, -1,1)),
        adjust="mvt", 
        level = 0.95, 
        infer = c(TRUE, TRUE)))
 contrast         estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
 Bas: ACT vs. CTL   -0.205 0.1710 Inf    -0.630    0.2207  -1.197  0.6430
 T1:  ACT vs. CTL   -0.342 0.1620 Inf    -0.745    0.0607  -2.113  0.1300
 T2:  ACT vs. CTL    0.060 0.0594 Inf    -0.088    0.2078   1.009  0.7700
 T3:  ACT vs. CTL    0.001 0.0055 Inf    -0.013    0.0142   0.113  1.0000
 T4:  ACT vs. CTL    0.001 0.0055 Inf    -0.013    0.0142   0.113  1.0000

Confidence level used: 0.95 
Conf-level adjustment: mvt method for 5 estimates 
P value adjustment: mvt method for 5 tests 
```

...and that's it.
