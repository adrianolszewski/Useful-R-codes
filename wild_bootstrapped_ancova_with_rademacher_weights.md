# Simple analysis of covariance (ANCOVA) with nonnormal residuals and heteroscedasticity using the wild-type bootstrap method with Rademacher (or Mammen) weights and the robust HC3 variance estimator.

## Description of an exemplary case
We will consider a simple, two-arm, randomized, controlled clinical trial (RCT) with a (anticipated!) Gaussian response.

The key analysis in such studies is the comparison between arms of the mean response after drug administration, adjusted for the baseline response.

OK, but why is this adjustment to the baseline even necessary in an RCT?
Well, this is a purely technical matter.
You see, randomization does not guarantee balance of covariates between arms in a given case.
A fairly large, statistically significant imbalance is not surprising.

Although such a difference at baseline is purely a sampling artifact, not accounting for it may noticeably reduce power.
To address this, you may consider the _constrained longitudinal data analysis (cLDA)_ or _analysis of covariance (ANCOVA)*_.

/* ANCOVA, BTW, generalizes to MMRM in longitudinal studies >2 time points) /

However, the situation can get worse for three reasons:

<ol>
<li>Although the population distribution of the analyzed response may be assumed by you (domain knowledge, experts opinion, educated guess)
to be well-approximated by a Gaussian distribution, the sample does not necessarily have t reflect this.
This means that the residuals from the covariance model (ANCOVA) do not have to be "nicely" normal. Even close to that.

<sub>/ Now, whether you should care about this is NOT the issue of this document. You may use QQ plots or (90+) formal tests for non-normality to diagnose this, and when it shows a noticeable problem, you may want to address it,
or do nothing and wishfully rely on (pray to...)  the Central Limit Theorem (and believe in _"N>30"_).
You may do the sensitivity analysis with other methods as well. **Anyway - it's all up to you.** / </sub>

**IF** you decide to address these non-Gaussian residuals, the **Freedman-Lane Permutation Covariance Analysis** is at your disposal.
In R you can do it with the [permuco](https://cran.r-project.org/web/packages/permuco/index.html) package (functions: `aovperm()`, `lmperm()`)</li>

<li>Post-randomization imbalance can also be reflected by different variances across groups (= heterogeneity of variances across groups = heteroscedasticity of residuals).
To account for this, **Wild Bootstrap with Rademacher or Mammen weights applied to ANCOVA** can reduce Type-1 error.</li>
, 
<li> Last, but not least, high-leverage observations (related to point 2) can occur as well, so heteroscedasticity-consistent standard errors (HC),
preferably HC2 or – more conservatively – HC3, can be very helpful in the bootstrap covariance analysis.</li>
</ol>

**In this paper, I will implement this kind of analysis from scratch in R.**

**Note:** We will need two separate approaches:

1. to obtain a confidence interval for the between-arm difference, we need the inference over the **full model**. Full = with both the treatment effect and the baseline covariate.
(under the standard treatment coding the treatment effect is just the value of $\beta$ coefficient for the treatment arm).

2. to obtain the p-value for the between-arm differene, we need the inferene over the **restricted model**. Restricted = without the treatment effect (just the baseline and intercept).

-----

# Math + programming
Before I go to the code, let's only briefly sketch the solution.

## 1. The General Linear Model (ANCOVA)
I'm going to follow textbooks and define the full model (used for the point estimate and the confidence interval) as:
$$y = X\beta + \epsilon$$

Where:
* $y$ is the $N \times 1$ vector of responses.
* $X$ is the $N \times 3$ design matrix (Intercept, Baseline, Arm).
* $\beta$ is the vector of coefficients $[\beta_0, \beta_{baseline}, \beta_{arm}]^T$.
* $\epsilon$ is the vector of errors, where we assume $E[\epsilon] = 0$ but allow for heteroscedasticity $Var(\epsilon) = \Sigma$.

The **Ordinary Least Squares (OLS)** estimator is: $$\hat{\beta} = (X^T X)^{-1} X^T y$$

PS: Yes, I know that naive calculations using the textbok formulas are inferior to QR factorization.
But still, compared to a loop-based bootstrap, even this naive appraoch will save a LOT of your time, so don't nag...

## 2. Residual Transformation with HC3 (robust) estimator
To account for heteroscedasticity and leverage, we will transform the raw OLS residuals $e_i = y_i - \hat{y}_i$.

The **HC3-adjusted and centered** residuals used in the bootstrap are defined as:
``` math
\begin{equation}
\tilde{e}_i = \left( \frac{e_i}{1 - h_i} \right) - \bar{e}_{adj}
\hspace{2cm} (2)
\end{equation}
```
Where $h_i$ are the diagonal elements of the hat matrix $H = X(X^TX)^{-1}X^T$, representing the leverage of each observation.

And after that, we will center these residuals. YOu may ask why since for the RAW residuals E[ε]=0.
OK, but here we divide by ($1−h_i$​) which "de-centers" ("un-centers?") them, and according to [Long & Ervin 2000] we should "recenter" them to restore E[ε]=0.

## 3. Wild bootstrap in action

Now, let's generate $R$ bootstrap samples of the response variable. The logic differs between the p-value and the confidence interval, as I mentioned before:

1. **for the p-value** we need to simulate the world where the null hypothesis $H_0: \beta_{arm} = 0$ is true. 
``` math
\begin{equation}
y^*_r = \hat{y}_{null} + \tilde{e}_{null} \cdot w_r
\end{equation}
```
Where:
* $\hat{y}_{null}$ are the fitted values from the model *without* the treatment arm.
* $w_r$ is a vector of **Rademacher weights** $\in \{-1, 1\}$.

PS: optionally you could try also the Mammen weights (but should give very close results):

``` math
\begin{equation}
v_b^* = 
\begin{cases} 
1 - \psi & \text{with probability } \frac{\psi}{\sqrt{5}} \\
\psi & \text{with probability } 1 - \frac{\psi}{\sqrt{5}}
\end{cases}
\end{equation}
```
where $\psi = \frac{1 + \sqrt{5}}{2} \approx 1.618$ is the _golden ratio_.
These weights ensure $\mathbb{E}[v\_b^{\ast}] = 0$ and $\mathbb{E}[(v\_b^{\ast})^2] = 1$, providing higher-order "refinements" over Rademacher weights (±1 with equal probability) in heteroskedastic settings.

2. **for the confidence interval** we need the full model (unrestricted) with which we simulate the world based on our observed effect:
``` math
\begin{equation}
y^*_r = \hat{y}_{full} + \tilde{e}_{full} \cdot w_r
\end{equation}
```
## 4. Bootstrapped coefficients (trick with manual algebra)
Instead of re-running the full optimization 10 000 times in a loop, I use pre-calculated projection matrix.
For each bootstrap iteration $r$:

```math
\begin{equation}
\hat{\beta}^*_r = \mathbb{P} y^*_r \quad \text{where} \quad \mathbb{P} = (X^T X)^{-1} X^T
\end{equation}
```

and then I extract the treatment effect:
```math
\begin{equation}
\hat{\beta}_{arm, r}^* = \hat{\beta}^*_r[index_{arm}]
\end{equation}
```

## 5. Inference
### p-value
Using the **Davison-Hinkley** correction I ensure the p-value is never exactly zero:
```math
\begin{equation}
p = \frac{\sum_{r=1}^R \mathbb{I}(|\hat{\beta}^*_{arm, r}| \ge |\hat{\beta}_{arm}|) + 1}{R + 1}
\end{equation}
```

Where $\mathbb{I}(\cdot)$ is the indicator function.

### confidence interval
I use the **percentile** approach on the unrestricted distribution:
```math
\begin{equation}
CI_{95\%} = [Q^*(0.025), Q^*(0.975)]
\end{equation}
```
Where $Q^{\ast} \text{ is the quantile function of the } \hat{\beta}^{\ast}_{\text{full}} \text{ distribution.}$

---------

# The code

## Generating sample data
```r
library(dplyr)

set.seed(12345)
N <- 100
df <-
  data.frame(baseline = rnorm(N, mean = 50, sd = 10),
             arm      = rbinom(N, 1, 0.5)) %>% 
  mutate(errors   = ifelse(arm == "Group B", rnorm(N, 0, 15), rnorm(N, 0, 5)),
         response = 10 + 0.8 * baseline + 1.4 * arm + errors,
         arm      = factor(arm, levels = c(0, 1), labels = c("Group A", "Group B")))
```

```r
> head(df)
  baseline     arm     errors response
1 55.85529 Group B -8.0966416 47.98759
2 57.09466 Group B  2.7419896 59.81772
3 48.90697 Group A  0.9764108 50.10198
4 45.46503 Group B -4.0324900 43.73953
5 56.05887 Group B -0.5431212 55.70398
6 31.82044 Group B -1.2547331 35.60162
```

## Fitting the required models
```r
fit_full <- lm(response ~ arm + baseline, data = df)
fit_null <- lm(response ~ baseline, data = df)

beta_hat <- coef(fit_full)["armGroup B"]
```

## Confidence interval
```r
R <- 10000
set.seed(12345)

X_full  <- model.matrix(fit_full)
XtX_inv <- solve(t(X_full) %*% X_full)
Xt      <- t(X_full)
arm_idx <- which(colnames(X_full) == "armGroup B")

# Rademacher weights (common to both bootstraps)
W <- matrix(sample(c(-1, 1), N * R, replace = TRUE), nrow = N, ncol = R)

h_null      <- hatvalues(fit_null)
res_null    <- residuals(fit_null) / (1 - h_null)
res_null    <- res_null - mean(res_null)

Y_star_null <- fitted(fit_null) + res_null * W
beta_boot_null <- (XtX_inv %*% Xt %*% Y_star_null)[arm_idx, ]

# The Davison-Hinkley finiteness correction
p_value     <- (sum(abs(beta_boot_null) >= abs(beta_hat)) + 1) / (R + 1)
```

## p-value
```r
h_full         <- hatvalues(fit_full)
res_full       <- residuals(fit_full) / (1 - h_full)
res_full       <- res_full - mean(res_full)
Y_star_full    <- fitted(fit_full) + res_full * W
beta_boot_full <- (XtX_inv %*% Xt %*% Y_star_full)[arm_idx, ]

CI_percentile <- quantile(beta_boot_full, probs = c(0.025, 0.975))
```

And the result:
```r
data.frame( estimate   = beta_hat,
            CI_lower   = CI_percentile[1],
            CI_upper   = CI_percentile[2],
            p_value    = p_value)

           estimate    CI_lower CI_upper   p_value
armGroup B 1.750053 -0.08306756 3.578945 0.0659934
```

