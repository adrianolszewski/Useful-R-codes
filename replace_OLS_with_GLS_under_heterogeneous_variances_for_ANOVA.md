**GLS-fit general linear model handles unequal variances.**
No need to switch to non-parametric methods, like Kruskal-Wallis (suffering from unequal variances even more - look at the botom!)

``` r
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(nlme)
library(emmeans)
```

# Generate data
``` r
set.seed(1000)

d <- data.frame(value = c(rnorm(n = 50, mean = 1, sd = 1),
                          rnorm(n = 50, mean = 1.92, sd = 10),
                          rnorm(n = 50, mean = 5, sd = 5)),
                   group = gl(3, 50))

ggplot(d, aes(y=value, x=group)) + 
  geom_boxplot() +
  theme_bw()
```
<img width="1185" height="707" alt="obraz" src="https://github.com/user-attachments/assets/ebdbe8c7-6775-47a0-b868-614ab8b9071d" />

# Check variances (Brown-Forsythe)
``` r
car::leveneTest(value ~ group, data=d) # by default uses medians as the centers -> Brown-Forsythe
```
```
## Levene's Test for Homogeneity of Variance (center = median)
##        Df F value    Pr(>F)    
## group   2  39.113 2.401e-14 ***
##       147                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
# Fit a linear model using GLS estimation
``` r
m_gls <- gls(value ~ group, 
             weights = varIdent(form = ~ 1 | group), 
             data = d)
```

# Do a pairwise comparison
No adjustments for multiple comparsons - we want to see the raw numbers to see what's going on
We also use the Satterthwaite degrees of freedom, just like the Welch t test does.

``` r
emm_result <- update(pairs(
  emmeans(m_gls, specs = ~group, mode = "satterthwaite"), adjust="none"), 
  infer = c(TRUE, TRUE)) %>% 
  data.frame()

emm_result %>% 
mutate(across(where(is.numeric), ~sprintf("%.3f", .))) %>% 
    select(-SE)
```
```
##          contrast estimate     df lower.CL upper.CL t.ratio p.value
## 1 group1 - group2   -2.993 49.760   -6.007    0.022  -1.994   0.052
## 2 group1 - group3   -4.501 53.101   -5.818   -3.185  -6.857   0.000
## 3 group2 - group3   -1.509 66.534   -4.758    1.740  -0.927   0.357
```
``` r
emm_result %>% 
 ggplot(aes(x = estimate, y = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL), color = "#2c3e50", size = 0.8) +
    labs(title = "Pairwise Group Comparisons (GLS)",
subtitle = "95% Confidence Intervals (Unadjusted)",
x = "Estimated Difference in Means, y = "Contrast) +
    theme_minimal(base_size = 13) +
    theme(
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.title = element_text(face = "bold"))
```
<img width="692" height="662" alt="obraz" src="https://github.com/user-attachments/assets/b3f1fa28-8687-4c14-9eb9-5f3d3d5afba4" />

# Compare with classic Welch (-Satterthwaite) version of the t test
``` r
combn(levels(d$group), 2, simplify = FALSE) %>% 
  map_dfr(~ {
    g1 <- .x[1]
    g2 <- .x[2]
    
    sub_data <- d %>% filter(group %in% c(g1, g2))
    
    t.test(value ~ group, data = sub_data, var.equal = FALSE) %>% 
      broom::tidy(conf.int = TRUE) %>% 
      mutate(contrast = sprintf("group%s - group%s", g1, g2))
  }) %>% 
  select(contrast, estimate, df = parameter, lower.CL = conf.low, upper.CL = conf.high, t.ratio = statistic, p.value) %>% 
  mutate(across(where(is.numeric), ~sprintf("%.3f", .)))
```

```
## # A tibble: 3 × 7
##   contrast        estimate df     lower.CL upper.CL t.ratio p.value
##   <chr>           <chr>    <chr>  <chr>    <chr>    <chr>   <chr>  
## 1 group1 - group2 -2.993   49.760 -6.007   0.022    -1.994  0.052  
## 2 group1 - group3 -4.501   53.102 -5.818   -3.185   -6.857  0.000  
## 3 group2 - group3 -1.509   66.534 -4.758   1.740    -0.927  0.357
```

# ANOVA using the GLS-fit model and Satterthwaite DF
``` r
# ANOVA
joint_tests(m_gls)
```
```
##  model term df1   df2 F.ratio p.value
##  group        2 49.76  25.263 <0.0001
```

# Can also use the naive approach from car::Anova() using either infinite (asymptotic) DF or residual DF

## Residual DF
``` r
joint_tests(m_gls, mode="df.error")
```
```
##  model term df1 df2 F.ratio p.value
##  group        2 145  25.263 <0.0001
```
``` r
car::Anova(m_gls, test.statistic = "F", error.df = 145)
```
```
## Analysis of Deviance Table (Type II tests)
## 
## Response: value
##            Df      F    Pr(>F)    
## group       2 25.263 3.862e-10 ***
## Residuals 145                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
## Asymptotic DF
``` r
joint_tests(m_gls, mode="asymptotic")
```
```
##  model term df1 df2 F.ratio  Chisq p.value
##  group        2 Inf  25.263 50.526 <0.0001
```
``` r
car::Anova(m_gls)
```
```
## Analysis of Deviance Table (Type II tests)
## 
## Response: value
##       Df  Chisq Pr(>Chisq)    
## group  2 50.527  1.067e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

# Check the residuals. Look what might happen when you used the OLS naively
``` r
rbind(data.frame(Residuals = "OLS applied naively: response residuals", res = residuals(m_gls, type="response")),
      data.frame(Residuals = "GLS in action: normalized residuals", res= residuals(m_gls, type="normalized"))) %>% 
    
    ggplot(aes(sample = res)) + 
    facet_wrap(~ Residuals, scales = "free") +
    qqplotr::stat_qq_band(bandType = "ts") + # Aldor-Noiman tail-sensitive simultaneous CI band
    qqplotr::stat_qq_point() +
    qqplotr::stat_qq_line(col='red') +
    theme_bw() +
    labs(caption = "Aldor-Noiman et al. tail-sensitive simultaneous confidence bands")
```
<img width="1185" height="707" alt="obraz" src="https://github.com/user-attachments/assets/6091a838-38a6-4494-bf74-bb3215b00432" />



---
Kruskal-Wallis and Mann-Whitney (Wilcoxon) are sensitive to unequal variances creating stochastic superiority:
<img width="899" height="738" alt="obraz" src="https://github.com/user-attachments/assets/29bbd1a3-6332-4f7b-9682-44ba7c98c154" />
**And remember - neither of them compares medians in general!** ==> https://www.researchgate.net/post/Mann-Whitney_Wilcoxon_rank_test_the_null_hypothesis_not_about_medians-in_case_you_needed_the_references?_init=1
