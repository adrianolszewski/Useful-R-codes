Despite a common claim that Mann-Whitney (Wilcoxon), Kruskal-Wallis, Fligner-Poicello and other rank-based tests compare medians,
**this is a widespread misconception when considered in general**, i.e. unless strong IID condition is met, resulting in the location-shift case, consistent with the proportional-odds model.

This  misconception has been spread so widely, so many textbooks, free and paid courses, even academic lecturers repeat it.

This is by-definition (by construction) a test of stochastic equivalence vs superiority, i.e. whether P(B>A) + 0.5*P(B=A) = 0.5.
Such superiority may occur for a difference in any properties of the compared samples: 
- **locations** (means, medians, Hodges-Lehman 2-sample pseudo-median and other measures),
- **dispersions** (variances, general scale parameter, range),
- **shapes** (skewness, kurtosis)

In other words, **when a difference is found with test between the compared samples, it cannot be simply attributed to either location, dispersion or shape, unless 2 of the 3 characteristics are comparable**. Otherwise each of them may "contribute" to stochastic superiority, confusing the analyst.

And then analysts may be surprised how is that possible that for exactly equal means or medians the test returns p-value < 0.0...01 (for quite small samples) or for much different means or medians the p-value > 0.999, while methods like Brown-Mood or quantile regression (under a variety of methods for obtaining standard errors) differ greatly in their findings, not to mention that visual assessments (e.g. box-plots) also do not support respective claims.

The following scenarios may easily happen:

1. Equal medians or means, H0 rejected (even with p<0.0...01 at small samples)
2. Different medians, H0 not rejected (even with p>0.999... at large samples)
3. Different entire empirical distributions (eCDFs), H0 not rejected (in the figure below it's the middle case)
<img width="823" height="360" alt="obraz" src="https://github.com/user-attachments/assets/107ceec2-d04a-4d08-b4bc-568b0836a5fa" />

There are several articles confirming that these tests fail as tests of medians*, but only a few books explain it. So, in case you need need, let me cite a few, with the most important one opening the list:

1.    Brunner, E., Bathke, A. C., & Konietschke, F. (2018). Rank and pseudo-rank procedures for independent observations in factorial designs: Using R and SAS. Springer.
2.    Hettmansperger, T. P., & McKean, J. W. (2010). Robust nonparametric statistical methods (2nd ed.). CRC Press. (chapter 2.4: Inference based on the Mann-Whitney-Wilcoxon)
3.    Thas, O. (2010). Comparing distributions. Springer. (chapter 9.3.3.1: Implied Null Hypothesis)
4.    Wilcox, R. R. (2021). Introduction to robust estimation and hypothesis testing (5th ed.). Academic Press. (chapter 5.7 Methods Based on Ranks and the Typical Difference)
5.    Fagerland MW, Sandvik L. The Wilcoxon-Mann-Whitney test under scrutiny. Stat Med. 2009 May 1;28(10):1487-97. doi: 10.1002/sim.3561. PMID: 19247980. [ https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/paranp?action=AttachFile&do=get&target=MannW.pdf ]

Plus the original paper of Mann and Whitney:

6.    Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. The Annals of Mathematical Statistics, 18(1), 50–60. https://doi.org/10.1214/aoms/1177730491

Plus some papers and discussions:
- Divine, G. W., Norton, H. J., Barón, A. E., & Juarez-Colunga, E. (2018). The Wilcoxon–Mann–Whitney Procedure Fails as a Test of Medians. The American Statistician, 72(3), 278–286. https://doi.org/10.1080/00031305.2017.1305291 [
Article The Wilcoxon–Mann–Whitney Procedure Fails as a Test of Medians
]
- Conroy, R. M. (2012). What hypotheses do “nonparametric” two-group tests actually test? The Stata Journal, 12(2), 182–190. https://doi.org/10.1177/1536867X1201200202 [
Article What Hypotheses do “Nonparametric” Two-Group Tests Actually Test?
]
- Hart A. (2001). Mann-Whitney test is not just a test of medians: differences in spread can be important. BMJ (Clinical research ed.), 323(7309), 391–393. https://doi.org/10.1136/bmj.323.7309.391 [
Article Mann-Whitney test is not just a test of medians: differences...
]
- Kleinman K, Example 2014.6: Comparing medians and the Wilcoxon rank-sum test [ http://proc-x.com/2014/06/example-2014-6-comparing-medians-and-the-wilcoxon-rank-sum-test/ ]
- https://stats.stackexchange.com/questions/363335/wilcoxon-signed-rank-test-null-hypothesis-statement
Plus some toy figures from my various presentations.

Follow the discussion at https://www.researchgate.net/post/Mann-Whitney_Wilcoxon_rank_test_the_null_hypothesis_not_about_medians-in_case_you_needed_the_references?_init=1

---

# Scenario 1: Equal medians, H0 rejected
```r
set.seed(1000)
#N <- 10  # p-value: 0.03247
#N <- 20  # p-value: 0.00625
#N <- 30  # p-value: 7.882e-05
N  <- 50  # p-value: insanely small

x1   <- runif(N, 2, 5)
x2   <- runif(N, 4, 5)
x3   <- runif(N, 3, 5)

x11  <- runif(N, 5.1, 6)
x21  <- runif(N, 5.1, 10)
x31  <- runif(N, 5.1, 8)

x111 <- sort(c(x1, 5, x11))
x222 <- sort(c(x2, 5, x21))
x333 <- sort(c(x3, 5, x31))

dat <-data.frame(response = c(x111, x222, x333), 
                 group=rep(LETTERS[1:3], each = 2*N+1))
```

```r
dat %>% 
    mutate(median = median(response), .by="group") %>% 
    ggplot(aes(x=response)) + 
    geom_histogram(aes(x=response), bins=nclass.Sturges) + 
    facet_wrap(~group) +
    theme_bw() + 
    geom_vline(aes(xintercept=median), col="red") +
    geom_rug(sides = "b", col="grey40") +
    geom_boxplot(aes(y=-5), width=3, show.legend = F) +
    ylab("count")
```
<img width="886" height="507" alt="obraz" src="https://github.com/user-attachments/assets/9059a259-24b9-423c-8a84-6c90fe05e039" />


Are all medians equal?
_Note: not just "approximate". `Equal` means `equal`_
```r
> unique(tapply(dat$response, dat$group, median))  # should see A: 5

A 
5 
```

Kruskal-Wallis (for N=50):
```r
> kruskal.test(response ~ group, data=dat)

	Kruskal-Wallis rank sum test

data:  response by group
Kruskal-Wallis chi-squared = 18.897, df = 2, p-value = 7.882e-05
```

If you have `tidyplots` installed:
https://tidyplots.org/

```r
library(dplyr)
library(emmeans)
library(quantreg)
library(patchwork)
library(tidyplots)

plot_superiority <- function(X, Y, name_X, name_Y) {
  
  expand_grid(X, Y) %>% 
    mutate(status = case_when(Y > X ~ "superior",
                              Y < X ~ "inferior",
                              .default = "equal")) -> superior_status
  
  prob_index <- (sum(superior_status$status == "superior") + 0.5*sum(superior_status$status == "equal")) / nrow(superior_status)
  
  med_X <- median(X)
  med_Y <- median(Y)
  
  rq_mod <- rq(response ~ group, data = dat, tau = 0.5)
  rq_p <- emmeans(rq_mod, specs = ~group, se = "boot") %>%
    pairs(reverse = TRUE) %>%
    summary(infer = TRUE) %>%
    pull(p.value) %>%
    scales::pvalue(add_p = TRUE)
  
  mwu_p <- wilcox.test(X, Y)$p.value %>% scales::pvalue(add_p = TRUE)
  
  plot_title <- sprintf("Stochastic superiority: P(%s>%s)+0.5*P(%s=%s) = %.1f%%", name_Y, name_X, name_Y, name_X, 100 * prob_index)
  plot_sub <- sprintf("quantile regression: %s | Mann-Whitney: %s", rq_p, mwu_p)
  
  ggplot(superior_status) +
    geom_segment(aes(x = Y, xend = X, y = 1, yend = 0, color = status), linewidth = 0.5, alpha=0.1) +
    geom_point(aes(x = X, y = 0), color = "black", size = 1) +
    geom_point(aes(x = Y, y = 1), color = "black", size = 1) +
    
    annotate(geom = "point", x = med_X, y = 0, color = "red", size = 3) +
    annotate(geom = "point", x = med_Y, y = 1, color = "red", size = 3) +
    
    annotate(geom = "segment",
             x = med_X, xend = med_Y, y = 0, yend = 1,
             color = "blue", linewidth = 1.5) +
    
    annotate(geom="text", 
             x = med_X, y = -0.05,
             label = sprintf("Median=%.1f", med_X),
             color = "blue", size = 4) +
    
    annotate(geom="text", 
             x = med_Y, y = 1.05,
             label = sprintf("Median=%.1f", med_Y),
             color = "blue", size = 4) +
    
    scale_color_manual(values = c("superior" = "#00BA38", "inferior" = "#F8766D", "equal" = "grey")) +
    scale_x_continuous(name = name_X, sec.axis = dup_axis(name = name_Y)) +
    labs(title = plot_title, subtitle = plot_sub) +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
    theme_bw() +
    theme(plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 13),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
}

patchwork::wrap_plots(
  plot_superiority(x111, x222, name_Y = "B", name_X = "A"), 
  plot_superiority(x222, x333, name_Y = "C", name_X = "B"),
  plot_superiority(x111, x333, name_Y = "C", name_X = "A"),

  # boxplots
  (dat %>% 
    mutate(median = median(response), .by="group") %>% 
    tidyplot(x = group, y = response) %>% 
    theme_tidyplot(fontsize = 11) %>% 
    add_boxplot() %>% 
    add_title(sprintf("Kruskal-Wallis p-value = %s", 
                      scales::pvalue(kruskal.test(response ~ group, data=dat)$p.value, 
                                     accuracy = 0.001, 
                                     add_p = T))) %>% 
    add_caption(paste("Each sample size =", N)) %>% 
    add_test_pvalue(method = "wilcox_test", label="Mann-Whitney (Wilcoxon) p-value={p.format}", label.size = 4, hide_info = TRUE) %>% 
    adjust_size(width = NA, height = NA) +
    geom_text(aes(y=median, label=sprintf("Med=%.5f", median)), check_overlap = TRUE, col="red") +
    scale_y_continuous(limits=c(0, 15))),

  ncol = 2)

```
<img width="1500" height="1028" alt="obraz" src="https://github.com/user-attachments/assets/9bdb7816-2fc3-41af-947c-10e5f2ba6805" />

Stochastic superiority: explains it all...
```r
> prob_index <- function(x, y) mean(outer(x, y, ">") + 0.5 * outer(x, y, "=="))

> setNames(combn(split(dat$response, dat$group), 2, function(x) prob_index(x[[1]], x[[2]])),
           nm = combn(c("A","B","C"), 2, paste, collapse=" vs. "))

  A vs. B   A vs. C   B vs. C 
0.3022743 0.3756985 0.6314577
```

If you have `brunnermunzel` installed:
https://cran.r-project.org/web/packages/brunnermunzel/

```r
> setNames(combn(split(dat$response, dat$group), 2, 
                 function(x) brunnermunzel::brunnermunzel.test(x[[2]], x[[1]])$estimate),
           nm = combn(c("A","B","C"), 2, paste, collapse=" vs. "))

  A vs. B   A vs. C   B vs. C 
0.3022743 0.3756985 0.6314577

# or, to match the plots descriptions:
> sprintf("%.1f%%", 100*(1-c(0.3022743, 0.3756985, 0.6314577)))
[1] "69.8%" "62.4%" "36.9%"
```

Explanation, books, papers, examples:
https://www.researchgate.net/post/Mann-Whitney_Wilcoxon_rank_test_the_null_hypothesis_not_about_medians-in_case_you_needed_the_references

---

# Scenario 2: Different medians, H0 not rejected
Under construction

In the meantime, look at this pic:
<img width="1349" height="867" alt="obraz" src="https://github.com/user-attachments/assets/4892128d-ff8a-41e9-99ff-57013a395922" />
