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
library(tidyplots)

dat %>% 
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
    scale_y_continuous(limits=c(0, 15))
```
<img width="742" height="737" alt="obraz" src="https://github.com/user-attachments/assets/e4b91439-ae0e-4348-97c2-9ef854b3eb2e" />

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
           nm = combn(c("A","B","C"), 2, paste, collapse="-"))

      A-B       A-C       B-C 
0.3022743 0.3756985 0.6314577 
```

Explanation, books, papers, examples:
https://www.researchgate.net/post/Mann-Whitney_Wilcoxon_rank_test_the_null_hypothesis_not_about_medians-in_case_you_needed_the_references
