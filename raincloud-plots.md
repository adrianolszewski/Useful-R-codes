These are the plots I use daily at work to investgate my data.

_I'm not happy to add the Shapiro-Wilk test added (I could the remaining 90+ and watch them contradicting each other), having the awesome Aldor-Noiman QQ bands, but -well- in my field statistical reviewers are not always familiar with the modern tools or don't accept QQ plots, and blindly stick to pointless normality tests at N > 1000...
I can do nothing with that, so I add the non-normality p-value and pretend I don't see it ;]_

# Experiment 1: 2-sample Raincloud plot with QQ plot.
## Version A) Jittered raw data. Suitable for continuous variables.
```r
library(ggpp)
library(gghalves)
library(ggplot2)
library(qqplotr)
library(patchwork)

stats <- d %>%
  group_by(Group) %>%
  summarize(across(
    MyColumn,
    list(Mean   = ~ mean(.x, na.rm = TRUE),
         SD     = ~ sd(.x, na.rm = TRUE),
         Min    = ~ min(.x, na.rm = TRUE),
         Max    = ~ max(.x, na.rm = TRUE),
         Median = ~ quantile(.x, type = 3, prob = 0.5, na.rm = TRUE),
         N      = ~ sum(!is.na(.x)),
         Sk     = ~ e1071::skewness(.x, na.rm = TRUE),
         SW     = ~ ifelse(n() > 3, 
                           sprintf("Shapiro-Wilk: %s", scales::pvalue(shapiro.test(.x)$p.value, add_p = TRUE)), 
                           "NA")),
    .names = "{.fn}")) %>%
  mutate("MyColumn" = NA) # needed later by the QQ plot

p_raincloud <-
  d %>%
  ggplot(aes(x = Group, y = MyColumn, group = Group)) +
  theme_bw() +
  geom_half_boxplot(nudge = 0.05, width = 0.3, colour = "grey20", outlier.colour = "red3", outlier.alpha = .3) +
  geom_half_violin(side = "r", nudge = 0.05, trim = TRUE, alpha = 0.7, color = "gray85", fill = "gray85", scale = "width") +
  geom_pointrange(data = stats,
                  aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD),
                  size = .3, color = "blue") +

  geom_jitter(alpha = 0.25, size = 1.1, color = "grey30",
              position = position_jitternudge(width = 0.05, height = 0, x = -0.25, nudge.from = "jittered", seed = 1000)) +

  labs(title = sprintf("Observed values of [Age at onset of symptoms]"),
       caption = "The blue points and whiskers represent mean ± SD") +
  
  theme(plot.caption = element_text(face = "italic")) +
  ylab("Age at onser of symptoms [years]") +
  xlab(NULL) +
  geom_label(data = stats,
             aes(y = 1.25 * max(Max),
                 label = sprintf("N         = %d\nMean (SD) = %.1f (%.1f)\nMedian    = %.1f\nSkewness  = %.2f",
                                 N, Mean, SD, Median, Sk)),
             size = 2.8, vjust = "top", hjust = "left",
             family = "mono", fontface = "bold", 
             position = position_nudge(x = -.25),
             border.color = NA) +
  
  theme(panel.grid.minor.y = element_blank())

p_qq <-
  d %>%
  ggplot(aes(sample = MyColumn)) +
  theme_bw() +
  qqplotr::stat_qq_band(bandType = "ts", B = 500) +
  qqplotr::stat_qq_point() +
  qqplotr::stat_qq_line(col = "red") +
  facet_wrap(~Group, ncol = 1, scales = "free") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  labs(
    title = "Normal distribution QQ plot",
    caption = "Aldor-Noiman tail-sensitive\nsimultaneous confidence bands") +
  theme(plot.caption = element_text(face = "italic", size = 8)) +
  geom_label(data = stats,
             aes(x = -Inf, y = Inf, label = SW),
             hjust = "left", vjust = "top",
             border.colour = NA, label.padding = unit(5, "mm"), fill = NA)

wrap_plots(list(p_raincloud, p_qq), nrow = 1, widths = c(1, 0.5))
```
<img width="1309" height="707" alt="obraz" src="https://github.com/user-attachments/assets/3a303fc3-1bdd-4ef1-adf4-44131211a576" />

---
## Version B) Dot-plot. Especially suitable for discrete data, like drug doses (5, 10, 20mg), scores (from questionnaires), stages, counts, any integers
Here let me just reuse the current data.

Replace 
```r
  geom_jitter(alpha = 0.25, size = 1.1, color = "grey30",
              position = position_jitternudge(width = 0.05, height = 0, x = -0.25, nudge.from = "jittered", seed = 1000)) +
```
with
```r
  geom_dotplot(binaxis = "y", stackdir = "down",
               method = "dotdensity", stackratio = .7,
               position = position_nudge(x = -0.3),
               dotsize = .2, fill="grey60", col="grey60") +
```
Remember to adjust the parameter to fit your data best. There are NO universally good settings, forget it.

<img width="1309" height="707" alt="obraz" src="https://github.com/user-attachments/assets/f5f066ca-f765-460f-b798-bc4b837292cf" />
---
## Version C) Hexagonal binning

Replace 
```r
  geom_jitter(alpha = 0.25, size = 1.1, color = "grey30",
              position = position_jitternudge(width = 0.05, height = 0, x = -0.25, nudge.from = "jittered", seed = 1000)) +
```
with
```r
  geom_hex(bins=30, position = position_jitternudge(width = 0.05, height = 0, x = -0.25, nudge.from = "jittered", seed = 1000)) +
  scale_fill_gradient2(low = "grey80",mid = "orange", high = "darkblue", midpoint = 10) +
```
As before, adjust the parameters, colours to your liking.

<img width="1309" height="707" alt="obraz" src="https://github.com/user-attachments/assets/e3145987-3aa7-4d3a-8d7d-63fd285f8dbb" />

---
Data for reproduction:
```r
d <- structure(list(MyColumn = c(NA, 133.248414966671, 65.6263807967717, 
68.5547320287609, 149.46477879798, 60.9254681269893, 29.4397483174843, 
84.2142679051708, 82.9325876522646, 291.606392707816, 110.723844495154, 
85.8759995002801, NA, NA, 118.859908763686, NA, 31.4959169118109, 
252.707537584091, 85.4064858800037, 33.7494535178028, 89.9411164626721, 
109.021430427276, 116.771159803082, 46.335882155018, 181.055867438463, 
43.673727532713, NA, 69.9615056070868, 52.0187131557384, 208.569164395173, 
96.7438877891421, 89.0270051666986, 259.98128038573, NA, 123.081904318808, 
70.5665395055162, 22.6560404575816, 190.064736475845, 68.8817841948936, 
29.6548768360108, 80.9902547783252, 84.592554350053, 139.871080707866, 
85.7276299226071, 293.985060011295, 24.1963990616517, 75.369696591928, 
44.7230898471558, 43.945749742496, NA, NA, 108.881669392248, 
NA, 155.8419489317, 179.232059893707, 41.4312575063794, 40.4479982713193, 
98.862016144744, NA, 23.0059288214805, 59.7977734374446, 56.2716580713624, 
121.10151194823, 124.244666674814, 210.253013809448, 44.463673847004, 
42.5155286044907, 29.9076588050531, 96.2772279041376, 223.985146111338, 
63.5606938495985, 125.208094878638, 221.169651533188, 141.146450914228, 
147.515411552122, 103.235049821985, 43.5392434117552, 214.115712918944, 
NA, 21.0415746147644, 153.907183331553, 83.7905962970129, 81.5985279422507, 
58.1415056551535, 184.165755059588, 68.5295600393381, NA, 84.9448297421627, 
NA, 165.709208692469, 145.491810854812, 29.9567440117256, 120.519686986377, 
51.3694338085935, 118.303336844932, 84.3670139223103, 30.972294020509, 
128.30797215838, 66.4039035280396, 13.6485889393066, 65.1444716284445, 
105.965016868487, 101.724422332994, 112.558051476091, 113.096992264619, 
91.1907156119112, 29.3904870778818, 37.6806552446004, 67.471117577829, 
NA, 183.380808613267, 61.9106042908996, 210.493742448012, 111.107362353707, 
228.075221278102, 79.2142117080102, NA, 176.606822109198, 80.7640671930192, 
46.4786583410588, 84.0177521127365, 121.134922696973, 124.97269135368, 
106.994413385049, 203.630519306911, 76.6797452733054, 36.1669863507099, 
82.4944981149849, 104.436315014658, 208.707706714518, 104.04516072467, 
69.8990693606931, 243.822822052523, 86.6904764525719, 202.274137330135, 
45.3054487493892, 28.6720934622898, 99.2361205675185, 65.554580038236, 
NA, 67.695592570692, 134.762601500383, 69.7236144740728, 134.27111606883, 
144.978933704201, 74.1862625594999, 31.0315021221422, 34.2093521060986, 
62.4895416142391, NA, 49.4675391135788, 68.0905696296147, 158.258905106886, 
155.234362188971, 78.985506298646, NA, 29.2405983157033, 132.738156112287, 
70.6051687518736, NA, 81.2347110848304, NA, 63.1821689823815, 
126.916403041528, 197.242944372875, 50.2165437445578, 41.8166347676294, 
19.6026530174222, 55.1871586503236, NA, NA, 61.9106042908996, 
304.576697455974, 130.533935807895, 66.8840425794455, 71.3678023808582, 
28.3878410355831, 140.742021789012, 27.3803752222594, NA, 74.4651518277612, 
52.1497058407017, 142.414191266191, 127.564998284227, 125.145010134945, 
43.7321148689866, 38.6331271583135, 58.0819348664361, 72.108198752308, 
365.676154531972, NA, 82.7083511982995, 89.8656094325536, 200.341601969171, 
145.742332439519, 70.8179996632033, 13.5487423124374, 99.2361205675185, 
47.2487727822757, 50.9168192930243, 87.8162270292033, 92.8701074010128, 
154.416260262285, NA, 172.419776634479, 25.965191895825, 62.2538915477445, 
87.4269249274705, NA, 263.180027729561, 125.16467889477, 34.7276328728523, 
180.393944598385, 200.503625100566, 94.4582056400157, 85.9226572707136, 
29.9830207722993, 152.713363317792, NA, 23.6717714416099, 51.392237359917, 
137.589083029979, NA, 62.4020028457266, 110.266770796474, 42.5403262991652, 
NA, 67.8651107705013, NA, 115.619862871237, 153.577384426825, 
62.2101717310168, 163.994495089441, 127.8038460439, 133.5534074188, 
56.642466028164, 7.40552374841297, NA, 60.7431212392032, 31.5623619221465, 
142.580707796808, 153.443176847772, 109.756413179703, 118.889293111244, 
155.382514595011, 51.0374010015819, 64.7631359400003, 55.5408502160295, 
58.8723463285753, NA, 146.583363287034, NA, 234.243422816953, 
172.967793920088, 92.1524076908906, NA, 25.077796329853, 206.574555344071, 
50.4646578407553, 19.6198466002532, NA, 118.213103496736, 35.6035478165478, 
84.43638985228, 154.2079167525, 53.7953440032151, 41.1823963077463, 
60.8136008656232, 107.541300942828, 92.2817794398209, 144.594312188318, 
74.5978401343842, 156.716561582166, 154.12450373892, 98.1951885230804, 
102.970803554585, NA, 202.223130216012, 37.7206825641397, 34.7020112875252, 
NA, NA, 137.327970149541, 77.3732471634386, 115.658734225905, 
38.3776526777709, NA, 49.1745256365287, 39.1636181304966, 384.190446298934, 
NA, NA, 237.966736482099, NA, NA, 55.291400435167, 29.199456517101, 
142.750371848117, 41.7364159576769, NA, 90.260790094256, 123.153838075256, 
NA, 92.1124252933962, 44.299118631839, 21.2598595020048, NA, 
60.441513470385, 65.2726968841611, 143.580221866242, 186.663551483505, 
53.8611162240451, 280.152241794967, 189.567063731444, 178.898116742454, 
38.7845327958083, 18.0874307309723, 70.8829432625132, 56.8016005777847, 
NA, 11.282598761782, 84.2560684536724, 126.752868744507, 89.3811847759758
), Group = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 
1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L), levels = c("Group A", 
"Group B"), class = "factor")), row.names = c(NA, -324L), class = "data.frame")
```
