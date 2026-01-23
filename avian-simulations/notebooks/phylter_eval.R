require(ggplot2)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)
require(vroom)


df <- vroom("../results/phylter_metrics.tsv", col_names = c("TN", "FP", "FN", "TP",  "p", "r", "type", "method", "edge","change"))
df <- df %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP),
  fpr=(FP)/(FP+TN),
  fnr=(FN)/(FN+TP),
  tpr=(TP / (TP + FN))
)
di <- vroom("../population-information.tsv")
di$LABEL <- sub("_", '', di$LABEL)
df$edge <- sub("_", '', df$edge)
df <- merge(df, di, by.x = "edge", by.y = "LABEL")
df$CU <- as.numeric(df$NGEN)/as.numeric(df$SIZE)/2

df$method[df$method == "phylter_k_300"] <- "k=3.00"
df$method[df$method == "phylter_k_050"] <- "k=0.50"
df$method[df$method == "phylter_k_155"] <- "k=1.55"
df$method[df$method == "phylter_k_020"] <- "k=0.20"

df$change[df$change == "10X_up"] <- "10x increase"
df$change[df$change == "10X_down"] <- "10x decrease"
df$change[df$change == "recombination_increase"] <- "10x increase"
df$change[df$change == "recombination_suppression"] <- "Suppression"

df$type[df$type == "population_size_change"] <- "Population size"
df$type[df$type == "recombination"] <- "Recombination rate"


df %>%
  mutate(p=cut(p, c(0, 0.05, 0.15, 0.25)), ) %>%
  group_by(method, p, change, type) %>%
  summarise(fpr=mean(fpr), tpr=mean(tpr)) %>%
  ggplot() +
  facet_grid(cols = vars(type, change), rows = vars(p)) +
  aes(x=fpr, y=tpr, color=method) +
  geom_point(size=3) +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  labs(x="FPR", y="TPR", color="")
ggsave2("../figures/phylter_k.pdf", width = 8.5, height = 5)
