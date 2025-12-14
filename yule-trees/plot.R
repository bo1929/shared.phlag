require(ggplot2)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)
require(vroom)

rename <- function(df) {
  df$method[df$method == "cbphmm"] <- "ConsistentBP+HMM"
  df$method[df$method == "pbphmm"] <- "PresentBP+HMM"
  df$method[df$method == "dtohmm"] <- "TopologyOrder+HMM"
  df$method[df$method == "mlthmm"] <- "BestTopology+HMM"
  df$method[df$method == "bcbhmm_nbin3"] <- "BarycentricBins(3)+HMM"
  df$method[df$method == "bcbhmm_nbin6"] <- "BarycentricBins(6)+HMM"
  df$method[df$method == "bcbhmm_nbin4"] <- "BarycentricBins(4)+HMM"
  df$method[df$method == "bcbhmm-ns-3"] <- "BarycentricBins(3)+HMM"
  df$method[df$method == "bcbhmm-ns-2"] <- "BarycentricBins(2)+HMM"
  df$method[df$method == "bcbhmm-ns-6"] <- "BarycentricBins(6)+HMM"
  df$method[df$method == "bcbhmm_nb6"] <- "BarycentricBins(6)+HMM"
  df$method[df$method == "bcbhmm_nb3"] <- "BarycentricBins(3)+HMM"
  df$method[df$method == "bcbhmm-ns-9"] <- "BarycentricBins(9)+HMM"
  df$method[df$method == "gqshmm"] <- "GaussianQQS+HMM"
  return(df) 
}

df_intdiff <- vroom("results/baseline_metrics-introgression-difficult.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "b", "D", "R", "ST", "method"))
head(df_intdiff)
di <- vroom("simulated_events/introgression/selected-difficult/all_selected.tsv", col_names = c("ST", "R", "D", "ul", "di" ))
df <- merge(df_intdiff, di, by = c("ST", "D", "R"))
df$P <- df$FP + df$TP
df <- df %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP)
)
df %>% group_by(method) %>% summarize(wp=sum(!is.na(f1))/n(), f1=mean(f1, na.rm = TRUE), recall=mean(recall, na.rm = TRUE))
df_s <- df %>% filter(!grepl("phylter", method))
df_s <- rename(df_s) %>% filter(method != "BarycentricBins(4)+HMM")

df_s %>% ggplot() +
  facet_wrap(~p) +
  stat_ecdf(aes(x = precision, color = method)) +
  theme_minimal_grid() +
  scale_color_brewer(palette = "Set1") +
  labs( y = "ECDF", color = "")

df_s %>% 
  mutate(
    ul_bin=cut(ul, c(0, 10000000, 100000000, 1000000000, Inf)),
    di_bin=cut(di, c(0, 10000, 25000, 50000, 75000, 100000, Inf))) %>%
  # group_by(method, di_bin) %>%
  # summarise(TP=sum(TP), FP=sum(FP), TN=sum(FN), FN=sum(FN),  recall = (TP) / (TP + FN), precision = (TP) / (TP + FP), f1 = 2 * precision * recall / (precision + recall)) %>%
  ggplot() +
  facet_wrap(~di_bin)+
  aes(x=method, y=f1, fill=method) +
  # geom_point() +
  # geom_line(aes(group = method)) +
  stat_summary(color="black", geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  theme_minimal_grid() +
  scale_fill_brewer(palette = "Set1") +
  # labs(x="Deviation from ultrametricity") +
  labs(x="Divergence time") +
  theme(axis.text.x = element_blank())

df_s %>%
  ggplot() +
  aes(x=p, y=f1, fill=method) +
  # facet_wrap(~cut_interval(di, 2)) +
  stat_summary(color="black", geom="bar", position = position_dodge2()) +
  theme_minimal_grid() +
  scale_fill_brewer(palette = "Set1")

df_s %>% ggplot() +
  aes(x=di, y=recall, color=method) +
  geom_point(alpha=0.1) +
  theme_minimal_grid() +
  scale_color_brewer(palette = "Set1") +
  # stat_smooth() +
  stat_summary_bin(breaks = c(1, 10000, 50000, 100000, 1000000)) +
  stat_summary_bin(breaks = c(1, 10000,  50000, 100000, 1000000), geom="line") +
  scale_x_continuous(transform = "log", breaks = c(10, 1000, 10000, 100000)) +
  coord_cartesian(xlim = c(1, 1000000)) +
  labs(x="# of generations until speciation")

df_s %>% ggplot() +
  aes(x=ul, y=f1, color=method) +
  geom_point(alpha=0.1) +
  theme_minimal_grid() +
  scale_color_brewer(palette = "Set1") +
  # stat_smooth() +
  stat_summary_bin(breaks = 100000*c(1, 10000, 50000, 100000, 1000000)) +
  stat_summary_bin(breaks = 100000*c(1, 10000,  50000, 100000, 1000000), geom="line") +
  scale_x_continuous(transform = "log", breaks = 10000*c(10, 1000, 10000, 100000)) +
  # coord_cartesian(xlim = c(1, 1000000)) +
  labs(x="Deviation from ultrametricity")


df_s %>% 
  mutate(ul_bin=cut_number(ul, 3), di_bin=cut(di, c(0, 10000, 25000, 50000, 75000, 100000, 150000, 250000, 600000))) %>%
  group_by(method, di_bin) %>%
  summarise(TP=sum(TP), FP=sum(FP), TN=sum(FN), FN=sum(FN),  recall = (TP) / (TP + FN), precision = (TP) / (TP + FP), f1 = 2 * precision * recall / (precision + recall)) %>%
  ggplot() +
  aes(x=di_bin, y=f1, color=method) +
  geom_point() +
  geom_line(aes(group = method)) +
  theme_minimal_grid() +
  scale_color_brewer(palette = "Set1") +
  labs(x="# of generations until speciation") +
  theme(axis.text.x = element_text(angle = 30))


df_s %>% 
  group_by(method) %>%
  summarise(wp = sum( (P>0) /n())) %>%
  ggplot() +
  aes(x=method, y=wp, fill=method) +
  geom_col() +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1") +
  labs(y="Any prediction")

df <- vroom("results/baseline_metrics-population_increase_10X.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "method"))
#df <- vroom("results/population_increase_10X/all_metrics-population_increase_10X.tsv", col_names = c("method", "TN", "FP", "FN", "TP", "p", "r"))
df$P <- df$FP + df$TP
# th <- 1000
# df <- df %>% mutate(FN = ifelse(P > th, TP+FN, FN), TN = ifelse(P > th, TN+FP, TN), TP=ifelse(P >th, 0, TP), FP=ifelse(P > th, 0, FP))

df <- df %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP)
)

df %>% group_by(method) %>% summarize(wp=sum(!is.na(f1))/n(), f1=mean(f1, na.rm = TRUE), recall=mean(recall, na.rm = TRUE))
df_s <- df # %>% filter(method %in% c("phlag_0010-c_min25max81", "phylter-k155-island_0", "BP-HMM"))
df_s$method[df_s$method == "phlag_0010-c_min9max81"] <- "phlag_0010"
df_s <- rename(df_s)
df_s %>% filter(r >0.9) %>%
  ggplot() +
  aes(x=method, y=f1, fill=method) +
  facet_wrap(~(r*p)) +
  stat_summary(geom="bar", color="black") +
  theme_minimal_grid() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_blank())
 

df_s %>% ggplot() +
  # facet_wrap(~p) +
  stat_ecdf(aes(x = f1, color = method)) +
  theme_minimal_grid() +
  scale_color_brewer(palette = "Set1") +
  labs( y = "ECDF", color = "")

df_s %>% 
  group_by(method) %>%
  summarise(wp = sum( (P>0) /n())) %>%
  ggplot() +
  aes(x=method, y=wp, fill=method) +
  geom_col() +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1") +
  labs(y="Any prediction")

df <- vroom("results/baseline_metrics-recombination_suppression_support-default.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "method"))
df$P <- df$FP + df$TP

df <- df %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN))
)

df %>% group_by(method) %>% summarize(wp=sum(!is.na(f1))/n(), f1=mean(f1, na.rm = TRUE))
df_s <- rename(df)

df_s %>% ggplot() +
  # facet_wrap(~p) +
  # facet_wrap(p~(round(r, 2))) +
  # facet_grid(p~(round(r, 2))) +
  facet_wrap(~(round(1-r, 2))) +
  stat_ecdf(aes(x = f1, color = method)) +
  theme_minimal_grid() +
  scale_color_brewer(palette = "Set1") +
  labs( y = "ECDF", color = "") +
  coord_cartesian(xlim = c(0.5, 1))

df_s %>%
  ggplot() +
  aes(x=p, y=f1, color=method) +
  facet_wrap(~(round(1-r, 2))) +
  geom_point(alpha=0.15) +
  stat_summary() +
  stat_summary(geom="line") +
  theme_minimal_grid() +
  scale_color_brewer(palette = "Set1")

df_s %>% 
  group_by(method) %>%
  summarise(wp = sum((P>0) /n())) %>%
  ggplot() +
  aes(x=method, y=wp, fill=method) +
  geom_col() +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1") +
  labs(y="Any prediction")

df_s %>%
  ggplot() +
  aes(x=method, y=f1, fill=method) +
  facet_wrap(~(round(1-r, 2))) +
  stat_summary(color="black", geom="bar") +
  theme_minimal_grid() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_blank())

ds <- vroom("../../phlag/notebooks/branch_summary_stats-recombination_suppression_default.tsv")
df <- vroom("../../phlag/notebooks/baselines_test-all_branches-recombination_suppression_default.tsv")
df <- merge(ds, df, by = "lbl")

df<-rename(df)
df %>% filter(method %in% c("GaussianQQS+HMM", "PresentBP+HMM", "TopologyOrder+HMM", "BestTopology+HMM", "BarycentricBins(6)+HMM")) %>%
  mutate(b=as.numeric(sub("I", "", lbl))) %>%
  ggplot() +
  facet_wrap(~method) +
  aes(x=cut_number(entropy_bcb, 8), y=f1, color=cut_number(blen, 4)) + # , shape=cut_number(blen, 4)
  stat_summary() +
  scale_color_brewer(palette="Paired") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x="Entropy", color="Branch length (CU)")

df %>% filter(method %in% c("GaussianQQS+HMM", "PresentBP+HMM", "TopologyOrder+HMM", "BestTopology+HMM", "BarycentricBins(6)+HMM")) %>%
  ggplot() +
  aes(x=f1, color=method) +
  stat_ecdf() +
  theme_cowplot() +
  scale_color_brewer(palette="Set2")


ds <- vroom("../../phlag/notebooks/branch_summary_stats-population_increase_10X.tsv")
df <- vroom("../../phlag/notebooks/baselines_test-all_branches-population_increase_10X.tsv")
df <- merge(ds, df, by = "lbl")

df<-rename(df)
df %>% filter(method %in% c("GaussianQQS+HMM", "PresentBP+HMM", "TopologyOrder+HMM", "BestTopology+HMM", "BarycentricBins(6)+HMM")) %>%
  mutate(b=as.numeric(sub("I", "", lbl))) %>%
  ggplot() +
  facet_wrap(~method) +
  aes(x=cut_number(entropy_bcb, 8), y=f1, color=cut_number(blen, 4)) + # , shape=cut_number(blen, 4)
  stat_summary() +
  scale_color_brewer(palette="Paired") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x="Entropy", color="Branch length (CU)")

df %>% filter(method %in% c("GaussianQQS+HMM", "PresentBP+HMM", "TopologyOrder+HMM", "BestTopology+HMM", "BarycentricBins(6)+HMM")) %>%
  ggplot() +
  aes(x=f1, color=method) +
  stat_ecdf() +
  theme_cowplot() +
  scale_color_brewer(palette="Set2")
