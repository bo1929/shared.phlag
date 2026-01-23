require(ggplot2)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)
require(vroom)

dfpop <- vroom("../population-sizes.tsv")
dfpop %>% filter(SIZE>1) %>%
  ggplot() +
  aes(x=SIZE) +
  stat_density() +
  scale_x_continuous(transform = "log", n.breaks = 10, breaks = c(10000, 100000, 1000000, 10000000)) +
  theme_minimal_grid() +
  labs(x="Population size")

dfgen <- vroom("../num-generations.tsv")
dfgen %>%
  ggplot() +
  aes(x=as.numeric(NGEN)) +
  stat_density() +
  scale_x_continuous(transform = "log", n.breaks = 10, breaks = c(10000, 100000, 300000, 1000000, 10000000)) +
  theme_minimal_grid() +
  labs(x="# of generations")

df_pw <- vroom("../results/pairwise_dist-simulated_vs_stiller.tsv")
names(df_pw)[1] <- "Leaf"
df_pw %>% pivot_longer(cols = -c("Leaf")) %>%
  ggplot() +
  aes(x=value) +
  # stat_ecdf() +
  stat_density() +
  theme_bw() +
  labs(x="Avg. branch length difference (simulated - Stiller2024)/(Stiller2024)")

df_pw <- vroom("../results/pairwise_dist-inferred_vs_stiller.tsv")
names(df_pw)[1] <- "Leaf"
df_pw %>% pivot_longer(cols = -c("Leaf")) %>%
  ggplot() +
  aes(x=value) +
  # stat_ecdf() +
  stat_density() +
  theme_bw() +
  labs(x="Avg. branch length difference (inferred - Stiller2024)/(Stiller2024)")

df_qqs_se <- vroom("../results/qqs_dist-simulated_inferred.tsv")
df_qqs_se$type <- "Inferred"
df_qqs_si <- vroom("../results/qqs_dist-inferred_stiller.tsv")
df_qqs_si$type <- "Stiller"
df_qqs_st <- vroom("../results/qqs_dist-simulated_true.tsv")
df_qqs_st$type <- "Simulated"
df_qqs <- rbind(df_qqs_se, df_qqs_si, df_qqs_st)
names(df_qqs)[1] <- "Clade"

df_qqs %>% pivot_longer(cols = -c("Clade", "type")) %>%
  ggplot() +
  aes(x=value, fill=type) +
  facet_wrap(~name) +
  stat_density(color="black", position = position_identity(), alpha=0.3) +
  theme_bw() +
  labs(x="QQS distribution")

df_rfi <- vroom("../results/rf_inferred_simulated_to_species_tree.tsv", col_names = c("n", "y", "yn"))
df_rfi$type <- "Inferred"
df_rft <- vroom("../results/rf_true_simulated_to_species_tree.tsv", col_names = c("n", "y", "yn"))
df_rft$type <- "True"
df_rfs <- vroom("../results/rf_inferred_stiller_to_species_tree.tsv", col_names = c("n", "y", "yn"))
df_rfs$type <- "Stiller"
df_rf <- rbind(df_rft, df_rfi, df_rfs)
df_rf %>% ggplot() +
  aes(x=yn, fill=type)+
  stat_density(position = position_identity(), alpha=0.5) +
  theme_bw() +
  labs(x="RF distance to the species tree")

df_r <- vroom("../results/zz", col_names = c("TN", "FP", "FN", "TP", "p", "r", "method", "instance", "ps"))
df_r <- df_r %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP),
  fpr=(FP)/(FP+TN),
  fnr=(FN)/(FN+TP)
)
# df_r$LABEL <- sub(pattern = "anomaly_case-", "", sub(pattern = "_10X_.*", "", df_r$instance))
df_r$LABEL <- sub(pattern = "anomaly_case-recombination_increase-", "", df_r$instance)
df_r <- merge(df_r, vroom("../population-information.tsv"))
# df_r$di <- sub(pattern = "anomaly_case-.*_10X_", "", df_r$instance)
df_r$di <- "recombination"
df_r$CU <- as.numeric(df_r$NGEN)/as.numeric(df_r$SIZE)/2
# df_r <- rbind(df_r, x)

# x <- df_r
df_r %>% filter(!method %in% c("mlt-hmm", "phlag")) %>%
  ggplot() +
  facet_wrap(~p) +
  # stat_ecdf(aes(x=f1, color=method, linetype="fpr")) +
  stat_ecdf(aes(x=fpr, color=method, linetype="fnr")) +
  theme_bw() +
  labs(x="Rate") +
  scale_color_brewer(palette = "Paired")

df_r %>% pivot_wider(id_cols =c("di", "p", "r", "instance", "ps", "CU"), names_from = c("method"), values_from = c("f1")) %>%
  ggplot() +
  facet_grid(cols = vars(di), rows=vars(cut(CU, 4))) +
  aes(x=p, y=`pbp-hmm`-`dto-hmm`)+
  stat_summary() +
  theme_bw()

df_r %>% # filter(method %in% c("bcb-hmm", "phlag_xq001", "phlag_xq005","phlag_xq", "dto-hmm", "pbp-hmm", "phlag_xq_bothfree", "dtohmm_en")) %>%
  ggplot() +
  # facet_wrap(~p) +
  stat_ecdf(aes(x=fnr, color=method, linetype="f1")) +
  theme_bw() +
  scale_color_brewer(palette = "Paired")

l <- unique((df_r %>% filter(method == "phlag_xq"))$LABEL)

df_r$expand <- grepl("_ex", x = df_r$method)
df_r$method <- sub("_ex", "", df_r$method)

df_r %>% filter(LABEL %in% l & p == 0.15) %>% mutate(CUb = cut(CU, c(0, 1, Inf))) %>%
  group_by(CUb, method, expand, di) %>%
  summarize(fpr=mean(fpr), fnr=mean(fnr), f1=mean(f1), recall=mean(recall), precision=mean(precision)) %>%
  pivot_wider(id_cols=c("CUb", "method", "di"), names_from = "expand", values_from = c("f1", "fpr", "fnr", "precision", "recall")) %>%
  ggplot() +
  facet_grid(cols=vars(di), rows = vars(CUb)) +
  geom_segment(aes(x=method, color=method, y=f1_FALSE, yend=f1_TRUE), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_dodge2()) +
  theme_bw() +
  scale_color_brewer(palette = "Paired")
  

df_r %>% filter(LABEL %in% l) %>% # filter(method %in% c("bcb-hmm", "phlag_xq", "phlag_xc", "dto-hmm", "pbp-hmm")) %>%
  # filter(p == 0.20) %>% # filter(di != "down") %>%
  filter(LABEL != "N36" & LABEL !="N722") %>%
  ggplot() +
  aes(
    x=cut(CU, c(0, 0.66, Inf)),
    # x=p, # cut(p, c(0, 0.10, 0.15, 0.25)),
    y=f1,
    fill=method,
    color=method
  ) +
  # facet_grid(cols = vars(di), rows = vars(cut(CU, c(0, 0.1, 1, 3)))) +
  facet_grid(cols = vars(di), rows = vars(cut_number(p, 4))) +
  # facet_wrap(~p, scale="free") +
  # stat_summary(aes(group=method), geom="line") +
  stat_summary(geom="bar", position = position_dodge2()) +
  theme_bw() +
  labs(x="Anomaly portion", y="F1-score") +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired")#  +
scale_x_continuous(n.breaks = 6)

# ADD ALSO THE PARENT?
df_r %>% # filter(di!="down") %>%
  ggplot() +
  aes(
      x=cut(as.numeric(NGEN), c(10000, 50000, 100000, 500000, 1000000, 5000000)),
      y=recall,
      fill=method
    ) +
  # facet_grid(cols = vars(type), rows = vars(method)) +
  facet_grid(cols = vars(di), rows = vars(p)) +
  geom_boxplot() +
  # stat_summary(geom="bar", position = position_dodge()) +
  # stat_summary(geom="errorbar", position = position_dodge()) +
  theme_bw() +
  labs(x="Number of generations", y="F1-score") +
  theme(axis.text.x = element_text(angle=30, hjust = 1))

df_r %>% # filter(di!="up") %>%
  ggplot() +
  aes(
      # x=cut(as.numeric(NGEN)/as.numeric(SIZE), c(0, 0.1, 0.33, 1, 3, 10)),
      x=cut(as.numeric(SIZE), c(10000, 100000, 500000, 1000000, 10000000)),
      y=f1,
      fill=method
    ) +
  facet_grid(cols = vars(di), rows = vars(p)) +
  geom_boxplot() +
  # stat_summary(geom="bar", position = position_dodge()) +
  # stat_summary(geom="errorbar", position = position_dodge()) +
  theme_bw() +
  labs(x="Population size", y="F1-score") +
  theme(axis.text.x = element_text(angle=30, hjust = 1))

  df_r %>% # filter(di!="up") %>%
  ggplot() +
  aes(
    # x=cut(as.numeric(NGEN)/as.numeric(SIZE), c(0, 0.1, 0.33, 1, 3, 10)),
    x=cut_number(CU, 5),
    # x=cut(as.numeric(NGEN)/as.numeric(SIZE)*2, c(0, 0.5, 1.0, 4.0, Inf)),
    y=f1,
    fill=method
  ) +
  facet_grid(cols = vars(di), rows = vars(p)) +
  stat_summary(geom="bar", position = position_dodge()) +
  # stat_summary(geom="errorbar", position = position_dodge()) +
  #geom_boxplot() +
  theme_bw() +
  labs(x="CU", y="F1-score") +
  theme(axis.text.x = element_text(angle=30, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

df_r %>% # filter(di!="down") %>%
  ggplot() +
  aes(
    x=cut(DEPTH_NGEN, c(10000, 50000, 100000, 1000000, 2500000, 5000000, 7500000, 10000000, 50000000)),
    y=f1,
    fill=method
  ) +
  facet_grid(cols = vars(di), rows = vars(p)) +
  stat_summary(geom="bar", position = position_dodge()) +
  stat_summary(geom="errorbar", position = position_dodge()) +
  theme_bw() +
  labs(x="# of generations from the root", y="F1-score") +
  theme(axis.text.x = element_text(angle=30, hjust = 1))

df_r %>% # filter(di!="up") %>%
  ggplot() +
  aes(
    x=cut_number(as.numeric(HEIGHT_NGEN), 4),
    y=f1,
    fill=method
  ) +
  facet_grid(cols = vars(di), rows = vars(p)) +
  stat_summary(color="black", geom="bar", position = position_dodge()) +
  stat_summary(geom="errorbar", position = position_dodge()) +
  theme_bw() +
  labs(x="Height in number of edges", y="F1-score") +
  theme(axis.text.x = element_text(angle=30, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

df_r %>% filter(method %in% c("pbp-hmm", "phlag", "bcb-hmm", "dtohmm_en")) %>%
  ggplot() +
  aes(
    x=as.numeric(CU),
    y=recall,
    fill=method
  ) +
  facet_grid(cols = vars(di)) +
  # geom_point() +
  stat_summary_bin(geom="bar", position = position_dodge(width = 0.5), bins = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=30, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

dfp <- vroom("../results/phlag_metrics-admixture.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change"))
df <- vroom("../results/metrics-admixture.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change", "CUE"))
df<-rbind(dfp, df)
di <- vroom("../population-information.tsv")
di$LABEL <- sub("_", '', di$LABEL)
df <- merge(df, di, by.x = "edge", by.y = "LABEL")
df$CU <- as.numeric(df$NGEN)/as.numeric(df$SIZE)/2
names(df)

l <- unique((df %>% filter(method == "phlag_dto_eap01_ena4_free_ex"))$edge)
df <- df %>% filter(edge %in% l)
unique(df$method)

df <- df %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP),
  fpr=(FP)/(FP+TN),
  fnr=(FN)/(FN+TP)
)
df <- df %>% filter(!grepl("_hmm", method))
df$expanded <- grepl("_ex", df$method)
df$ilr <- grepl("_ilr", df$method)
df$sc <- grepl("_sc", df$method)
# df <- df %>% filter(!grepl("phlag", method))

df <- df %>% filter(!grepl("phlag", method)|grepl("phlag_dto_eap01_ena4_esp001_ex", method)|grepl("phlag_dto_eap01_ena4_esp001", method))
df$method <- sub("_.*", "", df$method)


renamemethod <- function(df) {
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
  df$method[df$method == "bcbhmm"] <- "BarycentricBins(6)+HMM"
  df$method[df$method == "bcbhmm_nb3"] <- "BarycentricBins(3)+HMM"
  df$method[df$method == "bcbhmm-ns-9"] <- "BarycentricBins(9)+HMM"
  df$method[df$method == "gqshmm"] <- "GaussianQQS+HMM"
  df$method[df$method == "cqshmm"] <- "KmeansQQS+HMM"
  df$method[df$method == "phlag_dto_eap01_ena4_esp001_ex"] <- "zPhlag"
  df$method[df$method == "phlag"] <- "zPhlag"
  
  return(df) 
}

df %>% filter(method %in% c("bcbhmm")) %>%
  ggplot() +
  facet_wrap(change~type) +
  stat_summary(aes(x=method, shape=expanded, color=sc, y=f1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()

df <- df %>% filter(method != "bcbhmm" | sc)

df %>% filter(method %in% c("gqshmm")) %>%
  ggplot() +
  facet_wrap(change~type) +
  stat_summary(aes(x=method, shape=expanded, color=ilr, y=f1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()

df <- df %>% filter(method != "gqshmm" | !ilr)

df %>% filter(method %in% c("cqshmm")) %>%
  ggplot() +
  facet_wrap(change~type) +
  stat_summary(aes(x=cut(CU, c(0, 0.1, 1, Inf)), shape=expanded, color=ilr, y=f1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()

df <- df %>% filter(method != "cqshmm" | ilr)

unique(df$method)
df<-renamemethod(df)

df %>% # filter(p==0.15) %>% 
  mutate(CU_bin = cut(CU, c(0, 1, Inf))) %>%
  group_by(CU_bin, method, change, type, expanded) %>%
  summarize(fpr=mean(fpr), fnr=mean(fnr), f1=mean(f1), recall=mean(recall), precision=mean(precision)) %>%
  pivot_wider(id_cols=c("CU_bin", "method", "change", "type"), names_from = "expanded", values_from = c("f1", "fpr", "fnr", "precision", "recall")) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(CU_bin)) +
  # geom_segment(aes(x=method, color=method, y=f1_FALSE, yend=f1_TRUE), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_dodge2()) +
  # geom_segment(aes(x=method, color=method, y=fpr_FALSE, yend=fpr_TRUE), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_dodge2()) +
  geom_segment(aes(x=fpr_FALSE, xend=fpr_TRUE, color=method, y=fnr_FALSE, yend=fnr_TRUE), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_dodge2()) +
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle=90)) +
  labs(x="FPR change", y="FNR change")

df %>% # filter(p>0.02) %>%
  mutate(CU_bin = cut(CU, c(0, 0.1, 1, Inf))) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(CU_bin)) +
  aes(y=f1, x=expanded, fill=method) +
  # stat_ecdf(aes(x=f1, linetype=expanded, color=method)) +
  stat_summary(geom="bar", color="black", position = position_dodge2()) +
  stat_summary(geom="errorbar", position = position_dodge2()) +
  theme_bw() +
  scale_fill_brewer(palette = "Paired")

df %>% filter(edge != "N521") %>% # filter(p>0.02) %>%
  # filter((expanded & method != "bcbhmm") | (!expanded & method == "bcbhmm")) %>%
  filter(expanded | grepl("phlag", method)) %>%
  mutate(CU_bin = cut(CU, c(0, 0.1, 1, Inf))) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(CU_bin)) +
  aes(y=f1, x=p, fill=method) +
  stat_summary(geom="bar", color="black", position = position_dodge2()) +
  stat_summary(geom="errorbar", position = position_dodge2()) +
  theme_bw() +
  scale_fill_brewer(palette = "Paired")

df$expanded <- grepl("_ex", df$method)
df$derived <- grepl("_derived", df$method)
df$method <- sub("_.*", "", df$method)
df %>% filter(edge != "N251") %>% # filter(p>0.02) %>%
  # filter((expanded & method != "bcbhmm") | (!expanded & method == "bcbhmm")) %>%
  # filter(expanded | grepl("phlag", method)) %>%
  # filter(type == "recombination") %>%
  mutate(CU = cut(as.numeric(CU), c(0, 1, Inf))) %>%
  filter(!derived & !expanded) %>%
  # mutate(CU_bin = cut_number(CU, 4)) %>%
  ggplot() +
  facet_grid(cols=vars(1*pa, 3*derived+1, 2*expanded+1), rows=vars(edge)) + # rows = vars(cut(p, c(0, 0.05, 0.15, 0.25)))
  aes(x=change, y=recall, color=method, fill=method) +
  stat_summary(geom="bar", color="black", position = position_dodge2()) +
  stat_summary(geom="errorbar", color="black", position = position_dodge2()) +
  # stat_ecdf(aes(x=f1, color=method)) +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  labs(x="True coalescent unit", y="F1-score")# +
  theme(legend.position = "bottom", legend.direction = "horizontal")

df %>% filter(edge !="N251") %>%
  # filter((expanded & method != "bcbhmm") | (!expanded & method == "bcbhmm")) %>%
  # filter(expanded | grepl("phlag", method)) %>%
  # filter(type == "recombination") %>%
  filter(expanded & derived) %>%
  mutate(CU_bin = cut_number(as.numeric(NGEN), 5)) %>%
  # mutate(CU_bin = cut_number(CU, 4)) %>%
  ggplot() +
  facet_grid(cols=vars(interaction(pa, change)), rows=vars(cut(p, c(0, 0.05, 0.15, 0.20, 0.25)))) + # rows = vars()
  aes(x=method, y=f1, color=method, fill=method) +
  stat_summary(geom="bar", color="black", position = position_dodge2()) +
  stat_summary(geom="errorbar", position = position_dodge2()) +
  # stat_ecdf(aes(x=f1, color=method)) +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  labs(x="True coalescent unit", y="F1-score")# +
  theme(legend.position = "bottom", legend.direction = "horizontal")


## do no neighbor if len < 0.7
## do no neighbor first, if not finds below p^\prime, then add neighbors
df %>%
  # filter((expanded & method != "bcbhmm") | (!expanded & method == "bcbhmm")) %>%
  filter(!expanded | grepl("phlag", method)) %>%
  filter(type == "population_size_change") %>%
  mutate(CU_bin = cut(CU, c(0, 0.25, 1, 2.5, Inf))) %>%
  # mutate(CU_bin = cut_number(CU, 4)) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(cut(p, c(0, 0.05, 0.15, 0.25)))) +
  aes(y=f1, x=CU_bin, fill=method) +
  stat_summary(geom="bar", color="black", position = position_dodge2()) +
  stat_summary(geom="errorbar", position = position_dodge2()) +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(x="True coalescent unit", y="F1-score") +
  theme(legend.position = "bottom", legend.direction = "horizontal")

dfp <- vroom("../results/phlag-revert.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change"))
dfp <- merge(dfp, vroom("../population-information.tsv"), by.x = "edge", by.y = "LABEL")
dfp$CU <- as.numeric(dfp$NGEN)/as.numeric(dfp$SIZE)/2
names(dfp)

l <- unique((dfp %>% filter(method == "phlag_dto_eap01_ena4_esp001_ex"))$edge)
dfp <- dfp %>% filter(edge %in% l)
unique(dfp$method)

dfp <- dfp %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP),
  fpr=(FP)/(FP+TN),
  fnr=(FN)/(FN+TP)
)
dfp$expanded <- grepl("_ex", dfp$method)

dfp %>%
  # filter((expanded & method != "bcbhmm") | (!expanded & method == "bcbhmm")) %>%
  filter(expanded) %>%
  filter(type == "recombination") %>%
  # filter(type == "population_size_change") %>%
  mutate(CU_bin = cut(CU, c(0, 0.25, 1, Inf))) %>%
  # mutate(CU_bin = cut_number(CU, 4)) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(cut(p, c(0, 0.05, 0.15, 0.25)))) +
  aes(y=f1, x=CU_bin, fill=method, color=method) +
  stat_summary(geom="bar", color="black", position = position_dodge2()) +
  stat_summary(geom="errorbar", position = position_dodge2()) +
  # stat_summary(geom="point") +
  # stat_summary(geom="line", aes(group=method)) +
  # stat_ecdf(aes(x=fpr, color=method)) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  theme_bw() +
  labs(x="True coalescent unit", y="F1-score")

df %>%
  # filter((expanded & method != "bcbhmm") | (!expanded & method == "bcbhmm")) %>%
  # filter(expanded) %>%
  # filter(type == "recombination") %>%
  # filter(type == "population_size_change") %>%
  group_by(type, change, method, edge) %>% summarize(f1=mean(f1))
dft<-vroom("../results/admixture_timepast.tsv", col_names = c("edge", "a1", "a2", "timepast", "pa"))
df<-merge(dft, df, by="edge")
library(ggcorrplot)
ggcorrplot(cor(
  df %>%
    select(f1, CU, NGEN, timepast, pa, X13, p, change) %>%
    mutate(change=(change=="025")*1, NGEN=as.numeric(NGEN))
  ), method = "circle")
names(df)
cor(
  df %>%
    select(f1, CU, NGEN, timepast, pa, X13, p) %>%
    mutate(NGEN=as.numeric(NGEN))
)

dfp %>% filter(expanded) %>% filter(method %in% c(
  "phlag_dto_eap01_ena4_free_ex",
  # "phlag_dto_eap01_ena4_rev0001_ex",
  # "phlag_dto_eap005_ena2_esp001_ex",
  "phlag_dto_eap01_ena4_esp001_ex"
  # "phlag_dto_eap005_ena2_esp001_ex",
  # "phlag_mlt_eap01_ena4_free_ex",
  # "phlag_dto_eap01_ena4_esp001_ex"
  # "phlag_mlt_eap01_ena4_esp001_ex"
  )) %>%
  mutate(
    CU_bin = cut(CU, c(0, 0.1, 1, 5)),
    p_bin = cut(p, c(0, 0.05, 0.15, 0.25)),
    free=grepl("free", method),
    method=sub("_eap.*", "", method)) %>%
  group_by(CU_bin, p_bin, free, change, type, method) %>%
  summarise(f1=mean(f1)) %>%
  pivot_wider(id_cols = c("CU_bin", "method", "p_bin", "change", "type"), names_from = "free", values_from = c("f1")) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(p_bin)) +
  geom_segment(aes(x=CU_bin, y=`TRUE`, yend=`FALSE`, color=method), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_dodge2()) +
  theme_bw()

unique(df$method)
df<-renamemethod(df)

dfsub <- df %>% filter(method %in% c("dtohmm", "dtohmm_ex", "phlag_dto_eap01_ena4_free_ex", "phlag_dto_eap01_ena4_esp001_ex"))
dfsub$method[dfsub$method == "dtohmm"] <- "TopologyOrderSB+HMM"
dfsub$method[dfsub$method == "dtohmm_ex"] <- "TopologyOrder+HMM"
dfsub$method[dfsub$method == "phlag_dto_eap01_ena4_free_ex"] <- "TopologyOrder+HMM+TransitionPriors"
dfsub$method[dfsub$method == "phlag_dto_eap01_ena4_esp001_ex"] <- "TopologyOrder+HMM+TransitionPriors+ST-EM"

dfsub %>%  mutate(
  CU_bin = cut(CU, c(0, 0.1, 1, Inf)),
  p_bin = cut(p, c(0, 0.05, 0.15, 0.25)),
  trp=grepl("TransitionPriors", method),
  em=grepl("ST-EM", method),
  sb=grepl("SB", method)
) %>%
  group_by(CU_bin, p_bin, trp, sb, em, change, type, method) %>%
  summarise(f1=mean(f1)) %>%
  pivot_wider(
    id_cols = c("CU_bin", "p_bin", "change", "type"),
    names_from = c("sb", "trp",  "em"),
    values_from = c("f1")
    ) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(p_bin)) +
  geom_segment(aes(x=CU_bin, y=`TRUE_FALSE_FALSE`, yend=`FALSE_FALSE_FALSE`, color="Neighbors"), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_jitter(width = 0.1)) +
  geom_segment(aes(x=CU_bin, y=`FALSE_FALSE_FALSE`, yend=`FALSE_TRUE_FALSE`, color="Neighbors+TransitionPriors"), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_jitter(width = 0.1)) +
  geom_segment(aes(x=CU_bin, y=`FALSE_TRUE_FALSE`, yend=`FALSE_TRUE_TRUE`, color="Neighbors+TransitionPriors+ST-EM"), linewidth = 1, arrow = grid::arrow(length = unit(5, units="pt")), position = position_jitter(width = 0.1)) +
  theme_bw() +
  labs(x="CU bins", y="F1-score", color="Added component") +
  theme(legend.position = "bottom", legend.direction = "horizontal")




dfp <- vroom("../results/phylter_metrics.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change"))
dfp <- merge(dfp, vroom("../population-information.tsv"), by.x = "edge", by.y = "LABEL")
dfp$CU <- as.numeric(dfp$NGEN)/as.numeric(dfp$SIZE)/2
names(dfp)

unique(dfp$method)

dfp <- dfp %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP),
  fpr=(FP)/(FP+TN),
  fnr=(FN)/(FN+TP)
)

dfp %>%
  # filter((expanded & method != "bcbhmm") | (!expanded & method == "bcbhmm")) %>%
  # filter(type == "recombination") %>%
  # filter(type == "population_size_change") %>%
  mutate(CU_bin = cut(CU, c(0, 0.25, 1, Inf))) %>%
  # mutate(CU_bin = cut_number(CU, 4)) %>%
  ggplot() +
  facet_grid(cols=vars(change, type), rows = vars(cut(p, c(0, 0.05, 0.15, 0.25)))) +
  aes(y=precision, x=CU_bin, fill=method, color=method) +
  stat_summary(geom="bar", color="black", position = position_dodge2()) +
  stat_summary(geom="errorbar", color="black", position = position_dodge2()) +
  # stat_summary(geom="point") +
  # stat_summary(geom="line", aes(group=method)) +
  # stat_ecdf(aes(x=f1, color=method)) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  theme_bw() +
  labs(x="True coalescent unit", y="F1-score")

