require(ggplot2)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)
require(vroom)

#df1 <- vroom("../results/metrics-population_size_change.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change"))
#df2 <- vroom("../results/metrics-recombination_rate_change.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change"))
#df3 <- vroom("~/Downloads/metrics_admixture_final.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "time", "rate", "blen_est"))
#df<-df3
df <- rbind(df1, df2)
#dfp <- vroom("../results/phlag.tsv", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change"))
#df<-rbind(dfp, df)
df <- vroom("~/Desktop/phlag/shared.phlag/bird-trees/results/x", col_names = c("TN", "FP", "FN", "TP", "p", "r", "type", "method", "edge", "change"))
di <- vroom("~/Desktop/phlag/shared.phlag/bird-trees/population-information.tsv")
di$LABEL <- sub("_", '', di$LABEL)
df$edge <- sub("_", '', df$edge)
df <- merge(df, di, by.x = "edge", by.y = "LABEL")
df$CU <- as.numeric(df$NGEN)/as.numeric(df$SIZE)/2
df<-df %>% filter(!grepl("phlag_mlt", method))
df<-df %>% filter(!grepl("eap005_ena2", method))
df<-df %>% filter(!grepl("rev00", method))
df<-df %>% filter(!grepl("esp005", method))
df<-df %>% filter(!grepl("phlag_dto_noprior_noprior_free", method))
df$emission <- sub("_.*", '', df$method)
df$emission[df$emission == "phlag"] <- "dtohmm"
unique(df$emission)
df$emission[df$emission == "dtohmm"] <- "Topology order"
df$emission[df$emission == "pbphmm"] <- "Bipartition"
df$emission[df$emission == "mlthmm"] <- "Dominant topolgy"
df$emission[df$emission == "gqshmm"] <- "Gaussian QQS"
df$emission[df$emission == "bcbhmm"] <- "Barycentric QQS bins"
df$emission[df$emission == "cqshmm"] <- "K-means QQS bins"
df$expanded <- grepl("_ex", df$method)
df$ilr <- grepl("_ilr", df$method)
df$transition_prior <- grepl("eap01_ena4", df$method)
df$bestK <- grepl("_sc", df$method) |  grepl("cqshmm", df$method)
df$method <- sub("_ex", '', df$method)
unique(df$method)
df <- df %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP),
  fpr=(FP)/(FP+TN),
  fnr=(FN)/(FN+TP)
)
df$emission_prior <- "No prior"
df$emission_prior[df$method == "phlag_dto_eap01_ena4_updated"] <- "Updated branch lengths"
df$emission_prior[df$method == "phlag_dto_eap01_ena4_fixed"] <- "Fixed branch lengths"
df$emission_prior[df$method == "phlag_dto_eap01_ena4_esp001"] <- "Penalty-based prior"
df$emission_prior[df$method == "phlag_dto_eap01_ena4_fesp001"] <- "Penalty-based fixed prior"

write.table(df, "~/Downloads/metrics-final-popsize_recombrate.tsv", sep='\t', quote = FALSE,  row.names = FALSE)



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












phlag_default <- "phlag_dto_eap01_ena4_esp001"
phlag_default_ex <- "phlag_dto_eap01_ena4_esp001_ex"
phlag_no_prior_ex <- "phlag_dto_noprior_noprior_free_ex"
l <- unique((df %>% filter(method == phlag_default_ex))$edge)
# df <- df %>% filter(edge %in% l)
df$True <- df$TP+df$FP
# df<-df %>% mutate(
#   TN=if_else(True>1000, FP+TN, TN),
#   FN=if_else(True>1000, TP+FN, FN),
#   TP=if_else(True>1000, 0, TP),
#   FP=if_else(True>1000, 0, FP),
#  True=if_else(True>1000, 0, True))
names(df)
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
df$emission <- sub("_.*", "", df$method)
df$emission[df$emission == "phlag"] <- "dtohmm"
df <- df %>% filter(method != "dtohmm_xz_noex")
df$ilr <- grepl("_ilr", df$method)
df$sc <- grepl("_sc", df$method)
unique((df %>% filter(emission=="dtohmm"))$method)
unique(df$method)
unique((df %>% filter(emission=="dtohmm"))$method)
# df <- df %>% filter(!grepl("phlag", method))
df <- df %>% filter(!grepl("phlag", method) | grepl(phlag_default_ex, method) | grepl(phlag_default, method) | grepl(phlag_no_prior_ex, method))
df$transition <- F
df$transition <- grepl("eap01_ena4", df$method)
df$prior <- "free"
df$prior[grepl("esp", df$method)] <- "div_penalty"
unique(df$method)
df$method <- sub("_.*", "", df$method)
df$method[grepl("phlag", df$method)] <- "dtohmm"

df <- df %>% filter(method != "bcbhmm" | sc)
df <- df %>% filter(method != "gqshmm" | !ilr)
df <- df %>% filter(method != "cqshmm" | ilr)
table(df$method)

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
  df$method[df$method == "phlag_dto_eap01_ena4_esp001_ex"] <- "Phlag"
  df$method[df$method == "phlag"] <- "Phlag"
  
  return(df) 
}
df<-renamemethod(df)

df$change[df$change == "10X_up"] <- "10x increase"
df$change[df$change == "10X_down"] <- "10x decrease"
df$change[df$change == "recombination_increase"] <- "10x increase"
df$change[df$change == "recombination_suppression"] <- "Supression"

# df$type[df$type == "population_size_change"] <- "Population size"
# df$type[df$type == "recombination"] <- "Recombination rate"
df$type[df$type == "recombination"] <- "recombination_rate_change"

df2<-merge(df, vroom("../admixture_time.txt"))
df2 %>% # filter(grepl("dtohmm", method)) %>%
 #  mutate(method=sub("_tmptmp2", "", sub(".*hmm_", "", method)) ) %>%
  # filter((grepl("_ex", method) & CU>1) | (!grepl("_ex", method) & CU<1)) %>%
 filter(grepl("_ex", method)) %>%
 filter(divergence!=427476) %>%
 # filter(divergence!=427476) %>%
  filter(divergence!=408896) %>%
  # filter(type=="population_size_change") %>%
  mutate(
    method=sub("_.*", "", method),
    # cu_bin=cut(CU,c(0, 0.1, 1, 2.5, 4)),
    # cu_bin=cut(CU, c(0, 0.1, 1, Inf)),
    p_bin=cut(p, c(0, 0.05, 0.15, 0.25)),
    time_bin=cut(divergence, c(0, 350000, 650000))
  ) %>% #filter(expanded) %>%
  ggplot() +
  facet_grid(cols = vars(p_bin), rows = vars(time_bin)) +
  # facet_wrap(~divergence) +
  aes(x=method, y=f1, fill=method) +
  stat_summary(geom="bar", position = position_dodge2(), color="black") +
  stat_summary(geom="errorbar", color="black", position = position_dodge2()) +
  theme_bw() +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(angle = 90))

th <- 75
df %>%
  # filter(type=="population_size_change") %>%
  # filter(type=="recombination") %>%
  mutate(
    cu_bin=cut(CU, c(0, 0.1, 1, Inf)),
    p_bin=cut(p, c(0, 0.05, 0.15, 0.25))
  ) %>%
  pivot_wider(
    id_cols = c(change, cu_bin, p_bin, edge, type, change, method, p),
    names_from = expanded,
    values_from = c(f1, recall, precision, fpr, fnr, True)
    ) %>% mutate(f1=if_else(True_FALSE <= th, f1_TRUE, f1_FALSE), f1ex = f1_TRUE) %>%
  ggplot() +
  facet_grid(cols = vars(type, change), rows = vars(p_bin)) +
  aes(x=cu_bin, y=f1, fill=method) +
  stat_summary(geom="bar", position = position_dodge2(), color="black") +
  # stat_summary(aes(y=f1ex, color=method), geom="errorbar", position = position_dodge2()) +
  stat_summary(geom="errorbar", color="black", position = position_dodge2()) +
  theme_bw() +
  scale_fill_brewer(palette="Paired") +
  scale_color_brewer(palette="Paired") +
  labs(x="True branch length in CU", y="F1-score")
  
tt <- df %>%
  filter(type=="recombination_rate_change") %>%
  mutate(
    cu_bin=cut(CU, c(0, 0.1, 1, Inf)),
    p_bin=cut(p, c(0, 0.05, 0.15, 0.25))
  ) 
# write.table(tt, "~/Downloads/metrics-population_size_change.tsv", sep='\t', quote = FALSE,  row.names = FALSE)
write.table(dfx, "~/Downloads/metrics-admixture-edited-final.tsv", sep='\t', quote = FALSE,  row.names = FALSE)
