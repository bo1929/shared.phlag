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
dfgen %>% #filter(NGEN>1 & is.numeric(NGEN)) %>%
  ggplot() +
  aes(x=NGEN) +
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
df_rfs <- vroom("../results/rf_stiller_to_species_tree.tsv", col_names = c("n", "y", "yn"))
df_rfs$type <- "Stiller"
df_rf <- rbind(df_rft, df_rfi, df_rfs)
df_rf %>% ggplot() +
  aes(x=yn, fill=type)+
  stat_density(position = position_identity(), alpha=0.5) +
  theme_bw() +
  labs(x="RF distance to the species tree")

df_r <- vroom("../results/baselines-simulations/results.tsv", col_names = c("method", "type", "instance", "TN", "FP", "FN", "TP", "p", "r"))
df_r <- df_r %>% mutate(
  recall = (TP) / (TP + FN),
  precision = (TP) / (TP + FP),
  f1 = TP/ (TP+0.5*(FP+FN)),
  acc=(TP+TN)/(TP+TN+FN+FP)
)
df_r %>%
  ggplot() +
  aes(y=precision, x = p, color=instance) +
  facet_grid(cols = vars(type), rows = vars(method)) +
  geom_point() +
  geom_line() +
  theme_bw()


