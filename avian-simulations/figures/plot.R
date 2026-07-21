library(tidyverse)

# ---------------------------------------------------------------------------
# TNR-TPR figure (Phlag emission / prior comparison)
# ---------------------------------------------------------------------------
p_data <- read.csv("./metrics-final-popsize_recombrate.tsv", sep = "\t")

phlag_data <- p_data %>%
  filter(!emission %in% c("Barycentric QQS bins", "K-means QQS bins") | bestK == T) %>%
  filter(emission != "Bipartition") %>%
  filter((expanded == F & CU < 1) | (expanded == T & CU >= 1)) %>%
  mutate(
    cu_bin = cut(CU, c(0, 1, Inf)),
    p_bin = cut(p, c(0, 0.05, 0.15, 0.25)),
    TNR = TN / (TN + FP),
    TPR = TP / (TP + FN)
  ) %>%
  group_by(p_bin, method, type, emission_prior, emission, transition_prior, ilr) %>%
  summarise(TPR = mean(TPR), TNR = mean(TNR))

phlag_data$transition_prior <- factor(phlag_data$transition_prior, levels = c(TRUE, FALSE))
phlag_data$emission <- factor(
  phlag_data$emission,
  levels = c("Topology order", "Bipartition", "Barycentric QQS bins", "K-means QQS bins", "Dominant topolgy", "Gaussian QQS"),
  labels = c("Topology order", "Bipartition", "Barycentric bins", "K-means bins", "Dominant topolgy", "Gaussian")
)
phlag_data$emission_prior <- factor(
  phlag_data$emission_prior,
  levels = c("Penalty-based prior", "Penalty-based fixed prior", "Updated branch lengths", "Fixed branch lengths", "No prior"),
  labels = c("prior-updated", "prior-fixed", "MSC-updated", "MSC-fixed", "segment")
)

ggplot(data = phlag_data, aes(x = 1 - TNR, y = TPR, color = emission_prior,
                              group = interaction(emission, emission_prior, transition_prior, ilr))) +
  geom_path(aes(linetype = transition_prior)) +
  geom_point(aes(size = p_bin, shape = interaction(ilr, emission))) +
  facet_wrap(
    ~type,
    labeller = labeller(type = c(
      "population_size_change" = "population size change",
      "recombination_rate_change" = "recombination rate change"
    )),
    ncol = 1
  ) +
  scale_size_manual(values = c(1, 2, 3)) +
  scale_linetype(labels = c("Equation (1)", "None")) +
  scale_color_manual(values = c("#1f78b4", "#a6cee3", "#e31a1c", "#fb9a99", "#8073ac")) +
  scale_shape_manual(
    values = c(19, 3, 2, 17, 8, 0, 15),
    labels = c("Topology order", "Barycentric bins", "K-means bins", "K-means bins (ILR)",
               "Dominant topolgy", "Gaussian", "Gaussian (ILR)")
  ) +
  labs(x = "false positive rate (FPR)", y = "recall (TPR)",
       size = "outlier proportion", shape = "emission type",
       color = "emission parameterization", linetype = "transition prior") +
  theme_bw() +
  guides(
    color = guide_legend(order = 1, ncol = 1),
    shape = guide_legend(order = 2, ncol = 1),
    size = guide_legend(order = 4, ncol = 1),
    linetype = guide_legend(order = 3, ncol = 1)
  ) +
  theme(legend.position = "right", legend.box = "vertical")

# ---------------------------------------------------------------------------
# Compare methods: F1 by branch length (pop size / recomb)
# ---------------------------------------------------------------------------
MA_data <- read.csv("./metrics-MA-all.tsv", header = F, sep = " ")
names(MA_data) <- c("type", "dir", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1")

MA_data <- MA_data %>%
  filter(type != "admixture") %>%
  mutate(
    p = str_extract(dir, "p\\d{3}") %>% str_remove("p") %>% as.numeric() / 100,
    change = str_extract(dir, "\\d+X_.*$"),
    method = "MA_ws100"
  )
MA_data[MA_data$type == "population_size_change" & MA_data$change == "10X_up", ]$type <- "population_size_increase"
MA_data[MA_data$type == "population_size_change" & MA_data$change == "10X_down", ]$type <- "population_size_decrease"
MA_data <- MA_data[, c("p", "type", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1", "method")]

phylter_data <- read.csv("./metrics-phylter-popsize_recombrate.tsv", sep = "\t", header = F)
names(phylter_data) <- c("TN", "FP", "FN", "TP", "p", "X6", "type", "phylter", "edge", "change")
phylter_data <- phylter_data %>%
  mutate(
    k = str_extract(phylter, "k_\\d{3}") %>% str_remove("k_") %>% as.numeric() / 100,
    method = "phylter",
    acc = (TP + TN) / (TP + TN + FP + FN),
    f1 = 2 * TP / (2 * TP + FP + FN)
  )
phylter_data[phylter_data$type == "population_size_change" & phylter_data$change == "10X_up", ]$type <- "population_size_increase"
phylter_data[phylter_data$type == "population_size_change" & phylter_data$change == "10X_down", ]$type <- "population_size_decrease"
phylter_data[phylter_data$type == "recombination", ]$type <- phylter_data[phylter_data$type == "recombination", ]$change
phylter_data <- phylter_data[, c("p", "type", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1", "method")]

phlag_data <- p_data %>%
  filter((emission == "Topology order" & transition_prior == T & emission_prior == "Penalty-based prior") |
           emission == "Bipartition") %>%
  filter((expanded == F & CU < 1) | (expanded == T & CU >= 1)) %>%
  mutate(k = 1)
phlag_data[phlag_data$type == "population_size_change" & phlag_data$change == "10X_up", ]$type <- "population_size_increase"
phlag_data[phlag_data$type == "population_size_change" & phlag_data$change == "10X_down", ]$type <- "population_size_decrease"
phlag_data[phlag_data$type == "recombination", ]$type <- phlag_data[phlag_data$type == "recombination", ]$change
phlag_data <- phlag_data[, c("p", "type", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1", "method")]

cu_data <- p_data %>% group_by(edge) %>% summarise(CU = mean(CU))

all_data <- merge(rbind(MA_data %>% filter(k == 0.3), phylter_data %>% filter(k == 0.5), phlag_data), cu_data) %>%
  mutate(
    cu_bin = cut(CU, c(0, 0.1, 1, Inf)),
    p_bin = cut(p, c(0, 0.1, 0.25))
  )
all_data$method <- factor(
  all_data$method,
  levels = c("phylter", "MA_ws100", "pbphmm", "phlag_dto_eap01_ena4_esp001"),
  labels = c("PhylteR", "moving average", "Bipartition+HMM", "Phlag")
)
all_data$type <- factor(
  all_data$type,
  levels = c("recombination_increase", "recombination_suppression", "population_size_increase", "population_size_decrease", "admixture"),
  labels = c("recombination increase", "recombination suppression", "population size increase", "population size decrease", "admixture")
)

all_data %>%
  ggplot(aes(x = as.factor(cu_bin), y = f1, fill = method, group = interaction(method, p_bin))) +
  stat_summary(geom = "bar", position = position_dodge(width = 0.8), width = 0.8, color = "grey20", linewidth = 0.5) +
  stat_summary(geom = "errorbar", position = position_dodge(width = 0.8), width = 0.5, color = "grey20") +
  facet_grid(p_bin ~ type) +
  scale_fill_manual(values = c("#abd9e9", "#b2abd2", "#fdae61", "#66c2a5")) +
  labs(x = "branch length (CU)", y = "f1-score", fill = "") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  coord_cartesian(ylim = c(0, 1))

# ---------------------------------------------------------------------------
# ROC figure
# ---------------------------------------------------------------------------
all_data <- merge(rbind(MA_data, phylter_data, phlag_data), cu_data) %>%
  mutate(
    cu_bin = cut(CU, c(0, 0.1, 1, Inf)),
    p_bin = cut(p, c(0, 0.1, 0.25)),
    TNR = TN / (TN + FP),
    TPR = TP / (TP + FN)
  )
all_data$method <- factor(
  all_data$method,
  levels = c("phlag_dto_eap01_ena4_esp001", "pbphmm", "MA_ws100", "phylter"),
  labels = c("Phlag", "Bipartition+HMM", "moving average", "PhylteR")
)
all_data$type <- factor(
  all_data$type,
  levels = c("population_size_increase", "population_size_decrease", "recombination_increase", "recombination_suppression", "admixture"),
  labels = c("population size change", "population size change", "recombination", "recombination", "admixture")
)

all_data %>%
  group_by(method, k, type, p_bin) %>%
  summarise(TPR = mean(TPR), TNR = mean(TNR)) %>%
  ggplot(aes(x = TNR, y = TPR, color = method, group = interaction(method))) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.5) +
  facet_wrap(p_bin ~ type) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "specificity (TNR)", y = "sensitivity (TPR)", fill = "") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, vjust = 0.6))

# ---------------------------------------------------------------------------
# Admixture: F1 by divergence time
# ---------------------------------------------------------------------------
MA_data <- read.csv("./metrics-MA-all.tsv", header = F, sep = " ")
names(MA_data) <- c("type", "dir", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1")

MA_data <- MA_data %>%
  filter(type == "admixture") %>%
  mutate(
    p = str_extract(dir, "p\\d{3}") %>% str_remove("p") %>% as.numeric() / 100,
    method = "MA_ws100"
  )
MA_data <- MA_data[, c("p", "type", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1", "method")]

phylter_data <- read.csv("./phylter-admixture_metrics.tsv", sep = "\t", header = T) %>%
  mutate(method = "phylter", k = 0.5)
phylter_data <- phylter_data[, c("p", "type", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1", "method")]

ad_data <- read.csv("./metrics-admixture-edited-final.tsv", sep = "\t")
div_data <- ad_data %>% group_by(edge) %>% summarise(divergence = mean(divergence))

phlag_data <- ad_data %>%
  filter(method %in% c("phlag_dto_eap01_ena4_esp001", "pbphmm")) %>%
  filter(divergence != 427476, divergence != 408896) %>%
  mutate(k = 1)
phlag_data <- phlag_data[, c("p", "type", "edge", "k", "TP", "TN", "FP", "FN", "acc", "f1", "method")]

all_data_ad <- merge(rbind(MA_data %>% filter(k == 0.3), phylter_data, phlag_data), div_data) %>%
  mutate(
    div_bin = cut(divergence, c(0, 350000, 650000)),
    p_bin = cut(p, c(0, 0.1, 0.25))
  )
all_data_ad$method <- factor(
  all_data_ad$method,
  levels = c("phylter", "MA_ws100", "pbphmm", "phlag_dto_eap01_ena4_esp001"),
  labels = c("PhylteR", "moving average", "Bipartition+HMM", "Phlag")
)

all_data_ad %>%
  ggplot(aes(x = as.factor(div_bin), y = f1, fill = method, group = interaction(method, p_bin))) +
  stat_summary(geom = "bar", position = position_dodge(width = 0.8), width = 0.8, color = "grey20", linewidth = 0.5) +
  stat_summary(geom = "errorbar", position = position_dodge(width = 0.8), width = 0.5, color = "grey20") +
  facet_grid(p_bin ~ type) +
  scale_fill_manual(values = c("#abd9e9", "#b2abd2", "#fdae61", "#66c2a5")) +
  scale_x_discrete(labels = c("[0,3.5]", "(3.5,6.5]"), name = expression(divergence ~ time ~ (10^6 ~ years))) +
  labs(y = "f1-score", fill = "") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, vjust = 0.6)) +
  coord_cartesian(ylim = c(0, 1))
