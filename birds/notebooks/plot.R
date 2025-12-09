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
