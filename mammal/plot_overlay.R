require(ggplot2)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)
require(vroom)
require(tidyquant)
require(ggnewscale)
library(RColorBrewer)

format_mb <- function(l) {
  l <- paste0(round(l / 1e6), " Mb")
  return(l)
}

taxon <- vroom(
  "taxon_map_order.txt",
  comment = "#",
  delim = "\t",
  col_names = c("name", "taxon")
)
ttaxon <- vroom(
  "taxon_map.txt",
  comment = "#",
  delim = "\t",
  col_names = c("name", "staxon")
)

pos <- vroom("pos", col_names = c("pos"), delim = "\t")
distances <- vroom("./distances_chr3.txt", col_names = c("name", "distance"))

process_pred <- function(raw, pos_df) {
  raw$i <- seq_len(nrow(raw))
  raw$pos <- pos_df$pos
  raw$chr <- 3

  long <- raw %>%
    pivot_longer(cols = -c("i", "pos", "chr")) %>%
    group_by(name) %>%
    mutate(
      p = mean(value, na.rm = TRUE),
      value = ifelse(p > 0.5, 1 - value, value)
    ) %>%
    ungroup()

  long %>%
    arrange(name, i) %>%
    group_by(name) %>%
    mutate(
      value = (function(v) {
        rle_v <- rle(is.na(v))
        fill_runs <- rle_v$values & rle_v$lengths < 10
        fill_idx <- rep(fill_runs, rle_v$lengths)
        v[fill_idx] <- zoo::na.locf(v, na.rm = FALSE)[fill_idx]
        v
      })(value),
      new_seg = is.na(value) |
        is.na(lag(value)) |
        value != lag(value, default = first(na.omit(value))),
      run_id = cumsum(new_seg)
    ) %>%
    filter(!is.na(value)) %>%
    group_by(name, value, run_id, chr) %>%
    summarise(
      i_start = min(i),
      i_end = max(i),
      pos_start = pos[which.min(i)],
      pos_end = pos[which.max(i)],
      .groups = "drop"
    )
}

pred_raw_paper <- vroom(
  "preds-eap100_ep005_penalty15.txt",
  comment = "#",
  delim = "\t"
)
pred_raw_supp <- vroom(
  "preds-eap50_ep005_penalty15.txt",
  comment = "#",
  delim = "\t"
)

pred_merged_paper <- process_pred(pred_raw_paper, pos)
pred_merged_supp <- process_pred(pred_raw_supp, pos)

get_missing <- function(raw) {
  raw$i <- seq_len(nrow(raw))
  raw$pos <- pos$pos
  raw$chr <- 3
  long <- raw %>% pivot_longer(cols = -c("i", "pos", "chr"))
  (long %>%
    group_by(name) %>%
    summarize(x = sum(is.nan(value)) / n()) %>%
    filter(x > 0.1))$name
}

missing_edge <- union(get_missing(pred_raw_paper), get_missing(pred_raw_supp))


palreds <- brewer.pal(n = 9, name = "Reds")[3:9]
palblues <- brewer.pal(n = 9, name = "Blues")[3:9]

dist_breaks <- c(0, 0.25, 0.3, 0.4, 0.5, 0.75, 1)


annotate_merged <- function(pred_merged) {
  merge(merge(merge(pred_merged, distances), taxon), ttaxon) %>%
    filter(
      taxon %in%
        c("Carnivora", "Chiroptera", "Primates", "Artiodactyla", "Rodentia")
    ) %>%
    filter(!name %in% missing_edge) %>%
    mutate(ie = as.numeric(sub("I", "", name)))
}

df_paper <- annotate_merged(pred_merged_paper)
df_supp <- annotate_merged(pred_merged_supp)


ggplot() +

  geom_linerange(
    data = df_paper %>%
      distinct(staxon, ie, taxon) %>%
      mutate(xmin = 0, xmax = 200e6),
    aes(
      xmin = xmin,
      xmax = xmax,
      y = reorder(interaction(staxon, ie), ie)
    ),
    color = "darkgray",
    alpha = 0.4,
    linewidth = 0.05
  ) +

  geom_segment(
    data = df_paper %>% filter(value > 0),
    aes(
      x = pos_start,
      xend = pos_end,
      y = reorder(interaction(staxon, ie), ie),
      yend = reorder(interaction(staxon, ie), ie),
      color = cut(distance, dist_breaks)
    ),
    linewidth = 0.8,
    position = position_nudge(y = 0.2)
  ) +
  scale_color_manual(
    values = palreds,
    name = "Distance\n(paper, eap100)",
    drop = FALSE
  ) +

  new_scale_colour() +

  geom_segment(
    data = df_supp %>% filter(value > 0),
    aes(
      x = pos_start,
      xend = pos_end,
      y = reorder(interaction(staxon, ie), ie),
      yend = reorder(interaction(staxon, ie), ie),
      color = cut(distance, dist_breaks)
    ),
    linewidth = 0.8,
    position = position_nudge(y = -0.2)
  ) +
  scale_color_manual(
    values = palblues,
    name = "Distance\n(supp, eap50)",
    drop = FALSE
  ) +

  facet_wrap(~taxon, scales = "free_y", ncol = 1, space = "free_y") +

  scale_x_continuous(labels = format_mb) +

  theme_minimal_grid() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length.y = unit(0, "pt"),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.y = element_text(size = 6, margin = margin(r = -20))
  ) +
  guides(color = guide_legend(nrow = 1)) +
  labs(
    y = "Clade",
    x = "Coordinate",
    title = "Human chr3 — paper (red, upper) vs supp (blue, lower)"
  )
