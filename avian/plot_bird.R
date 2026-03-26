require(ggplot2)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)
require(vroom)
require(tidyquant)

pos <- vroom("positions-gene_trees-Stiller2024-concat-sorted.txt", col_names = c("chr", "pos"), delim = "\t")
pred_woZ <- vroom("preds-concat-eap100_ep005_penalty15-chr_1_5.txt", comment = "#", delim = "\t")
distances_woZ <- vroom("distances_concat-chr_1_5.txt", col_names = c("name", "distance"))
pred_wZ <- vroom("preds-concat-eap100_ep005_penalty15-chr_1_5_Z.txt", comment = "#", delim = "\t")
distances_wZ <- vroom("distances_concat-chr_1_5_Z.txt", col_names = c("name", "distance"))

#segment
pred_woZ <- vroom("eap40_ep005_noprior-concat_1_5/preds-concat-eap100_ep005_noprior.txt", comment = "#", delim = "\t")
distances_woZ <- vroom("eap40_ep005_noprior-concat_1_5/distances_concat.txt", col_names = c("name", "distance"))
pred_wZ <- vroom("eap40_ep005_noprior-concat_1_5_Z/preds-1_5_Z-eap100_ep005_noprior.txt", comment = "#", delim = "\t")
distances_wZ <- vroom("eap40_ep005_noprior-concat_1_5_Z/distances_chr1_5_Z.txt", col_names = c("name", "distance"))

pred_woZ$i <- 1:nrow(pred_woZ)
pred_woZ$pos <- (pos %>% filter(chr != "chrZ"))$pos
pred_woZ$chr <- (pos %>% filter(chr != "chrZ"))$chr
pred_woZ <- pred_woZ %>% pivot_longer(cols = -c("i", "pos", "chr"))

pred_wZ$i <- 1:nrow(pred_wZ)
pred_wZ$pos <- pos$pos
pred_wZ$chr <- pos$chr
pred_wZ <- pred_wZ %>% pivot_longer(cols = -c("i", "pos", "chr"))

pred_wZ$with_chrZ <- T
pred_woZ$with_chrZ <- F
pred <- rbind(pred_wZ, pred_woZ)

distances_wZ$with_chrZ <- T
distances_woZ$with_chrZ <- F
distances <- rbind(distances_wZ, distances_woZ)

ghdisc <- c(9, 11, 13, 21, 22, 23, 26, 27, 45, 44, 39, 43, 42, 37, 33, 5, 6, 62, 61, 60, 53, 49)

pred <- pred %>%
  group_by(name, with_chrZ) %>%
  mutate(p=mean(value, na.rm = TRUE), value=value) %>% ungroup() %>%
  mutate(
    value = ifelse(
      p > 0.5, 1-value, value
    )
  )

ann <- vroom("mapping_stiller_fig2a.tsv", col_names = c("name", "label_stiller"))

wa <- (
  pred %>% group_by(name) %>%
    summarise(avg_p=sum(value, na.rm=TRUE)/n()) %>%
    filter(avg_p > 0.05)
)$name

missing_edge <- (pred %>% group_by(name) %>%
                   summarize(x=sum(is.nan(value)/n())) %>%
                   filter(x>0.1))$name

pred_merged <- pred %>%
  arrange(name, i) %>%
  group_by(name) %>%
  mutate(
    value = (function(v) {
      rle_v   <- rle(is.na(v))
      fill_runs <- rle_v$values & rle_v$lengths < 10
      fill_idx  <- rep(fill_runs, rle_v$lengths)
      v[fill_idx] <- zoo::na.locf(v, na.rm = FALSE)[fill_idx]
      v
    })(value),
    new_seg = is.na(value) |
      is.na(lag(value)) |
      value != lag(value, default = first(na.omit(value))),
    run_id  = cumsum(new_seg)
  ) %>%
  filter(!is.na(value)) %>%
  group_by(name, value, run_id, chr) %>%
  summarise(
    i_start   = min(i),
    i_end     = max(i),
    pos_start = pos[which.min(i)],
    pos_end   = pos[which.max(i)],
    .groups   = "drop"
  )
anomalous_edge <- unique((pred_merged %>%
                            filter((i_end - i_start) > 2000 & value == 1))$name)

axis_df <- pred %>%
  arrange(i) %>%
  distinct(i, pos)
break_idx <- scales::pretty_breaks(6)(range(axis_df$i))
break_pos <- axis_df$pos[match(break_idx, axis_df$i)]
break_pos[1] <- 0
break_pos[length(break_pos)] <- max(axis_df$pos)

format_mb <- function(l) {
  # turn in to character string in scientific notation
  # l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  # l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  # l <- gsub("e", "%*%10^", l)
  # return this as an expression
  l <- paste0(round(l / 1e6), "")
  return(l)
}

library(RColorBrewer)
palreds <- brewer.pal(n = 6, name = "Reds")[1:6]
palblues <- brewer.pal(n = 6, name = "Blues")[1:6]

library(ggnewscale)
dfplot <- merge(pred_merged, ann, by.x = "name", by.y="name")
dfplot2<-merge(dfplot, distances)
dfplot2 <- dfplot2 %>% group_by(chr) %>% mutate(lim=max(pos_end)) %>% ungroup()
dfplot2<-dfplot2 %>% filter(!label_stiller %in% c("49", "53"))
dfplot2 %>% # filter(!name %in% anomalous_edge) %>%
  # filter(with_chrZ) %>%
  filter(!grepl("P", label_stiller)) %>%
  arrange(as.numeric(label_stiller)) %>%
  filter(label_stiller %in% ghdisc) %>% 
  arrange(as.numeric(label_stiller)) %>%
  filter(chr != "chrZ") %>%
  # filter(grepl("Homo", TaxonName) | grepl("Homi", TaxonName)) %>%
  ggplot() +
  # facet_wrap(~TaxonName, scales = "free_y") +
  geom_rect(aes(
    xmax=lim,xmin=min(pos_start),
    y=reorder(label_stiller, as.numeric(label_stiller))),fill=NA,height=1,color="darkgray",size=0.2)+
  # geom_segment(data= dfplot2%>%
  #                filter(!name %in% c("53", "49")) %>%
  #                filter(value>0) %>% filter(!grepl("P", label_stiller)) %>%
  #                arrange(as.numeric(label_stiller)) %>%
  #                filter(label_stiller %in% ghdisc) %>% filter(with_chrZ),
  #   aes(x = pos_start,
  #       xend = pos_end,
  #       y = reorder(label_stiller, as.numeric(label_stiller)),
  #       yend = reorder(label_stiller, as.numeric(label_stiller)),
  #       # y = name,
  #       # yend = name,
  #       color = cut(value*distance, c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1))),
  #   linewidth = 2, position = position_nudge(y = -0.2)) +
  # # scale_color_brewer(palette = "Reds") +
  # scale_color_manual(values = palblues) +
  # labs(y="Clade", color="with chrZ", x="Coordinate", title="Concatenation") +
  # new_scale_colour() +
  geom_segment(data= dfplot2%>%
                 filter(!name %in% c("53", "49")) %>%
                 filter(value>0) %>% filter(!grepl("P", label_stiller)) %>%
                 filter(label_stiller %in% ghdisc) %>% filter(!with_chrZ) %>% arrange(as.numeric(label_stiller)),
               aes(x = pos_start,
                   xend = pos_end,
                   y = reorder(label_stiller, as.numeric(label_stiller)),
                   yend = reorder(label_stiller, as.numeric(label_stiller)),
                   # y = name,
                   # yend = name,
                   color = cut(value*distance, c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1))),
               # linewidth = 2, position = position_nudge(y = 0.2)) +
               linewidth = 3) +
               # scale_color_brewer(palette = "Blues", ) +
  scale_color_manual(values = palreds) +
  #scale_color_manual(
  #  values = c("0" = "#d9dff5", "1" = "#ff9999"),
  #  guide = "none"
  # ) + 
  #scale_color_viridis_d( direction = -1,option = "G",begin  = 0.25) +
  facet_wrap(~chr, scales = "free_x", space = "free_x", shrink = TRUE) +
  scale_x_continuous(
    breaks = c(0, 5e7, 1e8, 1.5e8),
    # breaks = break_idx[1:length(break_idx)-1],
   # sec.axis = sec_axis(
   #   ~ .,
   #  name = "Coordinate (chr1)",
   #  breaks = break_idx[1:length(break_idx)-1],
   labels = format_mb
   #  )
  ) +
  # scale_y_discrete(name = "Clade", labels = y_labels) +
  theme_minimal_grid() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # axis.text.y = element_text(size = 10),
    axis.ticks.length.y.left = unit(0, "pt"),
    axis.ticks.length.y.right = unit(0, "pt"),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    panel.spacing.x = unit(0, "mm"),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(size = 10, margin = margin(r =-6)),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom"
  ) +
  # guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_alpha_continuous(guide = "none") +
  labs(y="Clade", color="Distance", x="Coordinate (Mb)", title="Concatenation & updated prior")
  # lbl <- vroom("tmp/mapping_stiller_fig2a.tsv", col_names = c("name", "stiller_label"))
# ggsave2("concat_birds-updated_prior-woZ.pdf", width = 12, height = 4.5)

library(ggnewscale)
dfplot <- merge(pred_merged, ann, by.x = "name", by.y="name")
dfplot2<-merge(dfplot, distances)
dfplot2 %>% # filter(!name %in% anomalous_edge) %>%
  # filter(with_chrZ) %>%
  filter(!grepl("P", label_stiller)) %>%
  arrange(as.numeric(label_stiller)) %>%
  filter(label_stiller %in% ghdisc) %>% 
  # filter(grepl("Homo", TaxonName) | grepl("Homi", TaxonName)) %>%
  ggplot() +
  geom_segment(data= dfplot2%>%filter(value>0) %>% filter(!grepl("P", label_stiller)) %>%
                 arrange(as.numeric(label_stiller)) %>%
                 filter(label_stiller %in% ghdisc) %>% filter(!with_chrZ),
               aes(x = pos_start,
                   xend = pos_end,
                   y = reorder(label_stiller, as.numeric(label_stiller)),
                   yend = reorder(label_stiller, as.numeric(label_stiller)),
                   # y = name,
                   # yend = name,
                   color = cut(value*distance, c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1))),
               linewidth = 2,
               # position = position_nudge(y = 0.2)
               ) +
  scale_color_brewer(palette = "Blues") +
  #scale_color_manual(
  #  values = c("0" = "#d9dff5", "1" = "#ff9999"),
  #  guide = "none"
  # ) + 
  #scale_color_viridis_d( direction = -1,option = "G",begin  = 0.25) +
  facet_wrap(~chr, scales = "free_x",space = "free_x") +
  scale_x_continuous(
    breaks = c(0, 5e7, 1e8, 1.5e8),
    # breaks = break_idx[1:length(break_idx)-1],
    # sec.axis = sec_axis(
    #   ~ .,
    #  name = "Coordinate (chr1)",
    #  breaks = break_idx[1:length(break_idx)-1],
    labels = format_mb
    #  )
  ) +
  # scale_y_discrete(name = "Clade", labels = y_labels) +
  theme_minimal_grid() +
  theme(
    # panel.grid.major.y = element_blank(),
    # panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.ticks.length.y.left = unit(0, "pt"),
    axis.ticks.length.y.right = unit(0, "pt"),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    panel.spacing.y = unit(10, "mm"),
    panel.border = element_blank(),
    axis.line = element_blank()
  ) +
  scale_alpha_continuous(guide = "none") +
  labs(y="Clade", color="without chrZ", x="Coordinate (Mb)", title="Concatenation")
# lbl <- vroom("tmp/mapping_stiller_fig2a.tsv", col_names = c("name", "stiller_label"))


library(ggnewscale)
dfplot <- merge(pred_merged, ann, by.x = "name", by.y="name")
dfplot2<-merge(dfplot, distances)
dfplot2 %>% # filter(!name %in% anomalous_edge) %>%
  # filter(with_chrZ) %>%
  filter(!grepl("P", label_stiller)) %>%
  arrange(as.numeric(label_stiller)) %>%
  filter(label_stiller %in% ghdisc) %>% 
  # filter(grepl("Homo", TaxonName) | grepl("Homi", TaxonName)) %>%
  ggplot() +
  # facet_wrap(~TaxonName, scales = "free_y") +
  geom_segment(data= dfplot2%>%filter(value>0) %>% filter(!grepl("P", label_stiller)) %>%
                 arrange(as.numeric(label_stiller)) %>%
                 filter(label_stiller %in% ghdisc) %>% filter(with_chrZ),
               aes(x = pos_start,
                   xend = pos_end,
                   y = reorder(label_stiller, as.numeric(label_stiller)),
                   yend = reorder(label_stiller, as.numeric(label_stiller)),
                   # y = name,
                   # yend = name,
                   color = cut(value*distance, c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1))),
               linewidth = 2) +
  scale_color_brewer(palette = "Reds") +
  labs(y="Clade", color="with chrZ", x="Coordinate", title="Concatenation") +
  #scale_color_manual(
  #  values = c("0" = "#d9dff5", "1" = "#ff9999"),
  #  guide = "none"
  # ) + 
  #scale_color_viridis_d( direction = -1,option = "G",begin  = 0.25) +
  facet_wrap(~chr, scales = "free_x",space = "free_x") +
  scale_x_continuous(
    breaks = c(0, 5e7, 1e8, 1.5e8),
    # breaks = break_idx[1:length(break_idx)-1],
    # sec.axis = sec_axis(
    #   ~ .,
    #  name = "Coordinate (chr1)",
    #  breaks = break_idx[1:length(break_idx)-1],
    labels = format_mb
    #  )
  ) +
  # scale_y_discrete(name = "Clade", labels = y_labels) +
  theme_minimal_grid() +
  theme(
    # panel.grid.major.y = element_blank(),
    # panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.ticks.length.y.left = unit(0, "pt"),
    axis.ticks.length.y.right = unit(0, "pt"),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    panel.spacing.y = unit(10, "mm"),
    panel.border = element_blank(),
    axis.line = element_blank()
  ) +
  scale_alpha_continuous(guide = "none") +
  labs(y="Clade", color="with chrZ", x="Coordinate (Mb)", title="Concatenation")
# lbl <- vroom("tmp/mapping_stiller_fig2a.tsv", col_names = c("name", "stiller_label"))

dfplot <- merge(pred, ann, by.x = "name", by.y="name")
dfplot <-merge(dfplot, distances)
dfplot <- dfplot %>% filter(label_stiller %in% ghdisc)
dfplot %>% # filter(!name %in% missing_edge & name %in% wa) %>%
  filter(grepl("N", name)) %>%
  mutate(ie=as.numeric(sub("N", "", name))) %>%
  mutate(pos_bin=as.numeric(cut(i, 1000))) %>%
  group_by(pos_bin, chr, label_stiller) %>%
  summarise(portion=mean(value, na.rm = TRUE)) %>%
  ggplot() +
  aes(x=pos_bin, y=label_stiller, fill=portion) +
  geom_tile() +
  facet_wrap(~chr, scale="free_x") +
  scale_x_continuous(n.breaks = 10) +
  theme(axis.text.y = element_text(size = 5, margin = margin(r =-25)))+
  scale_fill_viridis_b() +
  theme_minimal_grid() +
  labs(fill="Distance")

dfplot %>%
  group_by(pos_bin, label_stiller, chr) %>%
  mutate(pos_bin=as.numeric(cut(pos, 1000))) %>%
  summarise(portion=mean(value, na.rm = TRUE), dist=mean(distance*value, na.rm=TRUE)) %>%
  ggplot() +
  facet_wrap(~chr, scale="free_x") +
  aes(x=pos_bin, y=label_stiller, fill=cut(dist, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1))) +
  geom_tile() +
  scale_x_continuous(n.breaks = 10) +
  theme_minimal_grid() +
  theme(axis.text.y = element_text(size = 11)) +
  scale_fill_brewer(palette = "Reds") +
  labs(fill="Distance")
