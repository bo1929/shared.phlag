require(ggplot2)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)
require(vroom)
require(tidyquant)

taxon <- vroom("taxon_map_order.txt", comment = "#", delim = "\t", col_names = c("name", "taxon"))
ttaxon <- vroom("taxon_map.txt", comment = "#", delim = "\t", col_names = c("name", "staxon"))

pos <- vroom("pos", col_names = c("pos"), delim = "\t")
qqs <- vroom("bs-qqs.txt")

format_mb <- function(l) {
  # turn in to character string in scientific notation
  # l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  # l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  # l <- gsub("e", "%*%10^", l)
  # return this as an expression
  l <- paste0(round(l / 1e6), " Mb")
  return(l)
}

qqs<-qqs %>%
  pivot_longer(cols = -c("clade", "y"), names_to = "i") %>%
  mutate(i=as.numeric(sub("g_", "", i)))
pos$i <- 1:nrow(pos)
qqs<-merge(qqs, pos)
merge(qqs, ttaxon, by.x="clade", by.y="name") %>% 
  # mutate(pos=max(pos)-pos) %>%
  mutate(region=if_else(pos<55e6 & pos>35e6, "40Mb-55Mb", "rest")) %>%
  # mutate(region=if_else(pos<40e6 & pos>30e6, "30Mb-40Mb", region)) %>%
  ggplot() +
  aes(x=value, color=region) +
  facet_grid(rows = vars(interaction(clade, staxon)), cols=vars(y)) +
  stat_ecdf() + 
  scale_color_brewer(palette = "Set1") +
  theme_bw()

qqs$y <- factor(
  qqs$y,
  levels = c(1, 2, 3), 
  labels = c("Main", "Alternative 1", "Alternative 2")
)

merge(qqs, ttaxon, by.x="clade", by.y="name") %>% 
  filter(clade == "I160") %>%
  # mutate(pos=cut(pos, 400)) %>%
  # group_by(pos, clade, y) %>%
  # filter(y==1) %>%
  arrange(i) %>%
  # summarise(value=mean(value)) %>%
  # mutate(region=if_else(pos<55e6 & pos>35e6, "40Mb-55Mb", "rest")) %>%
  # mutate(region=if_else(pos<40e6 & pos>30e6, "30Mb-40Mb", region)) %>%
  ggplot() +
  aes(x=pos, y=value, color=as.factor(y)) +
  facet_wrap(~interaction(staxon, clade), ncol = 1) +
  geom_ma(n = 100,  na.rm=F, size=1, linetype = "solid") +
  # stat_summary_bin(bins = 500, size=0.05) + 
  scale_color_brewer(palette = "Set1", direction = -1) +
  scale_x_continuous( labels = format_mb) +
  theme(axis.text.x = element_text(angle=90)) +
  theme_bw() +
  xlim(c(25e6, 70e6)) +
  theme_cowplot() +
  labs(x="Coordinate", y="QQS", color="Topology")
ggsave("~/Downloads/ex.pdf", width = 12, height = 3)


predsu <- merge(merge(merge(pred_merged, distances), taxon), ttaxon)  %>% filter(name %in% unique(lqqs$clade)) %>% filter(!name %in% c("I219", "I225") )
names(predsu)[1] <- "clade"
merge(qqs%>% filter(!clade %in% c("I219", "I225") ), ttaxon, by.x="clade", by.y="name") %>% 
  arrange(i) %>%
  ggplot() +
  aes(x=pos, y=value, color=as.factor(y)) +
  geom_ma(n = 100,  na.rm=F, size=0.25, linetype = "solid") +
  # stat_summary_bin(bins = 500, size=0.05) + 
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous( labels = format_mb) +
  theme(axis.text.x = element_text(angle=90)) +
  # geom_tile(data=predsu, aes(x=pos, y=value), alpha=0.5, inherit.aes = F) +
  geom_rect(data=predsu %>% filter(value==1), aes(ymin=0, ymax=1, xmin=pos_start, xmax=pos_end, fill=cut(distance, c(0, 0.1, 0.2, 0.3, 0.4, 0.5) )), alpha =0.5, inherit.aes = F) +
  facet_wrap(~interaction(staxon, clade), ncol = 1, ) +
  theme_bw() +
  scale_fill_manual(values = c("#666666", "#bbbbbb")) +
  theme(legend.position = "bottom") +
  labs(x="Coordinate", y="QQS", color="Topology", fill="Distance")
ggsave("qqs-target_branches-highlighted.pdf", width = 10, height = 10)

merge(qqs, ttaxon, by.x="clade", by.y="name") %>% 
  mutate(region=if_else(pos<55e6 & pos>35e6, "35Mb-55Mb", "rest")) %>%
  ggplot() +
  aes(x=interaction(staxon, clade), y=value, color=region) +
  facet_wrap(~y) +
  geom_boxplot(outliers = FALSE) + 
  scale_color_brewer(palette = "Set1") +
  theme_bw()

merge(qqs, ttaxon, by.x="clade", by.y="name") %>% 
  mutate(region=if_else(pos<55e6 & pos>35e6, "35Mb-55Mb", "rest")) %>%
  # filter(region!="rest")%>%
  ggplot() +
  aes(x=i, y=value) +
  facet_grid(rows = vars(clade), cols=vars(y)) +
  geom_ma(na.rm=FALSE, ma_fun = SMA, n = 30, color = "black", size=0.15, linetype = "solid") +
  scale_color_brewer(palette = "Set1") +
  theme_bw()

pred <- vroom("preds-eap100_ep005_penalty15.txt", comment = "#", delim = "\t")
distances <- vroom("../mammals/distances_chr3.txt", col_names = c("name", "distance"))
pred$i <- 1:nrow(pred)
pred$pos <- pos$pos
pred$chr <- 3
pred <- pred %>% pivot_longer(cols = -c("i", "pos", "chr"))

lqqs <- merge(qqs, ttaxon, by.x="clade", by.y="name")
predsu <- pred %>% filter(name %in% unique(lqqs$clade))
predsu <- merge(predsu, distances)

ghdisc <- c("I227", "I228", "I229", "I230", "I231", "I226", "I225", "I217")

pred <- pred %>%
  # group_by(name, chr) %>%
  group_by(name) %>%
  mutate(p=mean(value, na.rm = TRUE), value=value) %>% ungroup() %>%
  mutate(
    value = ifelse(
      p > 0.5, 1-value, value
    )
  )

wa <- (
  pred %>% group_by(name) %>%
    summarise(avg_p=sum(value, na.rm=TRUE)/n()) %>%
    filter(avg_p > 0.10)
)$name

missing_edge <- (pred %>% group_by(name) %>%
                   summarize(x=sum(is.nan(value)/n())) %>%
                   filter(x>0.1))$name

wa<- (pred %>% group_by(name) %>%
  summarize(pval=mean(value, na.rm=TRUE)) %>% filter(pval > 0.05))$name

pred_merged <- pred %>% 
  arrange(name, i) %>%
  group_by(name) %>%
  mutate(
    value_filled = zoo::na.locf(value, na.rm = FALSE),
    new_seg = value_filled != lag(value_filled, default = first(value_filled)),
    run_id = cumsum(new_seg)
  ) %>%
  filter(!is.na(value_filled)) %>%
  group_by(name, value = value_filled, run_id, chr) %>%
  summarise(
    i_start   = min(i),
    i_end     = max(i),
    pos_start = pos[which.min(i)],
    pos_end   = pos[which.max(i)],
    .groups = "drop"
  )
anomalous_edge <- unique(
  (pred_merged %>% filter((i_end - i_start) > 3000 & value == 1))$name
  )

axis_df <- pred %>%
  arrange(i) %>%
  distinct(i, pos)
break_idx <- scales::pretty_breaks(6)(range(axis_df$i))
break_pos <- axis_df$pos[match(break_idx, axis_df$i)]
break_pos[1] <- 0
break_pos[length(break_pos)] <- max(axis_df$pos)

palreds <- brewer.pal(n = 9, name = "Reds")[3:9]
merge(merge(merge(pred_merged, distances), taxon), ttaxon) %>%
  # filter(!name %in% anomalous_edge) %>% 
  # filter(!name %in% anomalous_edge) %>% 
  filter(value>0) %>% 
  # filter(distance>0.25) %>% 
  # filter(name %in% c("I222", "I219", "I217", "I223")) %>% 
  filter(taxon %in% c("Carnivora", "Chiroptera", "Primates", "Artiodactyla", "Rodentia")) %>%
  # filter(name %in% wa) %>% 
  filter(!name %in% missing_edge) %>% 
  # filter(name %in% ghdisc) %>%
  # filter(!name %in% missing_edge) %>%
  # filter(distance > 0.2 & pval > 0.1) %>%
  # filter(distance > 0.2) %>%
  mutate(ie=as.numeric(sub("I", "", name))
         ) %>%
  ggplot() +
  aes(x = pos_start,
      xend = pos_end,
      y = reorder(interaction(staxon,ie), ie),
      yend = reorder(interaction(staxon,ie), ie),
      # y = nn,
      # yend = nn,
      color = cut(distance, c(0, 0.25, 0.3, 0.4, 0.5, 0.75, 1))) +
  facet_wrap(~taxon, scales = "free_y", ncol = 1, space = "free_y") +
  geom_linerange(aes(
    xmax=0,xmin=200e6, y=reorder(interaction(staxon,ie), ie)), color="darkgray",  alpha=0.5, linewidth = 0.05, position = position_identity())+
  geom_segment(linewidth = 1.5) +
  #scale_color_manual(
  #  values = c("#fb6a4a", "#ca0020"),
  #) + 
  #scale_color_viridis_d( direction = -1,option = "G",begin  = 0.25) +
  # scale_color_brewer(palette = "Reds") +
  scale_color_manual(values =palreds ) +
  # facet_wrap(~chr, scales = "free_x") +
  scale_x_continuous(
    # breaks = break_idx[1:length(break_idx)-1],
    # sec.axis = sec_axis(
    #   ~ .,
    #  name = "Coordinate (chr1)",
    #  breaks = break_idx[1:length(break_idx)-1],
    labels = format_mb
    #  )
  ) +
  theme_minimal_grid() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.text.y = element_text(size = 5),
    axis.ticks.length.y.left = unit(0, "pt"),
    axis.ticks.length.y.right = unit(0, "pt"),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.y = element_text(size=6, margin = margin(r=-20) )
  ) +
  scale_alpha_continuous(guide = "none") +
  guides(color = guide_legend(nrow = 1)) +
  labs(y="Clade", color="Distance", x="Coordinate", title="Human chr3")
# lbl <- vroom("tmp/mapping_stiller_fig2a.tsv", col_names = c("name", "stiller_label"))
ggsave("mammals-all-large.pdf", width = 11, height = 13)

ggsave("mammals-dist_gt025-large.pdf", width = 12, height = 8)
ggsave("mammals-dist_gt025-compact.pdf", width = 9, height = 6)

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

merge(merge(pred, distances), taxon) %>%
  filter(!name %in% missing_edge) %>%
  filter(!name %in% anomalous_edge) %>% 
  # filter(value>0) %>%
  # filter(distance > 0.3) %>%
  # filter(distance > 0.1 & pval > 0.1) %>%
  # filter(distance > 0.2), p) %>%
  # filter(value > 0 ) %>% 
  # filter(p > 0.1) %>%
  # filter(taxon %in% c("Carnivora", "Chiroptera", "Primates", "Artiodactyla", "Rodentia")) %>%
  mutate(ibin=cut(i, 100), ie=as.numeric(sub("I", "", name))) %>%
  # mutate(pos_bin = cut(pos, 100)) %>%
  group_by(ibin, taxon) %>%
  summarize(pos=mean(pos), p=(mean(value & distance > 0.2, na.rm=TRUE))) %>%
  # filter(p > 0.5) %>%
  # group_by(pos_bin, taxon) %>%
  # summarise(p=sum(n)/n()) %>%
  ggplot() +
  aes(x=pos, color=taxon, y=p) +
  # facet_wrap(~taxon, scale="free_y") +
  geom_point() +
  # scale_x_binned(n.breaks = 1000) +
  theme_minimal_vgrid() +
  theme(axis.text.y = element_text(size = 5)) +
  # scale_fill_brewer(palette = "Reds") +
  labs(fill="Distance")

# phylopic
get_phylopic_image <- function(taxon_name) {
  img <- tryCatch({
    print(taxon_name)
    uuid <- get_uuid(name = taxon_name, n = 1)
    img <- get_phylopic(uuid = uuid)
  }, warning = function(w) {
    img <- NULL
  }, error = function(e) {
    img <- NULL
  }, finally = {
    img <- NULL
  })
  return(img)
  
}
# install.packages("ggimage")
library(rphylopic)
library(ggimage)

taxon_image_mapping <- tibble(
  name = unique(dfplot$TaxonName),
  image = sapply(sub(" \\+.*", "", unique(dfplot$TaxonName)), get_phylopic_image)
)
names(taxon_image_mapping)[1] <- "TaxonName"
dfimg <- dfplot %>%
  distinct(TaxonName) %>%
  left_join(taxon_image_mapping, by = "TaxonName") %>%
  mutate(
    x = min(dfplot$i_start) - 0.05 * max(dfplot$i_end),
    y = TaxonName
  )


ggplot(dfplot,
       aes(x = i_start, xend = i_end, y = name, yend = name, color = factor(value))) +
  # geom_segment(linewidth = 1) +
  scale_color_manual(values = c("0" = "black", "1" = "red"), guide = "none") +
  scale_y_discrete(name = "Taxon / Label", labels = y_labels) +
  scale_x_continuous(
    name = "Index (i)",
    sec.axis = sec_axis(~ ., name = "Genomic position",
                        breaks = break_idx,
                        labels = scales::comma(break_pos))
  ) +
  geom_phylopic(
    img = dft$image, data = dft , aes(x = x, y = y), inherit.aes = FALSE) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )



dft<-subset(dfimg, !mapply(is.null, image)) %>% mutate(TaxonName=sub(" \\+.*", "", TaxonName))






dqqs %>%
  mutate(l=) %>%
  ggplot() +
  aes(x = i, y = I210) +
  facet_wrap(~y, nrow = 3) +
  # geom_hline(yintercept = 0.33) +
  geom_line(size=0.05) +
  # geom_area(aes(x = i, y = l), fill="red", alpha=1) + 
  # geom_line(aes(x = i, y = l), color="red", alpha=1) + 
  geom_ma(aes(x = i, y = l), ma_fun = SMA, n = 5, color = "red", size=0.25, linetype = "solid") +
  geom_ma(ma_fun = SMA, n = 50, color = "blue", size=0.25, linetype = "solid") +
  # geom_ma(aes(x = i, y = l), ma_fun = SMA, n = 50, color = "red", size=0.25, linetype = "solid") +
  theme_cowplot() +
  scale_x_continuous()# + 
coord_cartesian(ylim = c(0, 1), xlim = c(15000, 16000))

y<- vroom("/Users/asapci/Desktop/region-segmentation/CASTER-data/Biological/Mammal241/SlidingWindow/x", col_names = c("b", "chr", "pos", "a",  "top", "v1", "v2"))

df %>%
  pivot_longer(cols = -c("clade", "y")) %>%
  group_by(clade, y) %>%
  summarise(variance = var(value))