# Usage:
#   Rscript plot_clade.R <tree_file> <label1,label2,...> <output.pdf>

library(ape)
library(ggtree)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_clade.R <tree_file> <label1,label2,...> <title> <output.pdf>")
}

tree_path <- args[1]
labels    <- trimws(unlist(strsplit(args[2], ",")))
title    <- args[3]
out_file  <- args[4]

tree <- tryCatch(read.tree(tree_path), error = function(e) read.nexus(tree_path))

n_tips        <- length(tree$tip.label)
node_numbers  <- seq_along(tree$node.label) + n_tips
matched_nodes <- node_numbers[tree$node.label %in% labels]

if (length(matched_nodes) == 0) {
  stop("No internal node labels matched.\n",
       "Provided: ", paste(labels, collapse = ", "), "\n",
       "Node labels (first 20): ", paste(head(tree$node.label, 20), collapse = ", "))
}

all_desc_tips <- unique(unlist(lapply(matched_nodes, function(nd) {
  extract.clade(tree, nd)$tip.label
})))
mrca_node <- if (length(all_desc_tips) >= 2) getMRCA(tree, all_desc_tips) else matched_nodes[1]

parent_of <- function(tree, node) tree$edge[tree$edge[, 2] == node, 1]
parent     <- parent_of(tree, mrca_node)

outgroup_tip <- if (length(parent) > 0) {
  sister_tips <- setdiff(extract.clade(tree, parent)$tip.label, all_desc_tips)
  if (length(sister_tips) > 0) sister_tips[1] else NULL
} else NULL

keep_tips  <- c(all_desc_tips, outgroup_tip)
clade_tree <- keep.tip(tree, keep_tips)

ct_node_numbers <- seq_along(clade_tree$node.label) + length(clade_tree$tip.label)
highlight_nodes <- ct_node_numbers[clade_tree$node.label %in% labels]

n_plot_tips <- length(clade_tree$tip.label)
fig_h_mm    <- max(60, n_plot_tips * 3) * 1.25
fig_w_mm    <- 180

tip_pt   <- max(10, min(16, fig_h_mm / 18)) * 0.75
node_pt  <- max(8,  min(13, fig_h_mm / 20)) * 0.85
title_pt <- max(12, min(17, fig_h_mm / 15)) * 0.85

p_base <- ggtree(clade_tree,
                 layout        = "rectangular",
                 branch.length = "none",
                 size          = 0.5,
                 color         = "grey50")

td <- p_base$data   # columns: node, parent, x, y, label, isTip, ...

edges <- as.data.frame(clade_tree$edge)
colnames(edges) <- c("parent", "child")

xy        <- td[, c("node", "x", "y", "label", "isTip")]
edges     <- merge(edges, xy, by.x = "child",  by.y = "node") 
edges     <- merge(edges, xy[, c("node","x")], by.x = "parent", by.y = "node",
                   suffixes = c("", "_parent"))

edge_labs <- edges[!edges$isTip & !is.na(edges$label) & edges$label != "", ]
edge_labs$x_mid <- (edge_labs$x + edge_labs$x_parent) / 2

x_range <- max(td$x) - min(td$x)
x_left  <- min(td$x) - x_range * 0.12
x_right <- max(td$x) + x_range * 0.55

p <- p_base +
  geom_tree(aes(color = I(ifelse(node %in% highlight_nodes, "#090982", "grey50"))), size = 0.6) +
  geom_label(
    data     = edge_labs,
    aes(x = x_mid, y = y+0.7, label = label),
    size      = node_pt / .pt,
    color     = "black",
    fill      = "white",
    fontface  = "bold",
    alpha=0,
    label.size     = 0,        # no border
    label.padding  = unit(0.1, "lines"),
    label.r        = unit(0.1, "lines")
  ) +
  
  # Tip labels
  geom_tiplab(
    aes(
      label    = gsub("_", " ", label),
      fontface = ifelse(label == outgroup_tip, "italic", "italic"),
      color    = ifelse(label == outgroup_tip, "grey40", "grey15")
    ),
    size        = tip_pt / .pt,
    offset      = 0.1,
    show.legend = FALSE
  ) +
  
  coord_cartesian(xlim = c(x_left, x_right), clip = "off") +
  scale_color_identity() +
  theme_tree() +
  theme(
    plot.title      = element_text(size = title_pt, face = "bold",
                                   color = "grey15", margin = margin(b = 8)),
    plot.margin     = margin(t = 12, r = 4, b = 6, l = 10, unit = "mm"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  ggtitle(paste(title, paste(labels, collapse = ", "), sep=": "))

ggsave(out_file, p,
       width = fig_w_mm, height = fig_h_mm, units = "mm",
       device = cairo_pdf, limitsize = FALSE)

message(sprintf("Saved → %s  [%.0f x %.0f mm, %d tips, outgroup: %s]",
                out_file, fig_w_mm, fig_h_mm, n_plot_tips,
                if (!is.null(outgroup_tip)) outgroup_tip else "none"))