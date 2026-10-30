#!/usr/bin/env Rscript
# Step 8 - the six TP53 report figures (R edition).
# Backend: system R 4.6 + ggplot2/patchwork; follows the nature-figure contract
# (no "Figure N" prefix inside the artwork, no in-figure footnotes).

ROOT <- normalizePath(Sys.getenv("TP53_ROOT", "."), mustWork = FALSE)
local_lib <- file.path(ROOT, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(grid)
  library(ape); library(ggrepel); library(jsonlite); library(png); library(cowplot)
})
TBL  <- file.path(ROOT, "results", "tables")
OUT  <- file.path(ROOT, "assets", "r")
PDF  <- file.path(OUT, "pdf")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(PDF, recursive = TRUE, showWarnings = FALSE)

SKILL <- path.expand("~/.agents/skills/nature-figure/scripts")
source(file.path(SKILL, "panel_alignment.R"))

INK   <- "#1F2933"; MUTED <- "#7B8794"; GRIDC <- "#E4E7EB"
HOT   <- "#B3261E"; TEAL  <- "#2F7F73"; BLUE  <- "#335C81"
DOM   <- c(TAD = "#E8A87C", DBD = "#3E8E7E", TET = "#4F6D8C", Other = "#D9D9D9")
DOMAINS <- tibble::tribble(
  ~domain, ~start, ~end,
  "TAD", 1, 44, "Other", 45, 101, "DBD", 102, 292,
  "Other", 293, 324, "TET", 325, 356, "Other", 357, 393)
SPECIES <- c(P53_HUMAN = "Human", P53_BOVIN = "Cattle", P53_MOUSE = "Mouse",
             P53_CHICK = "Chicken", P53_XENLA = "Frog", P53_DANRE = "Zebrafish")
ORDER <- unname(SPECIES[c("P53_HUMAN", "P53_BOVIN", "P53_MOUSE",
                          "P53_CHICK", "P53_XENLA", "P53_DANRE")])
pt <- function(x) x / .pt

theme_nf <- function(base = 6.5) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = INK),
      axis.ticks = element_line(linewidth = 0.3, colour = INK),
      axis.text = element_text(size = base - 1.0, colour = INK),
      axis.title = element_text(size = base - 0.4, colour = INK),
      plot.title = element_text(size = base + 1, face = "bold", hjust = 0,
                                margin = margin(b = 1.5)),
      plot.title.position = "plot",
      plot.tag = element_text(size = base + 1, face = "bold", colour = INK),
      plot.tag.position = c(0.01, 0.99),
      legend.title = element_text(size = base - 0.9),
      legend.text = element_text(size = base - 1.3),
      legend.key.size = unit(2.4, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm")
    )
}
theme_set(theme_nf())

save_r <- function(p, name, w_mm, h_mm, multi = TRUE, exemptions = list()) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  if (multi) {
    # patchworkGrob needs a real device with the final physical size; a temporary
    # PDF probe avoids the stray Rplots.pdf that R's default device would create.
    probe <- tempfile(fileext = ".pdf")
    grDevices::pdf(probe, width = w, height = h, useDingbats = FALSE)
    tryCatch(require_patchwork_panel_alignment(
      p,
      manifest_path = file.path(OUT, paste0(name, ".alignment-layout.json")),
      report_path = file.path(OUT, paste0(name, ".alignment.json")),
      overlay_svg = file.path(OUT, paste0(name, ".alignment.svg")),
      width_in = w, height_in = h,
      exemptions = exemptions,
      audit_script = file.path(SKILL, "audit_panel_alignment.py"),
      python = Sys.which("python3"),
      tolerance_pt = 1.5, gutter_tolerance_pt = 1.5, strict = TRUE),
      error = function(e) message("patchwork alignment helper unavailable: ", conditionMessage(e)))
    grDevices::dev.off()
    unlink(probe)
  }
  ragg::agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
                units = "in", res = 600, background = "white")
  print(p); dev.off()
  grDevices::pdf(file.path(PDF, paste0(name, ".pdf")), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE)
  print(p); dev.off()
  message("wrote ", name)
}

parse_clustal <- function(path) {
  lines <- readLines(path, warn = FALSE)
  ids <- character(); seqs <- list()
  for (ln in lines) {
    if (!nzchar(trimws(ln))) next
    if (grepl("^(CLUSTAL|MUSCLE|#)", ln)) next
    if (grepl("^\\s", ln)) next
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(parts) < 2) next
    id <- parts[1]; s <- gsub("[0-9]", "", parts[2])
    if (is.null(seqs[[id]])) { seqs[[id]] <- ""; ids <- c(ids, id) }
    seqs[[id]] <- paste0(seqs[[id]], s)
  }
  seqs
}

load_data <- function() {
  cols <- read_csv(file.path(TBL, "alignment_columns.csv"), show_col_types = FALSE)
  cons <- read_csv(file.path(TBL, "conservation.csv"), show_col_types = FALSE)
  hot  <- read_csv(file.path(TBL, "hotspots.csv"), show_col_types = FALSE)
  cod  <- read_csv(file.path(TBL, "codon_mutation_counts.csv"), show_col_types = FALSE)
  pid  <- read_csv(file.path(TBL, "pairwise_identity.csv"), show_col_types = FALSE)
  list(cols = cols, cons = cons, hot = hot, cod = cod, pid = pid)
}

# ------------------------------------------------------------------ figure 1 --
figure1 <- function() {
  d <- load_data()
  seqs <- parse_clustal(file.path(ROOT, "results", "alignment", "tp53_clustalo.aln"))
  ncol <- nrow(d$cols)
  cols <- d$cols$human_pos
  pos_to_col <- setNames(seq_len(ncol), cols)
  gap <- lapply(names(SPECIES), function(key) {
    id <- grep(key, names(seqs), value = TRUE, fixed = TRUE)[1]
    s <- strsplit(seqs[[id]], "")[[1]]
    data.frame(species = SPECIES[[key]], col = seq_len(ncol),
               value = ifelse(s == "-", NA_real_, d$cols$conservation))
  })
  heat <- bind_rows(gap)
  heat$species <- factor(heat$species, levels = rev(ORDER))

  p1 <- ggplot(heat, aes(col, species, fill = value)) +
    geom_tile(height = 1) +
    scale_fill_gradientn(colors = c("#F5F7FA", "#D6E4EF", "#8FB8D6", "#4A7FA8", "#173F5F"),
                         limits = c(0.45, 1), na.value = "#E7E2DC",
                         name = "Conservation",
                         guide = guide_colourbar(barwidth = unit(2.4, "mm"),
                                                 barheight = unit(34, "mm"),
                                                 title.position = "top",
                                                 title.hjust = 0.5, ticks = FALSE)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, ncol + 0.5)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme(axis.title = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), axis.text.y = element_text(size = 6),
          axis.text.x = element_blank(),
          legend.position = "right", legend.justification = "center",
          legend.box.margin = margin(0, 0, 0, 2), plot.margin = margin(2, 1, 1, 1, "mm"))

  dom <- DOMAINS %>% mutate(xmin = pos_to_col[as.character(start)],
                            xmax = pos_to_col[as.character(end)] + 1)
  p2 <- ggplot(dom) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.05, ymax = 0.95, fill = domain),
              colour = "white", linewidth = 0.2) +
    geom_text(aes(x = (xmin + xmax) / 2, y = 0.5, label = ifelse(as.character(domain) == "Other", "", as.character(domain))),
              colour = "white", fontface = "bold", size = pt(4.6)) +
    scale_fill_manual(values = DOM, guide = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, ncol + 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void()

  prof <- data.frame(col = seq_len(ncol), cons = d$cols$conservation)
  ticks <- pos_to_col[as.character(seq(50, 350, 50))]
  p3 <- ggplot(prof, aes(col, cons)) +
    geom_area(fill = "#4A7FA8", alpha = 0.9) +
    geom_line(colour = "#173F5F", linewidth = 0.2) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, ncol + 0.5),
                       breaks = ticks, labels = names(ticks)) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(x = "Human p53 residue", y = "Conservation") +
    theme(axis.text = element_text(size = 5.6))

  fig <- (p1 / p2 / p3) + plot_layout(heights = c(4.1, 0.34, 1.1), guides = "collect") &
    theme(legend.position = "right")
  fig <- fig +
    plot_annotation(title = "Conservation across six p53 orthologues",
                    theme = theme(plot.title = element_text(size = 7.5, face = "bold")))
  save_r(fig, "fig1_alignment", 180, 108, multi = FALSE)
}

# ------------------------------------------------------------------ figure 2 --
figure2 <- function() {
  d <- load_data()
  tr <- read.tree(file.path(ROOT, "results", "alignment", "tp53_clustalo.nwk"))
  tr$tip.label <- sub("^.*\\|(P53_[A-Z]+)$", "\\1", tr$tip.label)
  n <- length(tr$tip.label); nt <- tr$Nnode
  dep <- node.depth.edgelength(tr)
  y <- numeric(n + nt)
  y[seq_len(n)] <- n:1
  # ape numbers internal nodes in pre-order (root first); children must be
  # resolved before parents, so walk the internal nodes in reverse order.
  for (node in rev((n + 1):(n + nt))) {
    ch <- tr$edge[tr$edge[, 1] == node, 2]
    y[node] <- mean(y[ch])
  }
  seg <- do.call(rbind, lapply(seq_len(nrow(tr$edge)), function(i) {
    p <- tr$edge[i, 1]; ch <- tr$edge[i, 2]
    data.frame(x = dep[p], xend = dep[ch], y = y[ch], yend = y[ch])
  }))
  vert <- do.call(rbind, lapply((n + 1):(n + nt), function(node) {
    ch <- tr$edge[tr$edge[, 1] == node, 2]
    data.frame(x = dep[node], xend = dep[node], y = min(y[ch]), yend = max(y[ch]))
  }))
  tips <- data.frame(x = max(dep) + max(dep) * 0.03, y = y[seq_len(n)],
                     label = unname(SPECIES[tr$tip.label]))
  p1 <- ggplot() +
    geom_segment(data = seg, aes(x, y, xend = xend, yend = yend),
                 colour = INK, linewidth = 0.28) +
    geom_segment(data = vert, aes(x, y, xend = xend, yend = yend),
                 colour = INK, linewidth = 0.28) +
    geom_text(data = tips, aes(x, y, label = label), hjust = 0, size = pt(5.8),
              colour = INK) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(-0.02 * max(dep), max(dep) * 1.4)) +
    scale_y_continuous(expand = expansion(add = 0.6)) +
    coord_cartesian(clip = "off") +
    labs(x = "Branch length", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

  pid <- d$pid %>% mutate(name = unname(SPECIES[species]),
                          name = factor(name, levels = ORDER[-1]))
  long <- pid %>%
    select(name, identity_full, identity_dbd) %>%
    pivot_longer(c(identity_full, identity_dbd), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(ifelse(metric == "identity_full", "Full length", "DBD (102-292)"),
                           levels = c("Full length", "DBD (102-292)")))
  p2 <- ggplot(long, aes(value, name, fill = metric)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.66) +
    geom_text(aes(label = sprintf("%.1f", value)), position = position_dodge(width = 0.72),
              hjust = -0.15, size = pt(5.0), colour = INK) +
    scale_fill_manual(values = c("Full length" = "#A9BFD2", "DBD (102-292)" = TEAL),
                      name = NULL) +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = c(0, 30, 60, 90)) +
    labs(x = "Identity to human (%)", y = NULL) +
    theme(legend.position = "right", legend.direction = "vertical",
          legend.key.height = unit(2.2, "mm"))

  cons <- d$cons %>% mutate(domain = factor(domain, levels = names(DOM)))
  p3 <- ggplot(cons, aes(domain, conservation, fill = domain)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.25,
                 colour = MUTED, alpha = 0.95) +
    geom_jitter(width = 0.14, height = 0, size = 0.25, alpha = 0.25,
                colour = "#20364B") +
    scale_fill_manual(values = DOM, guide = "none") +
    scale_y_continuous(limits = c(0, 1.03), expand = c(0, 0)) +
    labs(x = NULL, y = "Conservation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.2))

  fig <- (p1 | p2 | p3) +
    plot_layout(widths = c(0.95, 1.15, 1.0)) &
    theme(plot.margin = margin(7, 1.5, 1.5, 1.5, "mm"))
  fig <- fig +
    plot_annotation(title = "p53 conservation peaks in the DNA-binding domain",
                    tag_levels = "a")
  save_r(fig, "fig2_cross_species", 180, 64, multi = TRUE)
}

# ------------------------------------------------------------------ figure 3 --
figure3 <- function() {
  d <- load_data()
  subs <- read_csv(file.path(TBL, "hotspot_substitutions.csv"), show_col_types = FALSE)
  summ <- jsonlite::fromJSON(file.path(ROOT, "results", "hotspot_summary.json"))
  hot <- d$hot %>% mutate(codon = paste0(wt, pos))
  codon_levels <- hot$codon[order(-hot$n_missense)]
  aa_class <- c(A="Hydrophobic", V="Hydrophobic", L="Hydrophobic", I="Hydrophobic",
                M="Hydrophobic", F="Hydrophobic", W="Hydrophobic", Y="Hydrophobic",
                S="Polar", T="Polar", N="Polar", Q="Polar",
                K="Basic", R="Basic", H="Basic",
                D="Acidic", E="Acidic", C="Cysteine", G="Gly/Pro", P="Gly/Pro")
  class_col <- c(Hydrophobic="#C2703D", Polar="#2F7F73", Basic="#335C81",
                 Acidic="#C9A227", Cysteine="#7B6D8D", "Gly/Pro"="#7E9E72")

  sub <- subs %>%
    mutate(codon = paste0(wt, pos),
           codon = factor(codon, levels = codon_levels),
           cls = factor(aa_class[as.character(mut)], levels = names(class_col)))
  tot <- sub %>% group_by(mut) %>% summarise(total = sum(n), .groups = "drop")
  mut_levels <- sub %>% distinct(mut, cls) %>% left_join(tot, by = "mut") %>%
    arrange(cls, desc(total)) %>% pull(mut)

  grid <- expand.grid(codon = codon_levels, mut = mut_levels,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
    left_join(sub %>% select(codon, mut, n), by = c("codon", "mut")) %>%
    mutate(n = replace_na(n, 0)) %>%
    left_join(sub %>% distinct(mut, cls), by = "mut") %>%
    mutate(codon = factor(codon, levels = codon_levels),
           mut = factor(mut, levels = mut_levels),
           cls = factor(cls, levels = names(class_col)))

  xlab <- setNames(sprintf("%s\n%d (%.1f%%)", hot$codon, hot$n_missense,
                           100 * hot$n_missense / summ$n_mutations), hot$codon)

  # matrix: every cell is exactly one codon-to-residue substitution; colour
  # encodes only the count, so residue identity and class are not double-coded
  p1 <- ggplot(grid, aes(codon, mut, fill = n)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(data = subset(grid, n >= 25), aes(label = n),
              colour = "white", fontface = "bold", size = pt(4.2)) +
    geom_text(data = subset(grid, n > 0 & n < 25), aes(label = n),
              colour = INK, size = pt(3.9)) +
    facet_grid(cls ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient(low = "#F4F7FA", high = "#173F5F",
                        limits = c(0, max(grid$n)), guide = "none") +
    scale_x_discrete(position = "top", labels = xlab, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x.top = element_text(size = 5.1, face = "bold", lineheight = 0.9),
          axis.text.y = element_text(size = 5.4, face = "bold"),
          strip.background = element_blank(),
          strip.placement = "left",
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 5.0,
                                           face = "bold", colour = INK),
          panel.spacing.y = unit(1.0, "mm"))

  mut <- read_tsv(file.path(ROOT, "data", "mutations", "tp53_mutations_pancan.tsv"),
                  show_col_types = FALSE)
  cls <- c("Missense", "Nonsense", "Frameshift", "Splice site", "In-frame indel")
  mt <- mut$mutationType
  counts <- c(sum(grepl("Missense", mt)), sum(grepl("Nonsense", mt)),
              sum(grepl("Frame_Shift", mt)), sum(grepl("Splice", mt)),
              sum(grepl("In_Frame", mt)))
  other <- nrow(mut) - sum(counts)
  cdf <- data.frame(label = factor(c(cls, "Other"), levels = rev(c(cls, "Other"))),
                    n = c(counts, other))
  p2 <- ggplot(cdf, aes(n, label, fill = label)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%d (%.0f%%)", n, 100 * n / nrow(mut))),
              hjust = -0.08, size = pt(5.0), colour = INK) +
    scale_fill_manual(values = c("Missense" = TEAL, "Nonsense" = HOT,
                                 "Frameshift" = "#D08C60", "Splice site" = "#C9A227",
                                 "In-frame indel" = BLUE, "Other" = "#ADB5BD"),
                      guide = "none") +
    scale_x_continuous(limits = c(0, max(counts) * 1.3), expand = c(0, 0)) +
    labs(x = "TP53 mutations", y = NULL)

  fig <- (p1 | p2) + plot_layout(widths = c(1.18, 1))
  fig <- fig +
    plot_annotation(title = "Five missense hotspots and their substitution spectra",
                    tag_levels = "a")
  save_r(fig, "fig3_hotspots", 180, 84, multi = TRUE)
}

# ------------------------------------------------------------------ figure 4 --
figure4 <- function() {
  d <- load_data()
  cons <- d$cons; cod <- d$cod; hot <- d$hot

  p1 <- ggplot(cons, aes(human_pos, conservation)) +
    annotate("rect", xmin = 102, xmax = 292, ymin = -Inf, ymax = Inf,
             fill = TEAL, alpha = 0.07) +
    geom_area(fill = "#4A7FA8", alpha = 0.9) +
    geom_line(colour = "#173F5F", linewidth = 0.18) +
    geom_segment(data = hot, aes(x = pos, xend = pos, y = 0, yend = 0.72),
                 inherit.aes = FALSE, colour = HOT, linewidth = 0.22, linetype = "22") +
    geom_text(data = mutate(hot, vj = ifelse(pos %in% c(248, 282), 1.28, 1.10)),
              aes(x = pos, y = vj, label = paste0(wt, pos)),
              inherit.aes = FALSE, colour = HOT, fontface = "bold", size = pt(5.0),
              vjust = 0) +
    annotate("text", x = 191, y = 0.10, label = "DBD", colour = TEAL,
             fontface = "italic", size = pt(5.0)) +
    scale_x_continuous(limits = c(0, 393), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.40), expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Conservation") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(3, 1, 0, 1, "mm"))

  dens <- cod %>% filter(n_missense > 0)
  p2 <- ggplot(dens, aes(x = pos)) +
    geom_segment(aes(xend = pos, y = 1, yend = n_missense + 1),
                 colour = "#8FA3B5", linewidth = 0.18) +
    geom_segment(data = hot, aes(x = pos, xend = pos, y = 1, yend = n_missense + 1),
                 inherit.aes = FALSE, colour = HOT, linewidth = 0.55) +
    scale_y_log10(breaks = c(1, 10, 100), labels = c("1", "10", "100"),
                  limits = c(0.8, max(cod$n_missense) * 2.2), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 393), expand = c(0, 0)) +
    labs(x = NULL, y = "Missense per codon") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(colour = GRIDC, linewidth = 0.18))

  dom <- DOMAINS
  p3 <- ggplot(dom) +
    geom_rect(aes(xmin = start, xmax = end + 1, ymin = 0.08, ymax = 0.92, fill = domain),
              colour = "white", linewidth = 0.2) +
    geom_text(aes(x = (start + end + 1) / 2, y = 0.5, label = ifelse(as.character(domain) == "Other", "", as.character(domain))),
              colour = "white", fontface = "bold", size = pt(4.6)) +
    scale_fill_manual(values = DOM, guide = "none") +
    scale_x_continuous(limits = c(0, 394), expand = c(0, 0), breaks = seq(0, 350, 50)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Human p53 residue", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.title.y = element_blank())

  # one stacked figure with a shared x axis (no a/b/c tags, as in Figure 1)
  fig <- (p1 / p2 / p3) + plot_layout(heights = c(1.15, 1.45, 0.26)) +
    plot_annotation(title = "Conservation and somatic mutation density along p53")
  save_r(fig, "fig4_landscape", 180, 88, multi = TRUE)
}

# ------------------------------------------------------------------ figure 5 --
figure5 <- function() {
  d <- load_data()
  summ <- jsonlite::fromJSON(file.path(ROOT, "results", "hotspot_summary.json"))
  bg <- d$cons %>%
    left_join(d$cod %>% select(pos, n_missense), by = c("human_pos" = "pos")) %>%
    mutate(n_missense = coalesce(n_missense, 0L),
           dbd = domain == "DBD",
           grp = factor(ifelse(dbd, "DBD", "Outside DBD"),
                        levels = c("Outside DBD", "DBD", "Hotspot codon"))) %>%
    filter(!is.na(conservation))

  p1 <- ggplot(bg, aes(conservation)) +
    geom_histogram(data = subset(bg, conservation < 0.95), bins = 12,
                   boundary = 0.4, fill = "#C7D7E4", colour = "white", linewidth = 0.2) +
    annotate("rect", xmin = 0.95, xmax = 1.0, ymin = 0, ymax = 137,
             fill = "#E9A6A2", colour = "white", linewidth = 0.2) +
    annotate("text", x = 0.47, y = 148, label = "invariant\n(137)", colour = HOT,
             size = pt(5.0), hjust = 0, vjust = 1) +
    annotate("segment", x = 0.62, xend = 0.94, y = 146, yend = 138,
             colour = HOT, linewidth = 0.22) +
    scale_x_continuous(limits = c(0.38, 1.04), expand = c(0, 0)) +
    labs(x = "Conservation", y = "Residues")

  p2 <- ggplot(bg, aes(conservation, n_missense + 1, colour = grp)) +
    geom_point(size = 0.45, alpha = 0.55) +
    geom_point(data = mutate(d$hot, grp = "Hotspot codon"),
               aes(conservation, n_missense + 1, colour = grp),
               inherit.aes = FALSE, size = 1.5, stroke = 0.3,
               shape = 21, fill = HOT) +
    scale_colour_manual(values = c("Outside DBD" = "#ADB5BD", "DBD" = TEAL,
                                   "Hotspot codon" = HOT),
                        drop = FALSE, name = "Each point = one codon") +
    scale_y_log10(breaks = c(1, 10, 100), labels = c("1", "10", "100")) +
    annotate("text", x = 0.40, y = 80,
             label = sprintf("rho = %.2f\nP = %.0e", summ$spearman_conservation_vs_missense$rho,
                             summ$spearman_conservation_vs_missense$p),
             size = pt(5.0), colour = INK, hjust = 0, vjust = 1) +
    labs(x = "Conservation", y = "Missense per codon + 1") +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 5.2),
          legend.text = element_text(size = 5.0),
          legend.key.size = unit(2.4, "mm"))

  med <- median(bg$conservation)
  grp <- tibble(
    label = factor(c("Conserved\nhalf", "Variable\nhalf", "DBD", "Outside\nDBD"),
                   levels = c("Conserved\nhalf", "Variable\nhalf", "DBD", "Outside\nDBD")),
    n = c(sum(bg$n_missense[bg$conservation >= med]),
          sum(bg$n_missense[bg$conservation < med]),
          sum(bg$n_missense[bg$dbd]), sum(bg$n_missense[!bg$dbd])),
    size = c(sum(bg$conservation >= med), sum(bg$conservation < med),
             sum(bg$dbd), sum(!bg$dbd))) %>%
    mutate(dens = n / size)
  p3 <- ggplot(grp, aes(label, dens, fill = label)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.1f", dens)), vjust = -0.4,
              size = pt(5.0), colour = INK) +
    scale_fill_manual(values = c("#4A7FA8", "#ADB5BD", TEAL, "#A9BFD2"), guide = "none") +
    scale_y_continuous(limits = c(0, max(grp$dens) * 1.22), expand = c(0, 0)) +
    labs(x = NULL, y = "Missense per residue") +
    theme(axis.text.x = element_text(size = 4.9, lineheight = 0.9))

  fig <- (p1 | p2 | p3) + plot_layout(widths = c(1, 1, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
  fig <- fig +
    plot_annotation(title = "Conservation predicts where mutations accumulate",
                    tag_levels = "a")
  save_r(fig, "fig5_conservation_vs_mutation", 180, 72, multi = TRUE)
}

# ------------------------------------------------------------------ figure 6 --
figure6 <- function() {
  hs <- read_csv(file.path(TBL, "structure_hotspots.csv"), show_col_types = FALSE)
  img_a <- png::readPNG(file.path(ROOT, "results", "pymol", "fig6_view_a.png"))
  img_b <- png::readPNG(file.path(ROOT, "results", "pymol", "fig6_view_b.png"))
  grob_a <- grid::rasterGrob(img_a, interpolate = TRUE,
                             width = unit(54, "mm"), height = unit(54, "mm"))
  grob_b <- grid::rasterGrob(img_b, interpolate = TRUE,
                             width = unit(54, "mm"), height = unit(54, "mm"))

  role <- c("DNA contact" = HOT, "structural (buried)" = BLUE)
  hs2 <- hs %>%
    filter(structural_role %in% names(role)) %>%
    mutate(res = factor(paste0(wt, pos), levels = paste0(wt, pos)[order(rel_sasa)]))
  p3 <- ggplot(hs2, aes(rel_sasa, res, fill = structural_role)) +
    geom_col(width = 0.58) +
    geom_text(aes(label = sprintf("%.0f%%", rel_sasa)), hjust = -0.2,
              size = pt(5.0), colour = INK) +
    scale_fill_manual(values = role, name = "Hotspot class") +
    scale_x_continuous(limits = c(0, 112), expand = c(0, 0),
                       breaks = c(0, 25, 50, 75, 100)) +
    labs(x = "Relative solvent accessibility (%)", y = NULL) +
    theme(axis.text.y = element_text(face = "bold", size = 6.0),
          legend.position = "right", legend.direction = "vertical",
          legend.key.size = unit(2.4, "mm"),
          legend.title = element_text(size = 5.4),
          legend.text = element_text(size = 5.0))
  # the two PyMOL views are square rasters; a fixed-mm grid keeps them square
  # and aligned instead of letting patchwork rescale them
  draw <- function() {
    g3 <- ggplotGrob(p3)
    grid.newpage()
    pushViewport(viewport(width = unit(180, "mm"), height = unit(76, "mm")))
    pushViewport(viewport(layout = grid.layout(2, 1,
                                               heights = unit(c(8, 68), "mm"))))

    pushViewport(viewport(layout.pos.row = 1))
    grid.text("Hotspots at the DNA interface and in the buried core",
              x = unit(2, "mm"), just = "left",
              gp = gpar(fontsize = 7.5, fontface = "bold", col = INK))
    popViewport()

    pushViewport(viewport(layout.pos.row = 2))
    pushViewport(viewport(layout = grid.layout(
      1, 5, widths = unit(c(52, 3, 52, 3, 70), "mm"))))

    pushViewport(viewport(layout.pos.col = 1))
    grid.draw(grob_a)
    grid.text("a", x = unit(1, "mm"), y = unit(0.98, "npc"), just = c("left", "top"),
              gp = gpar(fontface = "bold", fontsize = 9, col = INK))
    popViewport()

    pushViewport(viewport(layout.pos.col = 3))
    grid.draw(grob_b)
    grid.text("b", x = unit(1, "mm"), y = unit(0.98, "npc"), just = c("left", "top"),
              gp = gpar(fontface = "bold", fontsize = 9, col = INK))
    popViewport()

    pushViewport(viewport(layout.pos.col = 5))
    grid.draw(g3)
    grid.text("c", x = unit(1, "mm"), y = unit(0.98, "npc"), just = c("left", "top"),
              gp = gpar(fontface = "bold", fontsize = 9, col = INK))
    popViewport()

    popViewport(2)
    popViewport(2)
  }

  ragg::agg_png(file.path(OUT, "fig6_structure.png"),
                width = 180 / 25.4, height = 76 / 25.4, units = "in",
                res = 600, background = "white")
  draw()
  dev.off()
  grDevices::pdf(file.path(PDF, "fig6_structure.pdf"),
                 width = 180 / 25.4, height = 76 / 25.4,
                 family = "Helvetica", useDingbats = FALSE)
  draw()
  dev.off()
  message("wrote fig6_structure")
}

args <- commandArgs(trailingOnly = TRUE)
figs <- if (length(args)) args else c("1", "2", "3", "4", "5", "6")
for (f in figs) {
  switch(f,
    "1" = figure1(), "2" = figure2(), "3" = figure3(),
    "4" = figure4(), "5" = figure5(), "6" = figure6())
}