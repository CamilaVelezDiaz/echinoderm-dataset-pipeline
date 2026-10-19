# ==============================================================================
# RAREFACTION AND RICHNESS ESTIMATION
# ==============================================================================
# Camila Velez | JCU Marine Biology MSc | Thesis Chapter — Analysis
#
# Addresses: given current sampling effort, how many species are likely to
# exist in each depth zone, and is sampling approaching completeness?
#
# Two things are produced from the SAME underlying analysis (not two
# separate analyses - "rarefied richness" is simply the standardized
# comparison you read off the rarefaction/extrapolation curve):
#   1. Sample-based rarefaction/extrapolation curves + Chao1 asymptotic
#      richness estimator, per depth zone (via the iNEXT package)
#   2. Rarefied richness at a common, standardized sample size - lets zones
#      with very different amounts of collecting effort be compared fairly,
#      rather than the Continental Shelf trivially "winning" on raw
#      richness just because it has vastly more records.
#
# Data type: ABUNDANCE. For each depth zone, each species' "abundance" is
# its number of records in that zone - this is what iNEXT calls
# individual-based abundance data (as opposed to incidence/presence-absence
# across multiple sampling sites).
#
# BOOTSTRAP REPLICATES: nboot = 999 throughout (increased from an earlier
# nboot = 100). iNEXT's own package default is 50, and a closely related
# package by the same authors (iNEXT.beta3D) suggests 100-200 as the
# "more accurate" tier - neither documents 999/1000 specifically. However,
# 999 matches the general bootstrap-methodology convention widely used
# elsewhere (phylogenetics, econometrics), so it was adopted here as the
# more conservative choice. Point estimates (Chao1, rarefied richness) are
# unaffected by nboot - only the width of the bootstrap CI/SE changes -
# and results were confirmed stable between nboot = 100 and nboot = 999
# (same point estimates, same pairwise significance conclusions, same
# compact letter groupings).
#
# CHANGE LOG (this version):
#   - Part 2 (rarefied richness at a standardized sample size): the
#     reference sample size is now chosen only among zones with at least
#     MIN_RECORDS_FOR_REFERENCE records, rather than always taking the
#     single smallest zone's total. Previously, a zone with an extremely
#     small species-level sample (e.g. Abyssal, once it had n=3 records
#     post-reclassification) would force EVERY zone down to that tiny
#     common sample size, producing a ceiling effect: well-sampled zones
#     trivially return rarefied_richness ~= reference_size regardless of
#     their true diversity, and the resulting near-zero standard errors
#     made trivial numeric differences appear spuriously "significant" in
#     the pairwise test. Such zones are now excluded from SETTING the
#     reference size, but are NOT excluded from the comparison itself -
#     they still receive an estimateD() value, now via extrapolation to
#     the larger, more informative reference size set by the remaining
#     zones, with correspondingly wider (and more honest) uncertainty.
#
# Input:  echino_wide.csv
# Output: rarefaction_extrapolation_curves.png
#         chao1_richness_by_zone.csv
#         rarefied_richness_comparison.csv
# ==============================================================================
# =============================================================================
# SETUP
# =============================================================================
setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")
# =============================================================================
# Library
# =============================================================================
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(ggplot2)
library(iNEXT)
library(patchwork)   # for inset panels on the rarefaction/extrapolation curves
# Fixed seed for REPRODUCIBLE bootstrap confidence intervals. iNEXT's CI
# calculation uses random resampling (nboot below) - without a fixed seed,
# every re-run of this script produces slightly different CI bounds (point
# estimates stay the same, only the CIs shift). Fixing the seed means this
# script always produces the exact same numbers, which matters for a
# result going into a thesis.
set.seed(42)
NBOOT <- 999  # single source of truth - change once here if this needs revisiting
# ------------------------------------------------------------------------------
# FIXED ZONE COLOUR PALETTE (used by both the species-level and genus-level
# rarefaction/extrapolation figures, so the two are visually consistent).
# Default ggplot2 hue colours put Continental Shelf (red) and Abyssal (purple)
# too close in perceived tone to distinguish at a glance, especially
# once the CI ribbon transparency is applied - this is a colorblind-friendly
# ColorBrewer "Dark2" palette instead, with each zone's colour chosen to be
# maximally separated from the others.
# ------------------------------------------------------------------------------
zone_colors <- c(
  "Continental Shelf" = "#1B9E77",  # teal-green
  "Upper Slope"        = "#D95F02",  # orange
  "Lower Slope"         = "#7570B3",  # blue-purple
  "Abyssal"     = "#E7298A"   # magenta/pink
)
# Depth-range legend labels (name on top, range below), matching the
# convention used across every other depth-zone figure in this thesis
# (Taxonomic Completeness Figure 1, TIa Figures 4-5, TRi by depth zone)
zone_labels_with_range <- c(
  "Continental Shelf" = "Continental Shelf\n(0 - 200 m)",
  "Upper Slope"        = "Upper Slope\n(200 - 1000 m)",
  "Lower Slope"         = "Lower Slope\n(1000 - 4000 m)",
  "Abyssal"     = "Abyssal\n(>4000 m)"
)
# ------------------------------------------------------------------------------
# LOAD DATA (defensive checks matching every other script in this project)
# ------------------------------------------------------------------------------
echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)
if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c("Species"="Species","Subspecies"="Species","Genus"="Genus","Subgenus"="Genus",
                "Family"="Family","Order"="Order","Class"="Class","Phylum"="Phylum")
  echino_wide <- echino_wide %>% mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}
if (!"crosses_zone_boundary" %in% names(echino_wide)) {
  echino_wide <- echino_wide %>% mutate(crosses_zone_boundary = str_starts(coalesce(depth_zone, ""), "Spans"))
}
zone_order <- c("Continental Shelf", "Upper Slope", "Lower Slope", "Abyssal")
# ------------------------------------------------------------------------------
# BUILD ABUNDANCE DATA: species x zone record counts
# ------------------------------------------------------------------------------
species_by_zone <- echino_wide %>%
  filter(rank6 == "Species", !is.na(accepted_name_final), !crosses_zone_boundary,
         depth_zone %in% zone_order) %>%
  count(depth_zone, accepted_name_final, name = "n_records")
cat("=== Records feeding rarefaction (species-level, non-boundary-crossing) ===\n")
print(species_by_zone %>% group_by(depth_zone) %>%
        summarise(n_species = n_distinct(accepted_name_final), n_records = sum(n_records), .groups = "drop") %>%
        mutate(depth_zone = factor(depth_zone, levels = zone_order)) %>% arrange(depth_zone))
abundance_list <- species_by_zone %>%
  filter(depth_zone %in% zone_order) %>%
  split(.$depth_zone) %>%
  lapply(function(df) {
    v <- df$n_records
    names(v) <- df$accepted_name_final
    sort(v, decreasing = TRUE)
  })
abundance_list <- abundance_list[intersect(zone_order, names(abundance_list))]
n_species_per_zone <- sapply(abundance_list, length)
cat("\nSpecies count per zone (feeding iNEXT):\n")
print(n_species_per_zone)
missing_zones <- setdiff(zone_order, names(abundance_list))
if (length(missing_zones) > 0) {
  cat("\n\u26A0\uFE0F  ZONE(S) EXCLUDED ENTIRELY:", paste(missing_zones, collapse = ", "), "\n")
  cat("These zones have ZERO species-level records and cannot contribute to\n")
  cat("rarefaction/Chao1 at all - there is no species identity to rarefy. This is a\n")
  cat("genuine data limitation to state explicitly in Results/Discussion, not a\n")
  cat("script error: rarefaction and richness estimation for this zone are simply\n")
  cat("not possible with the current data, and would require the deeper-water\n")
  cat("collections to be identified at least to species level before any estimate\n")
  cat("of true richness there could be attempted.\n")
}
too_small <- names(n_species_per_zone)[n_species_per_zone < 3]
if (length(too_small) > 0) {
  cat("\n\u26A0\uFE0F  WARNING:", paste(too_small, collapse = ", "),
      "have fewer than 3 species - Chao1/rarefaction results for these zones\n")
  cat("will be unstable or undefined, and should be interpreted with extreme caution\n")
  cat("(or potentially excluded from the comparison and noted as a data limitation).\n")
}
# ==============================================================================
# PART 1: RAREFACTION / EXTRAPOLATION CURVES + CHAO1
# ==============================================================================
cat("\n=== Running iNEXT (this may take a moment - nboot =", NBOOT, ") ===\n")
inext_result <- iNEXT(abundance_list, q = 0, datatype = "abundance",
                      endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = NBOOT)
chao1_table <- inext_result$AsyEst %>%
  filter(Diversity == "Species richness") %>%
  transmute(
    depth_zone = Assemblage,
    observed_richness = Observed,
    chao1_estimate = round(Estimator, 1),
    chao1_ci_lower = round(LCL, 1),
    chao1_ci_upper = round(UCL, 1)
  ) %>%
  mutate(depth_zone = factor(depth_zone, levels = zone_order)) %>%
  arrange(depth_zone)
cat("\n── Chao1 asymptotic species richness by depth zone ──\n")
cat("(observed = species actually recorded; chao1_estimate = predicted TRUE\n")
cat("richness accounting for likely undetected species; a wide gap between\n")
cat("observed and estimate suggests sampling is far from complete)\n\n")
print(chao1_table)
write_csv(chao1_table, "chao1_richness_by_zone.csv")
# --- Rarefaction/extrapolation curve plot ---
species_zone_order <- zone_order[zone_order %in% names(abundance_list)]
p_curves <- ggiNEXT(inext_result, type = 1, se = TRUE) +
  labs(x = "Number of records (sample size)", y = "Species richness",
       title = "Sample-based rarefaction and extrapolation by depth zone",
       subtitle = "Solid = rarefaction (interpolated); dashed = extrapolation; shading = 95% CI") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = zone_colors, limits = species_zone_order,
                      labels = zone_labels_with_range[species_zone_order]) +
  scale_fill_manual(values = zone_colors, limits = species_zone_order,
                    labels = zone_labels_with_range[species_zone_order]) +
  scale_shape_discrete(limits = species_zone_order,
                       labels = zone_labels_with_range[species_zone_order])
ggsave("rarefaction_extrapolation_curves.png", p_curves, width = 9, height = 6, dpi = 300)
cat("\nFigure saved: rarefaction_extrapolation_curves.png\n")


# ------------------------------------------------------------------------------
# INSET: zoomed view of the smaller zones (Continental Shelf excluded)
# ------------------------------------------------------------------------------
species_small_zones <- setdiff(names(abundance_list), "Continental Shelf")
small_zones_abundance <- abundance_list[species_small_zones]
if (length(small_zones_abundance) >= 2) {
  inext_small <- iNEXT(small_zones_abundance, q = 0, datatype = "abundance",
                       endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = NBOOT)
  
  small_zone_order <- species_zone_order[species_zone_order %in% names(small_zones_abundance)]
  
  p_inset <- ggiNEXT(inext_small, type = 1, se = TRUE) +
    scale_colour_manual(values = zone_colors, limits = small_zone_order,
                        labels = zone_labels_with_range[small_zone_order]) +
    scale_fill_manual(values = zone_colors, limits = small_zone_order,
                      labels = zone_labels_with_range[small_zone_order]) +
    scale_shape_discrete(limits = small_zone_order,
                         labels = zone_labels_with_range[small_zone_order]) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = "grey40"),
          plot.margin = margin(4, 8, 4, 4)) +
    labs(x = NULL, y = NULL)
  
  p_curves_with_inset <- p_curves +
    inset_element(p_inset, left = 0.50, bottom = 0.06, right = 0.97, top = 0.48)
  
  ggsave("rarefaction_extrapolation_curves_inset.png", p_curves_with_inset,
         width = 9, height = 6, dpi = 300)
  cat("Figure saved: rarefaction_extrapolation_curves_inset.png (with zoomed inset)\n")
} else {
  cat("Fewer than 2 smaller zones have usable species-level data - inset skipped.\n")
}

# Vector version for poster
ggsave("rarefaction_extrapolation_curves_inset.pdf", p_curves_with_inset,
       width = 9, height = 6, dpi = 300,
       device = cairo_pdf)
cat("rarefaction_extrapolation_curves_inset.pdf (vector)\n")

# ==============================================================================
# PART 2: RAREFIED RICHNESS AT A STANDARDIZED SAMPLE SIZE
# ==============================================================================
# The reference size is normally set by the smallest zone's total record
# count, so every zone can be fairly compared at a common effort level.
# However, a zone with an extremely small species-level sample (e.g.
# Abyssal, n=3) should NOT be allowed to set that reference: forcing every
# other zone down to n=3 produces a ceiling effect where any well-sampled
# zone trivially returns rarefied_richness ~ 3 regardless of its true
# diversity, and the resulting near-zero standard errors can make trivial
# numeric differences appear spuriously "significant" in the pairwise test
# below. MIN_RECORDS_FOR_REFERENCE excludes such zones from SETTING the
# reference size, without excluding them from the comparison itself - they
# still receive an estimateD() value, now via extrapolation to the larger,
# more informative reference size set by the remaining zones.
# ==============================================================================
MIN_RECORDS_FOR_REFERENCE <- 20
zones_for_reference <- names(abundance_list)[sapply(abundance_list, sum) >= MIN_RECORDS_FOR_REFERENCE]

if (length(zones_for_reference) < length(abundance_list)) {
  excluded_from_ref <- setdiff(names(abundance_list), zones_for_reference)
  cat(sprintf("\nNOTE: %s excluded from SETTING the standardized reference size ",
              paste(excluded_from_ref, collapse = ", ")))
  cat(sprintf("(fewer than %d records) - including it would force every zone down\n", MIN_RECORDS_FOR_REFERENCE))
  cat("to an uninformatively small common sample size. It is still INCLUDED in\n")
  cat("the standardized comparison itself, now via extrapolation to the larger\n")
  cat("reference size set by the remaining zones.\n\n")
}

reference_size <- min(sapply(abundance_list[zones_for_reference], sum))
cat(sprintf("\nStandardizing all zones to a common sample size of %d records\n", reference_size))
cat("(the smallest zone's total record count, among zones with a sufficient\n")
cat("species-level sample) for fair richness comparison.\n\n")

rarefied_estimates <- estimateD(abundance_list, datatype = "abundance", base = "size",
                                level = reference_size, q = 0, conf = 0.95, nboot = NBOOT)
rarefied_richness <- rarefied_estimates %>%
  transmute(
    depth_zone = Assemblage,
    reference_size = m,
    method = Method,
    rarefied_richness = round(qD, 1),
    ci_lower = round(qD.LCL, 1),
    ci_upper = round(qD.UCL, 1)
  ) %>%
  mutate(depth_zone = factor(depth_zone, levels = zone_order)) %>%
  arrange(depth_zone)
cat("── Rarefied species richness, standardized to", reference_size, "records per zone ──\n")
print(rarefied_richness)
write_csv(rarefied_richness, "rarefied_richness_comparison.csv")
cat("\nNOTE: 'method' indicates whether this standardized value was\n")
cat("interpolated (rarefaction, within the zone's own observed data range)\n")
cat("or extrapolated (beyond it) - extrapolated values carry more\n")
cat("uncertainty and should be flagged as such when reported.\n")
# ------------------------------------------------------------------------------
# PAIRWISE SIGNIFICANCE TEST
# ------------------------------------------------------------------------------
pairwise_rarefied_test <- function(richness_df, value_col = "rarefied_richness",
                                   lower_col = "ci_lower", upper_col = "ci_upper",
                                   group_col = "depth_zone") {
  df <- richness_df
  df$se <- (df[[upper_col]] - df[[lower_col]]) / (2 * 1.96)
  groups <- as.character(df[[group_col]])
  n <- length(groups)
  
  pairs <- combn(n, 2)
  results <- vector("list", ncol(pairs))
  for (i in seq_len(ncol(pairs))) {
    a <- pairs[1, i]; b <- pairs[2, i]
    diff <- df[[value_col]][a] - df[[value_col]][b]
    se_diff <- sqrt(df$se[a]^2 + df$se[b]^2)
    z <- diff / se_diff
    p_raw <- 2 * pnorm(-abs(z))
    results[[i]] <- tibble(group_a = groups[a], group_b = groups[b],
                           diff = round(diff, 1), z = round(z, 2), p_raw = p_raw)
  }
  results <- bind_rows(results) %>% mutate(p_holm = p.adjust(p_raw, method = "holm"))
  
  if (requireNamespace("multcompView", quietly = TRUE)) {
    pvec <- setNames(results$p_holm, paste(results$group_a, results$group_b, sep = "-"))
    raw_letters <- multcompView::multcompLetters(pvec)$Letters
    sorted_groups <- groups[order(-df[[value_col]])]
    new_letters <- setNames(letters[seq_len(n)], sorted_groups)
    letter_table <- tibble(!!group_col := names(raw_letters),
                           letter = new_letters[names(raw_letters)]) %>%
      left_join(df %>% dplyr::select(all_of(group_col), all_of(value_col)), by = group_col) %>%
      arrange(desc(.data[[value_col]]))
  } else {
    letter_table <- NULL
  }
  
  list(pairwise = results, letters = letter_table)
}
cat("\n── Pairwise significance test: species-level rarefied richness ──\n")
species_pairwise <- pairwise_rarefied_test(rarefied_richness)
print(species_pairwise$pairwise)
if (!is.null(species_pairwise$letters)) {
  cat("\nCompact letter grouping (a = highest richness; zones sharing a letter\n")
  cat("are NOT significantly different):\n")
  print(species_pairwise$letters)
}
write_csv(species_pairwise$pairwise, "rarefied_richness_pairwise_species.csv")
# ==============================================================================
# PART 3 (SUPPLEMENTARY): GENUS-LEVEL RAREFACTION, ALL 4 ZONES
# ==============================================================================
# Included to cover any zone whose species-level sample is too small to
# rarefy or compare meaningfully (see Part 2's reference-size note above).
# Genus-level richness is NOT directly comparable to species-level
# richness - reported separately and consistently across all 4 zones for
# that reason.
# ==============================================================================
cat("\n\n=== PART 3 (SUPPLEMENTARY): genus-level rarefaction, all 4 zones ===\n")
cat("Included to give every zone - including those with a very small\n")
cat("species-level sample - a sturdier basis for comparison. Genus-level\n")
cat("richness is NOT directly comparable to species-level richness - reported\n")
cat("separately and consistently across all 4 zones for that reason.\n\n")
genus_by_zone <- echino_wide %>%
  filter(rank6 %in% c("Species", "Genus"), !is.na(accepted_name_final),
         !crosses_zone_boundary, depth_zone %in% zone_order) %>%
  mutate(genus = word(accepted_name_final, 1)) %>%
  filter(!is.na(genus), genus != "") %>%
  count(depth_zone, genus, name = "n_records")
cat("Records feeding genus-level rarefaction (species- and genus-level,\n")
cat("non-boundary-crossing):\n")
print(genus_by_zone %>% group_by(depth_zone) %>%
        summarise(n_genera = n_distinct(genus), n_records = sum(n_records), .groups = "drop") %>%
        mutate(depth_zone = factor(depth_zone, levels = zone_order)) %>% arrange(depth_zone))
genus_abundance_list <- genus_by_zone %>%
  split(.$depth_zone) %>%
  lapply(function(df) {
    v <- df$n_records
    names(v) <- df$genus
    sort(v, decreasing = TRUE)
  })
genus_abundance_list <- genus_abundance_list[intersect(zone_order, names(genus_abundance_list))]
n_genera_per_zone <- sapply(genus_abundance_list, length)
cat("\nGenus count per zone (feeding iNEXT):\n")
print(n_genera_per_zone)
still_missing <- setdiff(zone_order, names(genus_abundance_list))
if (length(still_missing) > 0) {
  cat("\n\u26A0\uFE0F  Even at genus level,", paste(still_missing, collapse = ", "),
      "has zero usable records.\n")
  cat("No richness estimate of any kind is possible for this zone with current data.\n")
}
too_small_genus <- names(n_genera_per_zone)[n_genera_per_zone < 3]
if (length(too_small_genus) > 0) {
  cat("\n\u26A0\uFE0F  WARNING:", paste(too_small_genus, collapse = ", "),
      "have fewer than 3 genera - results unstable, interpret with caution.\n")
}
if (length(genus_abundance_list) >= 2) {
  inext_genus <- iNEXT(genus_abundance_list, q = 0, datatype = "abundance",
                       endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = NBOOT)
  
  chao1_genus_table <- inext_genus$AsyEst %>%
    filter(Diversity == "Species richness") %>%
    transmute(
      depth_zone = Assemblage,
      observed_genera = Observed,
      chao1_genus_estimate = round(Estimator, 1),
      chao1_ci_lower = round(LCL, 1),
      chao1_ci_upper = round(UCL, 1)
    ) %>%
    mutate(depth_zone = factor(depth_zone, levels = zone_order)) %>%
    arrange(depth_zone)
  
  cat("\n── Chao1 asymptotic GENUS richness by depth zone (all 4 zones) ──\n")
  print(chao1_genus_table)
  write_csv(chao1_genus_table, "chao1_genus_richness_by_zone.csv")
  
  cat("\n── Pairwise significance test: genus-level Chao1 estimates (all 4 zones) ──\n")
  cat("NOTE: a zone with a very small sample size will have a wide CI, so any\n")
  cat("comparison involving it may show NO significant difference from other\n")
  cat("zones simply because the CI is too wide to detect one - this itself is\n")
  cat("informative (confirms the zone is too sparsely sampled to draw firm\n")
  cat("conclusions), not a null result to be read as 'genuinely similar diversity'.\n\n")
  genus_pairwise <- pairwise_rarefied_test(chao1_genus_table,
                                           value_col = "chao1_genus_estimate",
                                           lower_col = "chao1_ci_lower",
                                           upper_col = "chao1_ci_upper")
  print(genus_pairwise$pairwise)
  if (!is.null(genus_pairwise$letters)) {
    cat("\nCompact letter grouping (a = highest estimated richness):\n")
    print(genus_pairwise$letters)
  }
  write_csv(genus_pairwise$pairwise, "chao1_genus_pairwise.csv")
  
  genus_zone_order <- zone_order[zone_order %in% names(genus_abundance_list)]
  
  p_genus_curves <- ggiNEXT(inext_genus, type = 1, se = TRUE) +
    labs(x = "Number of records (sample size)", y = "Genus richness",
         title = "Genus-level rarefaction and extrapolation by depth zone (supplementary)",
         subtitle = "Included for zones whose species-level sample is too small to compare directly") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom") +
    scale_colour_manual(values = zone_colors, limits = genus_zone_order,
                        labels = zone_labels_with_range[genus_zone_order]) +
    scale_fill_manual(values = zone_colors, limits = genus_zone_order,
                      labels = zone_labels_with_range[genus_zone_order]) +
    scale_shape_discrete(limits = genus_zone_order,
                         labels = zone_labels_with_range[genus_zone_order])
  
  ggsave("genus_rarefaction_curves_supplementary.png", p_genus_curves, width = 9, height = 6, dpi = 300)
  cat("\nFigure saved: genus_rarefaction_curves_supplementary.png\n")
  
  # ---- INSET: same rationale as the species-level plot above ----
  genus_small_zones <- setdiff(names(genus_abundance_list), "Continental Shelf")
  small_genus_abundance <- genus_abundance_list[genus_small_zones]
  
  if (length(small_genus_abundance) >= 2) {
    inext_genus_small <- iNEXT(small_genus_abundance, q = 0, datatype = "abundance",
                               endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = NBOOT)
    
    small_genus_zone_order <- genus_zone_order[genus_zone_order %in% names(small_genus_abundance)]
    
    p_genus_inset <- ggiNEXT(inext_genus_small, type = 1, se = TRUE) +
      scale_colour_manual(values = zone_colors, limits = small_genus_zone_order,
                          labels = zone_labels_with_range[small_genus_zone_order]) +
      scale_fill_manual(values = zone_colors, limits = small_genus_zone_order,
                        labels = zone_labels_with_range[small_genus_zone_order]) +
      scale_shape_discrete(limits = small_genus_zone_order,
                           labels = zone_labels_with_range[small_genus_zone_order]) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_blank(),
            plot.background = element_rect(fill = "white", colour = "grey40"),
            plot.margin = margin(4, 8, 4, 4)) +
      labs(x = NULL, y = NULL)
    
    p_genus_curves_with_inset <- p_genus_curves +
      inset_element(p_genus_inset, left = 0.50, bottom = 0.06, right = 0.97, top = 0.48)
    
    ggsave("genus_rarefaction_curves_supplementary_inset.png", p_genus_curves_with_inset,
           width = 9, height = 6, dpi = 300)
    cat("Figure saved: genus_rarefaction_curves_supplementary_inset.png (with zoomed inset)\n")
  } else {
    cat("Fewer than 2 smaller zones have usable genus-level data - inset skipped.\n")
  }
} else {
  cat("\nFewer than 2 zones have usable genus-level data - iNEXT requires at\n")
  cat("least 2 assemblages to run. Genus-level comparison not possible.\n")
}
cat("\n=== Rarefaction and richness estimation complete ===\n")
if (length(missing_zones) > 0) {
  cat("NOTE: species-level results (Parts 1-2) cover", length(abundance_list), "of",
      length(zone_order), "depth zones -", paste(missing_zones, collapse = ", "),
      "excluded due to zero species-level records (see warning above).\n")
  cat("The genus-level supplementary analysis (Part 3) covers all zones with\n")
  cat("any usable data, including", paste(missing_zones, collapse = ", "), ".\n")
}
cat("Outputs: rarefaction_extrapolation_curves.png, rarefaction_extrapolation_curves_inset.png,\n")
cat("chao1_richness_by_zone.csv, rarefied_richness_comparison.csv,\n")
cat("rarefied_richness_pairwise_species.csv, chao1_genus_richness_by_zone.csv,\n")
cat("chao1_genus_pairwise.csv, genus_rarefaction_curves_supplementary.png,\n")
cat("genus_rarefaction_curves_supplementary_inset.png\n")