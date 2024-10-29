source("03-model-info.r")   # load data
library("ggh4x")            # plotting help
library("Hmisc")            # capitalize
library("ggnewscale")       # fig 2
library("ggpubr")           # some plots are unpleasant

# ==============================================================================
# preamble
# ==============================================================================

if (!dir.exists("figures")) {
    dir.create("figures")
}

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------

# only aggregate trait stuff for main figs
main_filter <- function(x) {
    # easier to filter by what we don't want than what we do want
    dplyr::filter(x, 
        !grepl("size", tolower(resp)),
        !grepl("size", tolower(size_or_trait)),
        !grepl("size", tolower(PBLR_test_id)))
}

# only size stuff for supp figs
size_filter <- function(x) {
    dplyr::filter(x,
        !grepl("trait", tolower(resp)),
        !grepl("trait", tolower(size_or_trait)),
        !grepl("trait", tolower(PBLR_test_id)))
}

# auto tag facets with letters
tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf,
    y = Inf, hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(
        lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(
        data = tags, aes(x = x, y = y, label = label), ...,
        hjust = hjust, vjust = vjust, fontface = fontface, family = family,
        inherit.aes = FALSE)
}

# ------------------------------------------------------------------------------
# constants
# ------------------------------------------------------------------------------

# colour palette by treatment
treat_palette <- setNames(
    c("#999999", "#E69F00", "#56B4E9", "#009E73"),
    c("C", "T", "A", "AT")
)

# point shape by strain
pch_map <- c(
    # ciliates
    "Ciliates" = 15,
    # genera
    "Spiro" = 5, "Tetra" = 6, "Loxo" = 4,    
    # individual strains:
    "Spiro_C" = 0, "Spiro_D" = 5,
    "Tetra_1" = 2, "Tetra_2" = 6,
    "Loxo_1" = 3, "Loxo_2" = 4,
    # cyanobacteria
    "Cyanobacteria" = 16,
    # clades
    "VIII" = 9, "V" = 13,
    # individual strains:
    "VIII_2383" = 7, "VIII_2434" = 9,
    "V_2375" = 10, "V_2524" = 13
)

# point size by tax level
pointsize_map <- c(
    "system" = 4,
    "clade" = 2.5, 
    "strain" = 2
)

# alpha by tax level
alpha_map <- c(
    "system" = 0.5,
    "clade" = 0.75, 
    "strain" = 1
)

# ------------------------------------------------------------------------------
# plotting shorthands
# ------------------------------------------------------------------------------

# point types by strain etc.
pch_merged <- scale_shape_manual(values = pch_map)
pch_sys <- pch_map[c("Ciliates", "Cyanobacteria")]
pch_str <- pch_map[grep("_", names(pch_map))]

# removing a facet background
blank_strip <- element_rect(fill = "white", color = NA)

# double exp transform useful for fig 3 bc mainly interested in vals between .95 and 1
transform_exp2 <- scales::new_transform(
    "exp2", function(x) {exp(exp(x))}, function(x) {log(log(x))}
)

# ==============================================================================
# sorting out data to be ready to plot
# ==============================================================================

# ------------------------------------------------------------------------------
# main results
# ------------------------------------------------------------------------------

# merge main results and bootstrapping results
main_synth <- left_join(
    main_results %>%
        dplyr::filter(model_type == "augmented") %>%
        dplyr::select(system, clade, strain, treat, resp, pca_varexp, a_weights),
    data_synth %>% dplyr::select(strain, treat, resp, delta_error)) %>%
    ungroup

# summarise these by system
main_by_system <- main_synth %>%
    group_by(system, treat, resp) %>%
    summarise(across(pca_varexp:delta_error, mean)) %>%
    mutate(level = "system", id = system)

# summarise these by subsystem
main_by_clade <- main_synth %>%
    group_by(system, clade, treat, resp) %>%
    summarise(across(pca_varexp:delta_error, mean)) %>%
    mutate(level = "clade", id = clade)

# summarise these by strain i.e. no summarising
main_by_strain <- main_synth %>%
    mutate(level = "strain", id = strain)

# merge all these into a big data frame
main_long <- bind_rows(list(main_by_system, main_by_clade, main_by_strain)) %>%
    ungroup %>%
    mutate(stroke = ifelse(level == "clade", 1, 0.5))

# ------------------------------------------------------------------------------
# size results
# ------------------------------------------------------------------------------

# summarise these by system
size_by_system <- size_synth %>%
    group_by(system, treat, resp) %>%
    summarise(across(pca_varexp:delta_error, mean)) %>%
    mutate(level = "system", id = system)

# summarise these by subsystem
size_by_clade <- size_synth %>%
    group_by(system, clade, treat, resp) %>%
    summarise(across(pca_varexp:delta_error, mean)) %>%
    mutate(level = "clade", id = clade)

# summarise these by strain i.e. no summarising
size_by_strain <- size_synth %>%
    mutate(level = "strain", id = strain)

# merge all these into a big data frame
size_long <- bind_rows(list(size_by_system, size_by_clade, size_by_strain)) %>%
    ungroup %>%
    mutate(stroke = ifelse(level == "clade", 1, 0.5))

# ------------------------------------------------------------------------------
# principal component loadings
# ------------------------------------------------------------------------------

pc_data <- bind_rows(
    read.csv("data/cyanobacteria-pca-loadings.csv") %>% mutate(system = "Cyanobacteria"), 
    read.csv("data/ciliates-pca-loadings.csv") %>% mutate(system = "Ciliates")) %>%
    dplyr::select(system, strain, treat, trait, PC1) %>%
    mutate(treat = factor(treat, levels = c("C", "T", "A", "AT")),
        clade = str_split_i(strain, "_", 1))

# ==============================================================================
# Figure 1
# ==============================================================================

p <- ggplot(data_preds %>%
    # removing null models
    mutate(model_type = factor(model_type, levels = c("single", "augmented")))) +
    # elements
    aes(x = obs, y = pred, col = treat, pch = strain, lty = factor(aic_rank)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(alpha = 0.4, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
    # layout
    facet_nested(system + resp ~ model_type,
        scales = "free", independent = TRUE,
        strip = strip_nested(background_y = blank_strip),
        switch = "y",
        nest_line = element_line(),
        labeller = labeller(.default = capitalize)) +
    # legends
    pch_merged +
    scale_linetype(guide = "none") +
    scale_colour_manual(values = treat_palette) +
    labs(x = "Observed value",
        y = "Predicted value",
        pch = "Strain",
        col = "Treatment") +
    guides(pch = guide_legend(override.aes = list(size = 2, alpha = 1)),
        col = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
    # theming
    theme_bw() +
    theme(
        strip.placement.y = "outside",
	    strip.background.y = element_blank(),
        strip.text.y = element_text(colour = "black"))

p <- tag_facet(p)

ggsave("figures/fig-01.pdf", p, width = 6 * 0.833, height = 8 * 0.75)

# ==============================================================================
# Figure 2
# ==============================================================================

pch_traits <- c(
    "Aspect ratio" = 0,
    "Speed" = 1,
    "Linearity" = 2,
    "Size" = 3,
    "Chlorophyll" = 4,
    "Phycoerythrin" = 5,
    "Phycocyanin" = 6
)


# ------------------------------------------------------------------------------
# FIG S3
# ------------------------------------------------------------------------------

# Plot the PCA variance explained
p1 <- ggplot(data_synth) +
    aes(x = interaction(str_split_i(strain, "_", 2), clade), y = pca_varexp, fill = treat) +
    facet_grid2(. ~ system, scale = "free_x", space = "free_x") +
    scale_x_discrete(guide = "axis_nested") +
    geom_vline(xintercept = vline_pos, color = "grey92") +
    stat_summary(
        inherit.aes = FALSE,
        fun = mean,
        aes(x = 2, y = pca_varexp, yintercept = after_stat(y)),
        geom = "hline", lty = 2) +
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    labs(y = "Percentage variance explained\nby first principal component",
        x = "",
        fill = "Treatment") +
    theme_bw() +
    scale_fill_manual(values = treat_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
        labels = scales::percent) +
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), axis.title.x = element_blank())

p1 <- tag_facet(p1)

# ggsave("figures/fig-s03.pdf", p, width = 6, height = 2.5, device = "pdf")


p2 <- ggplot(pc_data) +
    aes(x = interaction(str_split_i(strain, "_", 2), clade),
        y = PC1, col = treat, pch = trait) +
    geom_hline(yintercept = 0) +
    geom_point(position = position_dodge(width = 0.5)) +
    facet_grid2(. ~ system, scale = "free_x", space = "free_x") +
    scale_x_discrete(guide = "axis_nested") +
    scale_color_manual(values = treat_palette, guide = "none") +
    scale_shape_manual(values = pch_traits) +
    theme_bw() +
    theme(legend.box = "vertical") +
    labs(x = "Strain", y = "PC1 loading", color = "Treatment", pch = "Trait") +
    guides(pch = guide_legend(order = 1))

p2 <- tag_facet(p2, tag_pool = letters[3:4])

p <- ggpubr::ggarrange(p1, p2, heights = c(0.4, 0.6), nrow = 2, align = "v")

ggsave("figures/fig-02.pdf", p, width = 7, height = 5.5, device = "pdf")

# ==============================================================================
# Figure 3
# ==============================================================================

p <- ggplot(main_long) +
    # elements
    aes(fill = treat, col = treat, x = a_weights, y = delta_error,
        size = level, alpha = level, stroke = stroke) +
    # strain
    geom_hline(yintercept = 0) +
    geom_point(aes(shape = id), data = main_long %>% dplyr::filter(level != "clade")) +
    scale_shape_manual(values = pch_map) + 
    # layout
    facet_grid2(resp ~ system,
        labeller = labeller(.default = capitalize)) +
    # legends +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = treat_palette) + 
    scale_fill_manual(values = treat_palette, guide = "none") +
    scale_size_manual(values = pointsize_map, guide = "none") + 
    scale_alpha_manual(values = alpha_map, guide = "none") + 
    labs(
        shape = "System & strain",
        x = "Augmented model Akaike weight",
        y = "Prediction accuracy difference",
        col = "Treatment") +
    theme_bw() +
    guides(shape = guide_legend(override.aes = list(size = c(rep(3, 2), rep(1.5, 10)))))

p <- tag_facet(p)

ggsave("figures/fig-03.pdf", p, width = 7, height = 5.5)
d
# ==============================================================================
# Figure 4
# ==============================================================================

bootstrapped_models <- readRDS("data/boot-test-results.RData") %>%
    dplyr::filter(size_or_trait == "trait")

boot_main_by_system <- bootstrapped_models %>%
    group_by(system, treat, resp, PBLR_test_id) %>%
    summarise(PBLR_test =  mean(PBLR_test)) %>%
    mutate(level = "system", id = system) %>%
    dplyr::select(-system)

boot_main_by_clade <- bootstrapped_models %>%
    mutate(clade = str_split_i(strain, "_", 1)) %>%
    group_by(system, clade, treat, resp, PBLR_test_id) %>%
    summarise(PBLR_test =  mean(PBLR_test)) %>%
    mutate(level = "clade", id = clade) %>%
    dplyr::select(-system, -clade)

boot_main_by_strain <- bootstrapped_models %>%
    mutate(level = "strain", id = strain)

boot_main_long <- bind_rows(
    list(boot_main_by_system, boot_main_by_clade, boot_main_by_strain)) %>%
    ungroup %>%
    mutate(stroke = ifelse(level == "clade", 1, 0.5)) %>%
    dplyr::select(system, clade, strain, treat, resp, PBLR_test_id, level, id, PBLR_test)

boot_mod_plot <- boot_main_long %>%
    mutate(
        PBLR_test_id = str_split_i(PBLR_test_id, " ~ ", 2),
        PBLR_test_id = factor(PBLR_test_id,
            levels = c("1", "density", "size", "trait"),
            labels = c("Constant", "Density-dependence", "Size-dependence", "Trait-dependence")),
        resp = factor(resp,
            levels = c("d_density", "d_size", "d_trait"),
            labels = c("Growth", "Size change", "Trait change"))) 

rect_df <- boot_mod_plot %>%
    dplyr::filter((resp == "Growth" & PBLR_test_id == "Density-dependence") |
        (resp == "Size change" & PBLR_test_id == "Size-dependence") |
        (resp == "Trait change" & PBLR_test_id == "Trait-dependence")) %>%
    dplyr::select(resp, PBLR_test_id) %>%
    distinct

p <- ggplot(boot_mod_plot %>% dplyr::filter(!is.na(strain))) +   
    # elements
    aes(x = system, y = PBLR_test, col = treat, pch = strain) +
    geom_point(position = position_dodge(width = 0.5)) +
    stat_summary(geom = "point", fun = "mean", aes(pch = system), size = 4, alpha = 0.5,
        position = position_dodge(width = 0.5)) +
    geom_rect(data = rect_df,
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, col = "black", fill = NA, inherit.aes = FALSE, lwd = 1.25) +
    geom_hline(yintercept = 0.95, lty = 2) +
    # layout
    facet_grid2(resp ~ PBLR_test_id,
        render_empty = FALSE) +
    # legend
    scale_y_continuous(transform = transform_exp2, breaks = c(0, 0.5, 0.95, 1)) +
    pch_merged +
    scale_color_manual(values = treat_palette) +
    labs(y = "Augmented model bootstrapped probability (1 - p)",
        x = "System", 
        pch = "Strain",
        color = "Treatment") +
    # theming
    theme_bw() +
    theme(legend.box = "horizontal") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)),
        pch = guide_legend(override.aes = list(alpha = 1)))

p <- tag_facet(p, hjust = -0.25)

ggsave("figures/fig-04.pdf", p, width = 8, height = 4)

# ==============================================================================
# Supplementary
# ==============================================================================

# some of the plots are separated by model system for size reasons
# switch between by uncommenting if needed:
# model_system <- "cilia"
# model_system <- "cyano"

# create bars to better separate strains
vline_pos <- (0.5):(length(unique(data_synth$strain)) - 0.5)

for (model_system in c("cilia", "cyano")) {

    if (model_system == "cilia") {
        trait_names <- c("size", "speed", "aspect_ratio", "linearity", "trait")
    } else {
        trait_names <- c("size", "chlorophyll", "phycocyanin", "phycoerythrin", "trait")
    }

    # load the raw data ----
    data <- as_tibble(read.csv("data/data_wide.csv")) %>%
        mutate(treat = factor(treat, levels = c("C", "T", "A", "AT"))) %>%
        dplyr::filter(grepl(model_system, tolower(system))) 

    # ------------------------------------------------------------------------------
    # Figures S1 and S2
    # ------------------------------------------------------------------------------

    # plot pops and traits ----
    p <- ggplot(data %>% mutate(density = log10(density)) %>%
            rename("density_log10" = density) %>%
            pivot_longer(all_of(c("density_log10", trait_names)))) +
        theme_bw() +
        scale_colour_manual(values = treat_palette) +
        aes(x = day, y = value, col = treat) +
        geom_point(size = 1) +
        facet_nested(name ~ species + str_split_i(strain, "_", 2), independent = "y", 
            nest_line = element_line(), scales = "free_y",
            switch = "y",
            labeller = labeller(.default = capitalize)) +
        labs(x = "Time (days)",
            y = NULL,
            col = "Treatment") +
        theme(
            legend.position = "bottom",
            strip.placement.y = "outside",
            strip.background.y = element_rect(fill = NA, color = NA))

    if (model_system == "cilia") {
        ggsave("figures/fig-s01.png", width = 8, height = 8, device = "png", dpi = 450)
    } else {
        ggsave("figures/fig-s02.png", width = 5, height = 7.5, device = "png", dpi = 450)
    }

    # --------------------------------------------------------------------------
    # Figures S4 and S5
    # --------------------------------------------------------------------------
    
    # plot pcgr vs. pop.
    p <- ggplot(data) +
        theme_bw() +
        scale_colour_manual(values = treat_palette) +
        aes(x = density, y = d_density, col = treat) +
        geom_hline(yintercept = 0) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
        facet_wrap2(.~strain, ncol = 2, scales = "free") +
        labs(x = "Density",
            y = "PCGR",
            col = "Treatment")

    p <- tag_facet(p)

    if (model_system == "cilia") {
        ggsave("figures/fig-s03.pdf", width = 5, height = 6, device = "pdf")
    } else {
        ggsave("figures/fig-s04.pdf", width = 5, height = 4, device = "pdf")
    }

    # --------------------------------------------------------------------------
    # Figures S6 and S7
    # --------------------------------------------------------------------------

    # plot dT vs. trait 
    p <- ggplot(data) +
        theme_bw() +
        scale_colour_manual(values = treat_palette) +
        aes(x = trait, y = d_trait, col = treat) +
        geom_hline(yintercept = 0) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
        facet_wrap(.~strain, ncol = 2, scales = "free") +
        labs(x = "Trait value",
            y = "Trait change",
            col = "Treatment")

    p <- tag_facet(p)

    if (model_system == "cilia") {
        ggsave("figures/fig-s05.pdf", width = 5, height = 6, device = "pdf")
    } else {
        ggsave("figures/fig-s06.pdf", width = 5, height = 4, device = "pdf")
    }

}

# ------------------------------------------------------------------------------
# Figure S8
# ------------------------------------------------------------------------------

# plot R squared
p <- ggplot(main_results) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = treat_palette) +
    aes(x = interaction(str_split_i(strain, "_", 2), clade), , y = R2, fill = treat) +
    scale_x_discrete(guide = "axis_nested") +
    geom_vline(xintercept = vline_pos, color = "grey92") +
    stat_summary(
        inherit.aes = FALSE,
        fun = mean,
        aes(x = 2, y = R2, yintercept = after_stat(y)),
        geom = "hline", lty = 2) +
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    geom_hline(yintercept = 0, col = "black") +
    facet_grid2(resp + model_type ~ system, scales = "free",
        labeller = labeller(.default = capitalize), space = "free_x",
        strip = strip_nested()) +
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
    labs(y = expression("R "^2),
        x = "Strain",
        fill = "Treatment")

p <- tag_facet(p)

ggsave("figures/fig-s07.pdf", p, width = 6, height = 6, device = "pdf")

# ------------------------------------------------------------------------------
# Figure S9
# ------------------------------------------------------------------------------

# form "pcgr ~ density"
p <- ggplot(size_preds %>%
    mutate(model_type = factor(model_type, levels = c("single", "augmented")))) +
    theme_bw() +#
    pch_merged +
    scale_colour_manual(values = treat_palette) +
    aes(x = obs, y = pred, col = treat, pch = strain, lty = factor(aic_rank)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(alpha = 0.4, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
    scale_linetype(guide = "none") +
    facet_nested(system + resp ~ model_type, scales = "free", independent = TRUE,
        strip = strip_nested(background_y = blank_strip), switch = c("y"),
        nest_line = element_line(),
        labeller = labeller(.default = capitalize)) +
    theme(strip.placement.y = "outside",
	    strip.background.y = element_blank(),
        strip.text.y = element_text(colour = "black")) +
    labs(x = "Observed value",
        y = "Predicted value",
        pch = "Strain",
        col = "Treatment") +
    theme(strip.placement.y = "outside",
	    strip.background.y = element_blank(),
        strip.text.y = element_text(colour = "black")) +
    guides(pch = guide_legend(override.aes = list(size = 2, alpha = 1)),
        col = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

p <- tag_facet(p)

ggsave("figures/fig-s08.pdf", p, width = 6 * 0.833, height = 8 * 0.75)

# ------------------------------------------------------------------------------
# Figure S10
# ------------------------------------------------------------------------------

p <- ggplot(size_long %>% dplyr::filter(level != "clade")) +
    # elements
    aes(fill = treat, col = treat, x = a_weights_augmented, y = delta_error,
        size = level, alpha = level, stroke = stroke) +
    # strain
    geom_hline(yintercept = 0) +
    geom_point(aes(shape = id)) +
    # layout
    facet_grid2(resp ~ system,
        labeller = labeller(.default = capitalize)) +
    # legends
    scale_shape_manual(values = pch_map) + 
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = treat_palette) + 
    scale_fill_manual(values = treat_palette, guide = "none") +
    scale_size_manual(values = pointsize_map, guide = "none") + 
    scale_alpha_manual(values = alpha_map, guide = "none") + 
    labs(
        shape = "System & strain",
        x = "Augmented model Akaike weight",
        y = "Prediction accuracy difference",
        col = "Treatment") +
    theme_bw() +
    guides(shape = guide_legend(override.aes = list(size = c(rep(3, 2), rep(1.5, 10)))))

p <- tag_facet(p)

ggsave("figures/fig-s09.pdf", p, width = 7, height = 5.5, device = "pdf")

# ------------------------------------------------------------------------------
# Figure S11
# ------------------------------------------------------------------------------

bootstrapped_models <- readRDS("data/boot-test-results.RData") %>%
    dplyr::filter(size_or_trait == "size")

boot_size_by_system <- bootstrapped_models %>%
    group_by(system, treat, resp, PBLR_test_id) %>%
    summarise(PBLR_test =  mean(PBLR_test)) %>%
    mutate(level = "system", id = system) %>%
    dplyr::select(-system)

boot_size_by_clade <- bootstrapped_models %>%
    mutate(clade = str_split_i(strain, "_", 1)) %>%
    group_by(system, clade, treat, resp, PBLR_test_id) %>%
    summarise(PBLR_test =  mean(PBLR_test)) %>%
    mutate(level = "clade", id = clade) %>%
    dplyr::select(-system, -clade)

boot_size_by_strain <- bootstrapped_models %>%
    mutate(level = "strain", id = strain)

boot_size_long <- bind_rows(
    list(boot_size_by_system, boot_size_by_clade, boot_size_by_strain)) %>%
    ungroup %>%
    mutate(stroke = ifelse(level == "clade", 1, 0.5)) %>%
    dplyr::select(system, clade, strain, treat, resp, PBLR_test_id, level, id, PBLR_test)

boot_mod_plot <- boot_size_long %>%
    mutate(
        PBLR_test_id = str_split_i(PBLR_test_id, " ~ ", 2),
        PBLR_test_id = factor(PBLR_test_id,
            levels = c("1", "density", "size", "trait"),
            labels = c("Constant", "Density-dependence", "Size-dependence", "Trait-dependence")),
        resp = factor(resp,
            levels = c("d_density", "d_size", "d_trait"),
            labels = c("Growth", "Size change", "Trait change"))) 

rect_df <- boot_mod_plot %>%
    dplyr::filter((resp == "Growth" & PBLR_test_id == "Density-dependence") |
        (resp == "Size change" & PBLR_test_id == "Size-dependence") |
        (resp == "Trait change" & PBLR_test_id == "Trait-dependence")) %>%
    dplyr::select(resp, PBLR_test_id) %>%
    distinct

p <- ggplot(boot_mod_plot %>% dplyr::filter(!is.na(strain))) +   
    aes(x = system, y = PBLR_test, col = treat, pch = strain) +
    geom_point(position = position_dodge(width = 0.5)) +
    stat_summary(geom = "point", fun = "mean", aes(pch = system), size = 4, alpha = 0.5,
        position = position_dodge(width = 0.5)) +
    geom_rect(data = rect_df,
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, col = "black", fill = NA, inherit.aes = FALSE, lwd = 1.25) +
    geom_hline(yintercept = 0.95, lty = 2) +
    facet_grid2(resp ~ PBLR_test_id, render_empty = FALSE) +
    scale_y_continuous(transform = transform_exp2, breaks = c(0, 0.5, 0.95, 1)) +
    theme_bw() +
    pch_merged +
    scale_color_manual(values = treat_palette) +
    labs(y = "Augmented model bootstrapped probability (1 - p)",
        x = "System", 
        pch = "Strain",
        color = "Treatment") +
    theme(legend.box = "horizontal") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)),
        pch = guide_legend(override.aes = list(alpha = 1)))

p <- tag_facet(p, hjust = -0.25)

ggsave("figures/fig-s10.pdf", p, width = 8, height = 4)
