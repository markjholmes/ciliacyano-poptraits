library("tidyverse") # general

# functions to compute trait change and growth rate
pcgr <- function(density, day) {
    log(lead(density, 1) / density) / (lead(day, 1) - day)
}

trait_change <- function(x, day) {
    (lead(x, 1) - x) / (lead(day, 1) - day)
}

# ==============================================================================
# Import and make the ciliate data uniform
# ==============================================================================

# genera
cilia_sp_levels <- c("Spiro", "Tetra", "Loxo", "Para")
cilia_sp_labels <- c("Spirostomum", "Tetrahymena", "Loxocephalus", "Paramecium")

# traits
cilia_trait_levels <- c("size", "speed", "aspect_ratio", "linearity", "trait")
cilia_trait_labels <- c("Size", "Speed", "Aspect ratio", "Linearity", "Trait")

# actual data
data_cilia <- read_csv("data/ciliates.csv") %>% %>%
    # PCA
    group_by(strain, treat) %>%
    mutate(day = as.numeric(day - first(day))) %>%
    nest %>%
    rowwise %>%
    mutate(
        pca_data = list(dplyr::select(data, aspect_ratio, size, speed, linearity)),
        pca = list(prcomp(pca_data, center = TRUE, scale = TRUE)),
        data = list(cbind(data, pca1 = pca$x[, 1])),
        pca_varexp = (pca$sdev^2)[1] / sum(pca$sdev^2)) %>%
    unnest(data) %>%
    # epxand pca and compute dT and clean up stuff
    rename(trait = pca1) %>%
    mutate(
        d_density = pcgr(density, day),
        across(all_of(cilia_trait_levels), ~ trait_change(.x, day),
            .names = "d_{.col}")) %>%
    ungroup %>%
    dplyr::filter(!is.na(d_density)) %>%
    # link to species
    mutate(species = gsub(".{2}$", "", strain))

# ------------------------------------------------------------------------------
# save the PCA loadings
# ------------------------------------------------------------------------------

cilia_pca <- data_cilia %>%
    dplyr::select(strain, treat, pca) %>%
    distinct %>%
    rowwise %>%
    mutate(pca = list(as.data.frame(pca$rotation) %>% rownames_to_column())) %>%
    unnest(pca) %>%
    rename(trait = rowname) %>%
    mutate(trait = factor(trait, levels = cilia_trait_levels, labels = cilia_trait_labels))

write.csv(cilia_pca, "data/ciliates-pca-loadings.csv", row.names = FALSE)

data_cilia <- dplyr::select(data_cilia, -pca, -pca_data)

# ==============================================================================
# Import and make the cyano data uniform
# ==============================================================================

# traits
cyano_trait_levels <- c("size", "chlorophyll", "phycocyanin", "phycoerythrin", "trait")
cyano_trait_labels <- c("Size", "Chlorophyll", "Phycocyanin", "Phycoerythrin", "Trait")

# actual data
data_cyano <- read_csv("data/cyanobacteria.csv") %>%
    # PCA
    group_by(species, strain, treat) %>%
    nest %>%
    rowwise %>%
    mutate(
        pca_data = list(dplyr::select(data, size:phycocyanin)),
        pca = list(prcomp(pca_data, center = TRUE, scale = TRUE)),
        data = list(cbind(data, pca1 = pca$x[, 1])),
        pca_varexp = (pca$sdev^2)[1] / sum(pca$sde^2)) %>%
    dplyr::select(species, strain, treat, data, pca, pca_varexp) %>%
    unnest(data) %>%
    rename(trait = pca1) %>%
    group_by(strain, treat, repl) %>%
    mutate(
        d_density = pcgr(density, day),
        across(all_of(cyano_trait_levels), ~ trait_change(.x, day),
            .names = "d_{.col}")) %>%
    dplyr::filter(!is.na(d_density)) %>%
    mutate(treat = factor(treat, levels = c("C", "T", "A", "AT")),
        strain = paste(species, "_", strain, sep = ""))

# ------------------------------------------------------------------------------
# save the PCA loadings
# ------------------------------------------------------------------------------

cyano_pca <- data_cyano %>%
    ungroup(repl) %>%
    dplyr::select(pca) %>%
    distinct %>%
    rowwise %>%
    mutate(pca = list(as.data.frame(pca$rotation) %>% rownames_to_column())) %>%
    unnest(pca) %>%
    rename(trait = rowname) %>%
    mutate(trait = factor(trait, levels = cyano_trait_levels, labels = cyano_trait_labels))

write.csv(cyano_pca, "data/cyanobacteria-pca-loadings.csv", row.names = FALSE)

data_cyano <- dplyr::select(data_cyano, -pca)

# ==============================================================================
# merge the data and save
# ==============================================================================

# merged data
data_merged <- bind_rows(
    mutate(data_cilia, system = "Ciliates"),
    mutate(data_cyano, system = "Cyanobacteria")) %>%
    relocate(system, species, strain, treat, repl, day, density,
        all_of(unique(c(cyano_trait_levels, cilia_trait_levels))),
        d_density,
        all_of(paste0("d_", c(cyano_trait_levels, cilia_trait_levels))))

write.csv(data_merged, "data/data_wide.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# pivot_data longer and save for density over time plots
# ------------------------------------------------------------------------------

long_data_cilia <- data_cilia %>%
    mutate(density = log10(density)) %>%
    pivot_longer(all_of(c("density", cilia_trait_levels))) 

long_data_cyano <- data_cyano %>%
    mutate(density = log10(density)) %>%
    pivot_longer(all_of(c("density", cyano_trait_levels)))

# long merged data
long_data_merged <- bind_rows(
    mutate(long_data_cilia, system = "Ciliates"),
    mutate(long_data_cyano, system = "Cyanobacteria")
    ) %>%
    relocate(system, species, strain, treat, repl, day,
        name, value,
        d_density,
        all_of(paste0("d_", c(cyano_trait_levels, cilia_trait_levels))))

write.csv(long_data_merged, "data/data_long.csv", row.names = FALSE)
