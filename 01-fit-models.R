library("tidyverse") # general

# ==============================================================================
# load in the data
# ==============================================================================

# function to do the modelling
modelling <- function(dat = dat, var_to_nest_by = "strain", formulas) {
    
    var_to_nest_by <- c(var_to_nest_by, "formula")

    test <- expand_grid(ungroup(dat), formula = formulas) %>%
        nest_by(across({{var_to_nest_by}})) %>%
        rename(dat = data) %>%
        ungroup %>% rowwise %>%
        mutate(model = list(lm(as.formula(formula), dat = dat))) %>%
        ungroup

    return(test)
}

data_merged <- as_tibble(read.csv("data/data_wide.csv")) %>%
    mutate(treat = factor(treat, levels = c("C", "T", "A", "AT")))

# ==============================================================================
# model fitting of size, not aggregate trait
# ==============================================================================

responses <- c("d_density", "d_size", "d_trait")
intercept <- c("1")
predictors <- c("density", "size", "trait")

null_formulae <- expand_grid("resp" = responses, "pred" = intercept) %>%
    rowwise %>%
    mutate(formula = paste0(resp, " ~ ", pred),
        model_type = "null")

single_formulae <- expand_grid("resp" = responses, "pred" = predictors) %>% 
    rowwise %>%
    mutate(formula = paste0(resp, " ~ ", pred),
        model_type = "single")

aug_formulae <- expand_grid(
    "resp" = responses,
    "pred1" = predictors,
    "pred2" = predictors) %>% 
    dplyr::filter(pred1 != pred2) %>%
    # remove reverses
    rowwise %>%
    mutate(preds = paste0(sort(c(pred1, pred2)), collapse = " ")) %>%
    distinct(resp, preds, .keep_all = TRUE) %>%
    mutate(formula = paste0(resp, " ~ ", pred1, " + ", pred2),
        model_type = "augmented") %>%
    mutate(pred = paste0(pred1, " + ", pred2)) %>%
    dplyr::select(resp, pred, formula, model_type)

all_formulae <- bind_rows(list(null_formulae, single_formulae, aug_formulae)) %>% ungroup

# safest way to duplicate the d_dens ~ dens formula in both is to split into not-trait and not-size and then bind them back together

size_formulae <- all_formulae %>%
    dplyr::filter(!grepl("trait", formula)) %>%
    mutate(size_or_trait = "size")

trait_formulae <- all_formulae %>%
    dplyr::filter(!grepl("size", formula)) %>%
    mutate(size_or_trait = "trait")

formula_df <- bind_rows(size_formulae, trait_formulae) %>%
    rowwise %>%
    mutate(self = ifelse(grepl(str_split_i(resp, "_", 2), pred), TRUE, FALSE)) %>%
    ungroup

# ==============================================================================
# fit models
# ==============================================================================

fitted_size <- modelling(
        dat = data_merged,
        var_to_nest_by = c("system", "strain", "treat"),
        formulas = size_formulae$formula)

fitted_trait <- modelling(
        dat = data_merged,
        var_to_nest_by = c("system", "strain", "treat"),
        formulas = trait_formulae$formula)

fitted_models <- bind_rows(fitted_size, fitted_trait) %>%
    left_join(formula_df, by = "formula", relationship = "many-to-many") %>%
    ungroup %>%
    dplyr::select(-dat, -formula)

# save
saveRDS(fitted_models, "data/fitted-models.RData")

