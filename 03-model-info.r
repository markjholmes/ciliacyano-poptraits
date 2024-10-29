library("tidyverse") # general
library("qpcR") # for akaike weights

fitted_models <- readRDS("data/fitted-models.RData") %>%
    mutate(model_type = factor(model_type, levels = c("null", "single", "augmented")), 
        clade = str_split_i(strain, "_", 1))

bootstrapped_models <- readRDS("data/boot-test-results.RData") %>%
    dplyr::select(system, strain, treat, resp, size_or_trait, self, contains("PBLR"))

data_merged <- read.csv("data/data_wide.csv") %>%
    dplyr::select(system, strain, treat, pca_varexp) %>%
    distinct

common_results <- fitted_models %>%
    rowwise %>%
    mutate(
        obs = list(model$model[,1]),
        pred = list(as.numeric(predict(model))),
        aic = AICc(model),
        R2 = abs(summary(model)$r.squared)
        ) %>%
    left_join(data_merged) %>%
    relocate(system, .before = strain) %>%
    mutate(
        resp = factor(resp,
            levels = c("d_density", "d_size", "d_trait"),
            labels = c("Growth", "Size change", "Trait change")), 
        treat = factor(treat,
            levels = c("C", "T", "A", "AT"),
            labels = c("C", "T", "A", "AT"))) %>%
    arrange(size_or_trait, system, strain, treat, resp)

all_results <- common_results %>%
    group_by(system, strain, treat, resp, size_or_trait) %>%
    mutate(aic_rank = order(aic),
        R2_rank = order(R2, decreasing = TRUE), 
        a_weights = akaike.weights(aic)$weights) %>%
    dplyr::select(-model, -obs, -pred)

main_results <- common_results %>%
    dplyr::filter(size_or_trait == "trait", model_type != "null", self) %>%
    dplyr::select(-model) %>%
    ungroup %>%
    distinct %>%
    group_by(system, strain, treat, resp) %>%
    mutate(aic_rank = order(aic),
        R2_rank = order(R2, decreasing = TRUE), 
        a_weights = akaike.weights(aic)$weights) %>%
    dplyr::select(-size_or_trait)

# join models to data and make predictions ----
data_preds <- main_results %>%
    ungroup %>%
    unnest_longer(obs:pred)

# Now error plot
data_synth <- main_results %>%
    relocate(pca_varexp, .before = aic) %>%
    rowwise %>%
    mutate(error = sum(abs(pred - obs))) %>%
    dplyr::select(-obs, -pred) %>%
    group_by(system, strain, treat, resp) %>%
    pivot_wider(names_from = model_type, values_from = aic:error) %>% 
    mutate(delta_error = (error_single - error_augmented) / error_single)

# find the best model for each predictor on the basis of AIC
data_best_mods <- main_results %>%
    group_by(system, strain, treat, resp) %>%
    mutate(a_weights = akaike.weights(aic)$weights) %>%
    dplyr::filter(aic == min(aic))

# ------------------------------------------------------------------------------
# ALTERNATIVE RESULTS
# ------------------------------------------------------------------------------

size_results <- common_results %>%
    dplyr::filter(size_or_trait == "size", model_type != "null", self) %>%
    dplyr::select(-model) %>%
    ungroup %>%
    distinct %>%
    group_by(system, strain, treat, resp) %>%
    mutate(aic_rank = order(aic),
        R2_rank = order(R2, decreasing = TRUE), 
        a_weights = akaike.weights(aic)$weights) %>%
    dplyr::select(-size_or_trait)

# join models to data and make predictions ----
size_preds <- size_results %>%
    ungroup %>%
    unnest_longer(obs:pred)

# Now error plot
size_synth <- size_results %>%
    relocate(pca_varexp, .before = aic) %>%
    rowwise %>%
    mutate(error = sum(abs(pred - obs))) %>%
    dplyr::select(-obs, -pred) %>%
    group_by(system, strain, treat, resp) %>%
    pivot_wider(names_from = model_type, values_from = aic:error) %>% 
    mutate(delta_error = (error_single - error_augmented) / error_single)

# find the best model for each predictor on the basis of AIC
size_best_mods <- size_results %>%
    group_by(system, strain, treat, resp) %>%
    mutate(a_weights = akaike.weights(aic)$weights) %>%
    dplyr::filter(aic == min(aic))
