library("tidyverse") # general

# ==============================================================================
# fit models
# ==============================================================================

# all functions for use inside mutate

# likelihood ratio of y : x
compute_LR <- function(x, y) {
    2 * (x - y)
} 

extract_p <- function(x, y) {
    ecdf(x)(y)
}

# 5000 pretty slow to replicate main text, default is 100
boot_test <- function(mod, nsim = 100L) {

    mod_formula <- format(formula(mod))

    response <- str_split_i(format(formula(mod)), " ~ ", 1)
    predictors <- setdiff(names(coef(mod)), "(Intercept)")

    sub_formula <- paste0(response, " ~ ", c("1", predictors))
    sub_models <- lapply(sub_formula, lm, mod$model)

    mod_LL <- logLik(mod)
    sub_LLs <- lapply(sub_models, logLik)
    obs_LR <- lapply(sub_LLs, compute_LR, x = mod_LL)

    sim_sub_LLs <- setNames(vector("list", length(sub_formula)), sub_formula)
    sim_full_LLs <- setNames(vector("list", length(sub_formula)), sub_formula)
    sim_PBLRS <- setNames(vector("list", length(sub_formula)), sub_formula)

    for (i in 1:length(sub_formula)) {
        form <- sub_formula[[i]]
        sub_model <- sub_models[[i]]
        sims <- as.list(simulate(sub_model, nsim))

        for (sim in sims) {
            sim_df <- mod$model 
            sim_df[,response] <- sim
            # this is the slow step
            sim_mod <- lm(sub_model, sim_df)
            full_mod <- lm(mod_formula, sim_df)
            # put the LLS on their lists
            sim_sub_LLs[[i]] <- append(sim_sub_LLs[[i]], logLik(sim_mod))
            sim_full_LLs[[i]] <- append(sim_full_LLs[[i]], logLik(full_mod))
        }
    }

    PBLR <- mapply(compute_LR, sim_full_LLs, sim_sub_LLs, SIMPLIFY = FALSE)
    PBLR_test <- mapply(extract_p, PBLR, obs_LR, SIMPLIFY = FALSE)

    print(paste0(mod_formula, " done"))

    return(list(PBLR_test))

}

source("01-fit-models.R")

models_for_boot <- modelling(
        dat = data_merged,
        var_to_nest_by = c("system", "strain", "treat"),
        formulas = formula_df %>%
            dplyr::filter(model_type == "augmented") %>%
            pull(formula)) %>%
    left_join(formula_df, by = "formula", relationship = "many-to-many")
    
bootstrapped_models <- models_for_boot %>%
    rowwise %>%
    mutate(PBLR_test = boot_test(model, nsim = 5000L)) %>% 
    unnest_longer(PBLR_test) %>%
    dplyr::select(-formula) %>%
    left_join(models_for_boot)

saveRDS(bootstrapped_models, "data/boot-test-results.RData")
