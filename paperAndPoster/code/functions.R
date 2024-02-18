
library(ggplot2)
library(FSA) # for the vbStart function to find good parameters for the von Bertalanffy function
library(minpack.lm) # a more robust version of nls for non-linear curve fitting to avoid problems with e.g. very small sample sizes
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("constants.R")

custom_theme <- function() {
  theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = alpha("black", 0.5), linetype = "dotted", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),  # Adjust the font size of axis labels
      axis.title = element_text(size = 16),  # Adjust the font size of axis titles
      plot.title = element_text(size = 20, hjust = 0.5),  # Adjust the font size and center the plot title
      plot.subtitle = element_text(size = 16, hjust = 0.5)  # Adjust the font size and center the plot subtitle
    )
}

paper_colours <- function() {
    scale_colour_grey() +
    scale_fill_grey()
}

poster_colours <- function() {

}

generalized_von_bertalanffy <- function(age, L, A, K) {
    L * (1 - A * exp(-K * age))
}

generalized_logistic <- function(age, L, q, k) {
    L / (1 + exp(q + k * age))
}

generate_initial_vb_params <- function(df, debug = FALSE) {
    initial_params_raw <- tryCatch({vbStarts(df$length ~ df$age)},
        error = function(e) {
            if (debug) {
                cat("Error getting initial von Bertalanffy parameters for data", conditionMessage(e))
            }
            NULL
        }
    )

    if (is.null(initial_params_raw)) return(NULL)

    # now convert to my parameterization
    # I use y = Linf * (1 - A e ^ (-Kt) )
    # they use y = Linf * (1 - e^(-K (t - t0)))
    # -> A = e^K t0
    initial_params <- list(
        L = initial_params_raw$Linf,
        K = initial_params_raw$K,
        A = exp(initial_params_raw$K * initial_params_raw$t0)
    )
    initial_params
}

generate_vb_curve_fit <- function(data) {
    start <- generate_initial_vb_params(data)
    fit <- nlsLM(length ~ generalized_von_bertalanffy(age, L, A, K), data = data, start = start)
    fit
}

generate_logistic_curve_fit <- function(df, initial_params, debug = FALSE) {
    fit <- tryCatch({ nlsLM(mass ~ generalized_logistic(age, L, q, k), data = df, start = initial_params) }, error = function(e) {
            cat("Error in generate_curve_fit", conditionMessage(e))
            print(df)
            NULL
        }
    )
    fit
}

generate_generate_alligator_samples <- function(params) {
    function(ages) {
        means <- generalized_von_bertalanffy(ages, L = params[['L']], A = params[["A"]], K = params[['K']])
        rnorm(length(ages), mean = means, sd = params[['sigma']])
    }
}
generate_female_sample <- generate_generate_alligator_samples(actual_female_params)
generate_male_sample <- generate_generate_alligator_samples(actual_male_params)

generate_combined_alligator_sample <- function(male_sampler, female_sampler, num_samples) {
    male_ages <- runif(num_samples, 1, ALLIGATOR_MAX_LIFESPAN)
    male_lengths <- male_sampler(male_ages)
    female_ages <- runif(num_samples, 1, ALLIGATOR_MAX_LIFESPAN)
    female_lengths <- female_sampler(female_ages)
    data.frame(
        age = c(male_ages, female_ages),
        length = c(male_lengths, female_lengths),
        sex = rep(c("M", "F"), each = num_samples)
    )
}

generate_predicted_sex <- function(combined_data, fit = NULL) {
    if (is.null(fit)) {
        fit <- generate_curve_fit(combined_data)
    }
    residuals <- resid(fit)
    ifelse(residuals < 0, "F", "M")
}

generate_fake_male_params <- function(E) {
    list(
        L = actual_female_params[['L']] + E * (actual_male_params[['L']] - actual_female_params[['L']]),
        A = actual_female_params[['A']] + E * (actual_male_params[['A']] - actual_female_params[['A']]),
        K = actual_female_params[['K']] + E * (actual_male_params[['K']] - actual_female_params[['K']]),
        sigma = actual_female_params[['sigma']] + E * (actual_male_params[['sigma']] - actual_female_params[['sigma']])
    )
}

ggsave_with_defaults <- function(filename, plot, width = 15, height = 8, dpi = 600, ...) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
  sprintf("Done saving %s", filename)
}

generate_filename <- function(name) sprintf("../images/%s", name)

paper_colours <- function() scale_colour_grey()
poster_colours <- function() scale_fill_viridis_d()