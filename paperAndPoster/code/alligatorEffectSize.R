source("functions.R")

set.seed(2024)
alpha <- 0.05
num_samples <- 1e3 / 2
effect_sizes <- seq(0, 2, length.out = 10)
sample_df <- data.frame(E = numeric(), sample_diff = numeric(), actual_diff = numeric(), sex_prediction_accuracy = numeric())
for (E in effect_sizes) {
    # generate the fake males
    fake_male_params <- generate_fake_male_params(E)
    generate_fake_male_sample <- generate_generate_alligator_samples(fake_male_params)
    
    # samples
    female_age <- runif(num_samples, 1, ALLIGATOR_MAX_LIFESPAN)
    female_length <- generate_female_sample(female_age)
    male_age <- runif(num_samples, 1, ALLIGATOR_MAX_LIFESPAN)
    male_length <- generate_fake_male_sample(male_age)
    combined_data <- data.frame(
        age = c(female_age, male_age),
        length = c(female_length, male_length),
        sex = rep(c("F", "M"), each = num_samples)
    )
    combined_data <- subset(combined_data, length > 0)

    population_fit <- generate_vb_curve_fit(combined_data)
    population_params <- coef(population_fit)

    prior_L <- c(population_params[['L']], 0.5)
    prior_K <- c(population_params[['K']], 0.025)
    prior_A <- c(population_params[['A']], 0.25)
    prior_sigma <- c(sigma, 0.05)
    
    combined_data$predicted_sex <- generate_predicted_sex(combined_data, fit = population_fit)
    male_data <- subset(combined_data, predicted_sex == "M")
    female_data <- subset(combined_data, predicted_sex == "F")
    sex_prediction_accuracy <- sum(combined_data$sex == combined_data$predicted_sex) / nrow(combined_data)

    data <- list(
            N_male = nrow(male_data),
            age_male = male_data$age,
            length_male = male_data$length,
            N_female = nrow(female_data),
            age_female = female_data$age,
            length_female = female_data$length,
            prior_L = prior_L,
            prior_A = prior_A,
            prior_K = prior_K,
            prior_sigma = prior_sigma
        )
    model <- stan("vb_dimorphism_model.stan", data = data, iter = 8000)
    samples <- extract(model)
    sample_df <- rbind(
        sample_df,
        data.frame(E = E, sample_diff = samples[['diff']], actual_diff = fake_male_params[['L']] - actual_female_params[['L']],
        sex_prediction_accuracy = sex_prediction_accuracy
        )
    )
}

sample_df$error <- sample_df$sample_diff - sample_df$actual_diff

p <- ggplot(sample_df, aes(error, group = E, colour = E)) +
    geom_density(aes(y = after_stat(density))) +
    # scale_color_distiller(
    #     type = "seq",
    #     direction = -1,
    #     palette = "Greys"
    # ) +
    labs(x = "Error (Predicted - Actual)", y = "Density", title = "Effect Size on Dimorphism Error", colour = "Effect Size") +
    scale_y_continuous(labels = NULL)
paper_image <- p +
    scale_color_distiller(
        type = "seq",
        direction = -1,
        palette = "Greys"
    ) +
    paper_theme(legend.position = c(0.2, 0.8), legend.background = element_rect(fill = 'white'))
poster_image <- p +
    scale_color_distiller(
        type = "seq",
        direction = -1,
        # palette = "Reds"
        palette = "Greys"
    ) +
    poster_theme(legend.position = c(0.2, 0.8), legend.background = element_rect(fill = 'white'))
ggsave_with_defaults(generate_filename("alligatorEffectSize.png"), plot = paper_image)
ggsave_with_defaults(generate_filename("alligatorEffectSize-poster.png"), plot = poster_image)

# now the error by sex prediction error
error_by_error_df <- data.frame(
    E = numeric(),
    mean_error = numeric(),
    lower = numeric(),
    upper = numeric(),
    sex_prediction_accuracy = numeric()
)
for (E in unique(sample_df$E)) {
    relevant <- sample_df[abs(sample_df$E - E) < 0.000001, ]
    mean_error <- mean(relevant$error)
    bounds <- quantile(relevant$error, c(alpha / 2, 1 - alpha / 2))
    sex_prediction_accuracy <- sum(sample_df$sex == sample_df$predicted_sex)
    error_by_error_df <- rbind(
        error_by_error_df,
        data.frame(
            E = E,
            mean_error = mean_error,
            lower = bounds[1],
            upper = bounds[2],
            sex_prediction_accuracy = unique(relevant$sex_prediction_accuracy)
        )
    )
}

p <- ggplot(error_by_error_df, aes(sex_prediction_accuracy, mean_error, colour = E)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    geom_smooth(method = 'lm', se = FALSE, colour = 'black') +
    # scale_color_distiller(
    #     type = "seq",
    #     direction = -1,
    #     palette = "Greys"
    # ) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Sex Prediction Accuracy", y = "Error (Predicted - Actual)", colour = "Effect Size", title = "Relationship Between Dimorphism Error and Sex Prediction Accuracy")

paper_image <- p +
    scale_color_distiller(
        type = "seq",
        direction = -1,
        palette = "Greys"
    ) +
    paper_theme(legend.position = c(0.9, 0.8), legend.background = element_rect(fill = 'white'))
poster_image <- p +
    scale_color_distiller(
        type = "seq",
        direction = -1,
        # palette = "Reds"
        palette = "Greys"
    ) +
    poster_theme(legend.position = c(0.8, 0.8), legend.background = element_rect(fill = 'white'))
ggsave_with_defaults(generate_filename('alligatorSexPredictionAccuracyDimorphismError.png'), plot = paper_image)
ggsave_with_defaults(generate_filename('alligatorSexPredictionAccuracyDimorphismError-poster.png'), plot = poster_image)