source("functions.R")

# for reproducibility
set.seed(2024)

sample_sizes <- c(5, 10, 20, 50, 100)
sample_df <- data.frame(sample_size = integer(), diff = numeric())
for (sample_size in sample_sizes) {
    male_ages <- runif(sample_size, 1, ALLIGATOR_MAX_LIFESPAN)
    male_length <- generate_male_sample(male_ages)
    female_ages <- runif(sample_size, 1, ALLIGATOR_MAX_LIFESPAN)
    female_length <- generate_female_sample(female_ages)
    combined_data <- data.frame(
        age = c(male_ages, female_ages),
        length = c(male_length, female_length),
        sex = rep(c("M", "F"), each = sample_size)
    )
    population_fit <- generate_vb_curve_fit(combined_data)
    population_params <- coef(population_fit)
    prior_L <- c(population_params[['L']], 0.5)
    prior_K <- c(population_params[['K']], 0.025)
    prior_A <- c(population_params[['A']], 0.25)
    prior_sigma <- c(sigma, 0.05)
    combined_data$predicted_sex <- generate_predicted_sex(combined_data, fit = population_fit)
    male_data <- subset(combined_data, predicted_sex == "M")
    female_data <- subset(combined_data, predicted_sex == "F")

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
        data.frame(sample_size = sample_size * 2, diff = samples[['diff']])
    )
}

alpha <- 0.05
summary_df <- data.frame(
    sample_size = integer(),
    mean = numeric(),
    lower = numeric(),
    upper = numeric()
)

for (sample_size in unique(sample_df$sample_size)) {
    relevant <- sample_df[abs(sample_df$sample_size - sample_size) < 0.001, ]
    sample_mean <- mean(relevant$diff)
    bounds <- quantile(relevant$diff, c(alpha / 2, 1 - alpha / 2))
    summary_df <- rbind(
        summary_df,
        data.frame(sample_size = sample_size, mean = sample_mean, lower = bounds[1], upper = bounds[2])
    )
}

zero_sample <- rnorm(8000, 0, sqrt(2) * prior_L[2])
zero_mean <- mean(zero_sample)
bounds <- quantile(zero_sample, c(alpha / 2, 1 - alpha / 2))

summary_df <- rbind(
    summary_df,
    data.frame(sample_size = 0, mean = zero_mean, lower = bounds[1], upper = bounds[2])
)
actual_diff <- actual_male_params[['L']] - actual_female_params[['L']]
summary_df$sample_size <- factor(summary_df$sample_size)
p <- ggplot(summary_df, aes(sample_size, mean)) +
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    geom_hline(yintercept = actual_diff, linewidth = 2, linetype = 'dashed') +
    annotate('label', label = "Actual Level", x = 4.5, y = 0.7, size = unit(8, 'pt'), fill = 'white', label.size = NA) +
    labs(x = "Total Sample Size", y = "Dimorphism Level (m)", title = "Dimorphism Level (Mean and 95% Credible Interval) and Sample Size") +
    custom_theme()
ggsave("../images/alligatorSampleSize.png", plot = p, width = 15, height = 8)