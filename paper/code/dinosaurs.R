source("functions.R")
set.seed(2024)

#### Maiasaura
data <- read.csv('data/Maiasaura.csv')
head(data)

maiasaura <- data.frame(
    age = data$Age..years. + data$Age.shift..years.,
    mass = data$Fraction.multiplied.to.ROM.44770.body.mass..kg.
)

maiasaura_plot <- ggplot(maiasaura, aes(age, mass)) +
    geom_point() +
    labs(x = "Age (years)", y = "Mass (kg)", title = "Maiasaura peeblesorum") +
    custom_theme()

initial_maiasaura_params <- list(
    L=2500,
    k=0.5,
    q=10
)

population_fit <- generate_logistic_curve_fit(maiasaura, initial_params = initial_maiasaura_params)
population_params <- coef(population_fit)

# generate the priors plot
num_samples <- 1e3
Ls <- rnorm(num_samples, population_params[['L']], sd = 500)
ks <- rnorm(num_samples, population_params[['k']], sd = 0.0125)
qs <- rnorm(num_samples, population_params[['q']], sd = 0.4)
sigmas <- rnorm(num_samples, 200, 50)

curve_data <- data.frame(age = numeric(), length = numeric())
curve_ages <- seq(min(maiasaura$age), max(maiasaura$age), by = 0.5)
for (i in 1:1e3) {
    index <- sample(1:num_samples, size = 1)
    L <- Ls[index]
    k <- ks[index]
    q <- qs[index]
    sigma <- sigmas[index]

    curve_mass_mean <- generalized_logistic(curve_ages, L = L, k = k, q = q)
    curve_mass <- rnorm(length(curve_mass_mean), mean = curve_mass_mean, sd = sigma)
    curve_data <- rbind(
        curve_data,
        data.frame(age = curve_ages, mass = curve_mass)
    )
}

range_df <- data.frame(age = numeric(), lower = numeric(), upper = numeric())
for (age in unique(curve_data$age)) {
    relevant <- curve_data[abs(curve_data$age - age) < 0.1, ]
    masses <- relevant$mass
    lower <- quantile(masses, c(0.025))
    upper <- quantile(masses, c(0.975))
    range_df <- rbind(
        range_df,
        data.frame(
            age = age,
            lower = lower,
            upper = upper
        )
    )
}
p <- ggplot() +
    geom_point(data = maiasaura, mapping = aes(age, mass)) +
    geom_ribbon(data = range_df, mapping = aes(x = age, ymin = lower, ymax = upper), alpha = 0.2) +
    labs(x = "Age (years)", y = "Mass (kg)", title = "Prior 95% Quantiles") +
    custom_theme()

ggsave("../images/maiasauraPrior.png", plot = p, width = 15, height = 8)

# now do the regression
maiasaura$residual <- resid(population_fit)
maiasaura$sex <- ifelse(maiasaura$residual < 0, "S", "L")
large_sex_data <- subset(maiasaura, sex == "L")
small_sex_data <- subset(maiasaura, sex == "S")

model_data <- list(
    N_large = nrow(large_sex_data),
    age_large = large_sex_data$age,
    mass_large = large_sex_data$mass,
    N_small = nrow(small_sex_data),
    age_small = small_sex_data$age,
    mass_small = small_sex_data$mass,
    prior_L_mean = population_params[['L']],
    prior_L_sd = 500,
    prior_k_mean = population_params[['k']],
    prior_k_sd = 0.0125,
    prior_q_mean = population_params[['q']],
    prior_q_sd = 0.4,
    prior_sigma_mean = 200,
    prior_sigma_sd = 50
)

maiasaura.model <- stan("logistic_dimorphism_model.stan", data = model_data)
maiasaura_samples <- extract(maiasaura.model)
lapply(maiasaura_samples, head)

p <- ggplot() +
    custom_theme() +
    geom_point(maiasaura, mapping = aes(age, mass, colour = sex), size = 3) +
    scale_colour_grey() +
    labs(x = "Age (years)", y = "Mass (kg)", colour = "Sex")

# draw in samples from the male and female curves
for (i in 1:100) {
    index <- sample(1:length(maiasaura_samples[[ 'L_small' ]]), size = 1)
    L_small <- maiasaura_samples[['L_small']][index]
    k_small <- maiasaura_samples[['k_small']][index]
    q_small <- maiasaura_samples[['q_small']][index]
    curve_mass_small <- generalized_logistic(curve_ages, L = L_small, k = k_small, q = q_small)

    L_large <- maiasaura_samples[['L_large']][index]
    k_large <- maiasaura_samples[['k_large']][index]
    q_large <- maiasaura_samples[['q_large']][index]
    curve_mass_large <- generalized_logistic(curve_ages, L = L_large, k = k_large, q = q_large)

    p <- p + 
        geom_line(data = data.frame(age = curve_ages, mass = curve_mass_small, sex = "S"), mapping = aes(age, mass, colour = sex), alpha = 0.1) +
        geom_line(data = data.frame(age = curve_ages, mass = curve_mass_large, sex = "L"), mapping = aes(age, mass, colour = sex), alpha = 0.1)

}
p <- p + labs(title = "Maiasaura peeblesorum Data and Posterior Growth Curves")
ggsave('../images/maiasauraPosterior.png', plot = p, width = 15, height = 8)

# actual amount of dimorphism
p <- ggplot(data.frame(x = maiasaura_samples[['difference']]), aes(x)) +
    geom_density(aes(y = after_stat(density))) +
    labs(x = expression(L[l] - L[s]), y = "Density", title = "Posterior Distribution for Maiasaura Sexual Dimorphism") +
    custom_theme()
ggsave('../images/maiasauraDimorphism.png', plot = p, width = 15, height = 8)






########### Psittacosaurus
data <- read.csv('data/Psittacosaurus.csv')
head(data)

psittacosaurus <- data.frame(
    age = data$Age.estimate.shifted..years.,
    mass = data$Estimated.body.mass..kg.
)

initial_psittacosaurus_params <- list(
    L=37.38,
    k=0.55,
    q=10
)

fit.population <- generate_logistic_curve_fit(psittacosaurus, initial_params = initial_psittacosaurus_params)
psittacosaurus_population_params <- coef(fit.population)

# Find priors
num_samples <- 1e3
prior_L_sd <- 3
Ls <- rnorm(num_samples, psittacosaurus_population_params[['L']], sd = prior_L_sd)
prior_k_sd <- 0.00125
ks <- rnorm(num_samples, psittacosaurus_population_params[['k']], sd = prior_k_sd)
prior_q_sd <- 0.0004
qs <- rnorm(num_samples, psittacosaurus_population_params[['q']], sd = prior_q_sd)
prior_sigma_mean <- 3
prior_sigma_sd <- 0.1
sigmas <- rnorm(num_samples, prior_sigma_mean, prior_sigma_sd)

curve_data <- data.frame(age = numeric(), length = numeric())
curve_ages <- seq(min(psittacosaurus$age), 2 * max(psittacosaurus$age), by = 0.5)
for (i in 1:1e3) {
    index <- sample(1:num_samples, size = 1)
    L <- Ls[index]
    k <- ks[index]
    q <- qs[index]
    sigma <- sigmas[index]

    curve_mass_mean <- generalized_logistic(curve_ages, L = L, k = k, q = q)
    curve_mass <- rnorm(length(curve_mass_mean), mean = curve_mass_mean, sd = sigma)
    curve_data <- rbind(
        curve_data,
        data.frame(age = curve_ages, mass = curve_mass)
    )
}

range_df <- data.frame(age = numeric(), lower = numeric(), upper = numeric())
for (age in unique(curve_data$age)) {
    relevant <- curve_data[abs(curve_data$age - age) < 0.1, ]
    masses <- relevant$mass
    lower <- quantile(masses, c(0.025))
    upper <- quantile(masses, c(0.975))
    range_df <- rbind(
        range_df,
        data.frame(
            age = age,
            lower = lower,
            upper = upper
        )
    )
}
p <- ggplot() +
    geom_point(data = psittacosaurus, mapping = aes(age, mass)) +
    geom_ribbon(data = range_df, mapping = aes(x = age, ymin = lower, ymax = upper), alpha = 0.2) +
    labs(x = "Age (years)", y = "Mass (kg)", title = "Prior 95% Quantiles") +
    custom_theme()

ggsave("../images/psittacosaurusPrior.png", plot = p, width = 15, height = 8)

# now fit the psittacosaurus data
psittacosaurus$residual <- resid(fit.population)
psittacosaurus$sex <- ifelse(psittacosaurus$residual < 0, "S", "L")
large_sex_data <- subset(psittacosaurus, sex == "L")
small_sex_data <- subset(psittacosaurus, sex == "S")

model_data <- list(
    N_large = nrow(large_sex_data),
    age_large = large_sex_data$age,
    mass_large = large_sex_data$mass,
    N_small = nrow(small_sex_data),
    age_small = small_sex_data$age,
    mass_small = small_sex_data$mass,
    prior_L_mean = psittacosaurus_population_params[['L']],
    prior_L_sd = prior_L_sd,
    prior_k_mean = psittacosaurus_population_params[['k']],
    prior_k_sd = prior_k_sd,
    prior_q_mean = psittacosaurus_population_params[['q']],
    prior_q_sd = prior_q_sd,
    prior_sigma_mean = prior_sigma_mean,
    prior_sigma_sd = prior_sigma_sd
)

psittacosaurus.model <- stan("logistic_dimorphism_model.stan", data = model_data)
psittacosaurus.model

psittacosaurus_samples <- extract(psittacosaurus.model)

p <- ggplot() +
    custom_theme() +
    scale_colour_grey() +
    geom_point(psittacosaurus, mapping = aes(age, mass, colour = sex), size = 3) +
    labs(x = "Age (years)", y = "Mass (kg)", colour = "Sex")

# draw in samples from the large and small curves
for (i in 1:100) {
    index <- sample(1:length(psittacosaurus_samples[[ 'L_small' ]]), size = 1)
    L_small <- psittacosaurus_samples[['L_small']][index]
    k_small <- psittacosaurus_samples[['k_small']][index]
    q_small <- psittacosaurus_samples[['q_small']][index]
    curve_mass_small <- generalized_logistic(curve_ages, L = L_small, k = k_small, q = q_small)

    L_large <- psittacosaurus_samples[['L_large']][index]
    k_large <- psittacosaurus_samples[['k_large']][index]
    q_large <- psittacosaurus_samples[['q_large']][index]
    curve_mass_large <- generalized_logistic(curve_ages, L = L_large, k = k_large, q = q_large)

    p <- p + 
        geom_line(data = data.frame(age = curve_ages, mass = curve_mass_small, sex = "S"), mapping = aes(age, mass, colour = sex), alpha = 0.1) +
        geom_line(data = data.frame(age = curve_ages, mass = curve_mass_large, sex = "L"), mapping = aes(age, mass, colour = sex), alpha = 0.1)

}
p <- p + labs(title = "Psittacosaurus lujiatunensis Data and Posterior Growth Curves")
ggsave("../images/psittacosaurusPosterior.png", plot = p, width = 15, height = 8)

p <- ggplot(data.frame(x = psittacosaurus_samples[['difference']]), aes(x)) +
    geom_density(aes(y = after_stat(density))) +
    labs(x = expression(L[l] - L[s]), y = "Density", title = "Posterior Distribution for Psittacosaurus Sexual Dimorphism") +
    custom_theme()
ggsave('../images/psittacosaurusDimorphism.png', width = 15, height = 8)







########### Tyrannosaurus
data <- read.csv('data/Tyrannosaurus.csv')
head(data)

tyrannosaurus <- data.frame(
    age = data$Age..years.,
    mass = data$Mass..kg.
)

initial_tyrannosaurus_params <- list(
    L=5000,
    k=0.5,
    q=10
)

fit.population <- generate_logistic_curve_fit(tyrannosaurus, initial_params = initial_tyrannosaurus_params)
tyrannosaurus_population_params <- coef(fit.population)

# Find priors
num_samples <- 1e3
prior_L_sd <- 1000
Ls <- rnorm(num_samples, tyrannosaurus_population_params[['L']], sd = prior_L_sd)
prior_k_sd <- 0.0125
ks <- rnorm(num_samples, tyrannosaurus_population_params[['k']], sd = prior_k_sd)
prior_q_sd <- 0.004
qs <- rnorm(num_samples, tyrannosaurus_population_params[['q']], sd = prior_q_sd)
prior_sigma_mean <- 100
prior_sigma_sd <- 1
sigmas <- rnorm(num_samples, prior_sigma_mean, prior_sigma_sd)

curve_data <- data.frame(age = numeric(), length = numeric())
curve_ages <- seq(min(tyrannosaurus$age), 2 * max(tyrannosaurus$age), by = 0.5)
for (i in 1:1e3) {
    index <- sample(1:num_samples, size = 1)
    L <- Ls[index]
    k <- ks[index]
    q <- qs[index]
    sigma <- sigmas[index]

    curve_mass_mean <- generalized_logistic(curve_ages, L = L, k = k, q = q)
    curve_mass <- rnorm(length(curve_mass_mean), mean = curve_mass_mean, sd = sigma)
    curve_data <- rbind(
        curve_data,
        data.frame(age = curve_ages, mass = curve_mass)
    )
}

range_df <- data.frame(age = numeric(), lower = numeric(), upper = numeric())
for (age in unique(curve_data$age)) {
    relevant <- curve_data[abs(curve_data$age - age) < 0.1, ]
    masses <- relevant$mass
    lower <- quantile(masses, c(0.025))
    upper <- quantile(masses, c(0.975))
    range_df <- rbind(
        range_df,
        data.frame(
            age = age,
            lower = lower,
            upper = upper
        )
    )
}
p <- ggplot() +
    geom_point(data = tyrannosaurus, mapping = aes(age, mass)) +
    geom_ribbon(data = range_df, mapping = aes(x = age, ymin = lower, ymax = upper), alpha = 0.2) +
    labs(x = "Age (years)", y = "Mass (kg)", title = "Prior 95% Quantiles") +
    custom_theme()
ggsave('../images/tyrannosaurPrior.png', plot = p, width = 15, height = 8)

# now fit the tyrannosaurus data
tyrannosaurus$residual <- resid(fit.population)
tyrannosaurus$sex <- ifelse(tyrannosaurus$residual < 0, "S", "L")
large_sex_data <- subset(tyrannosaurus, sex == "L")
small_sex_data <- subset(tyrannosaurus, sex == "S")

model_data <- list(
    N_large = nrow(large_sex_data),
    age_large = large_sex_data$age,
    mass_large = large_sex_data$mass,
    N_small = nrow(small_sex_data),
    age_small = small_sex_data$age,
    mass_small = small_sex_data$mass,
    prior_L_mean = tyrannosaurus_population_params[['L']],
    prior_L_sd = prior_L_sd,
    prior_k_mean = tyrannosaurus_population_params[['k']],
    prior_k_sd = prior_k_sd,
    prior_q_mean = tyrannosaurus_population_params[['q']],
    prior_q_sd = prior_q_sd,
    prior_sigma_mean = prior_sigma_mean,
    prior_sigma_sd = prior_sigma_sd
)

tyrannosaurus.model <- stan("logistic_dimorphism_model.stan", data = model_data)
tyrannosaurus.model

tyrannosaurus_samples <- extract(tyrannosaurus.model)
p <- ggplot() +
    custom_theme() +
    scale_colour_grey() +
    geom_point(tyrannosaurus, mapping = aes(age, mass, colour = sex), size = 3) +
    labs(x = "Age (years)", y = "Mass (kg)", colour = "Sex")

# draw in samples from the large and small curves
for (i in 1:100) {
    index <- sample(1:length(tyrannosaurus_samples[[ 'L_small' ]]), size = 1)
    L_small <- tyrannosaurus_samples[['L_small']][index]
    k_small <- tyrannosaurus_samples[['k_small']][index]
    q_small <- tyrannosaurus_samples[['q_small']][index]
    curve_mass_small <- generalized_logistic(curve_ages, L = L_small, k = k_small, q = q_small)

    L_large <- tyrannosaurus_samples[['L_large']][index]
    k_large <- tyrannosaurus_samples[['k_large']][index]
    q_large <- tyrannosaurus_samples[['q_large']][index]
    curve_mass_large <- generalized_logistic(curve_ages, L = L_large, k = k_large, q = q_large)

    p <- p + 
        geom_line(data = data.frame(age = curve_ages, mass = curve_mass_small, sex = "S"), mapping = aes(age, mass, colour = sex), alpha = 0.1) +
        geom_line(data = data.frame(age = curve_ages, mass = curve_mass_large, sex = "L"), mapping = aes(age, mass, colour = sex), alpha = 0.1)

}
p <- p + labs(title = "Tyrannosaurus rex Data and Posterior Growth Curves")
ggsave('../images/tyrannosaurPosterior.png', plot = p, width = 15, height = 8)

p <- ggplot(data.frame(x = tyrannosaurus_samples[['difference']]), aes(x)) +
    geom_density(aes(y = after_stat(density))) +
    labs(x = expression(L[l] - L[s]), y = "Density", title = "Posterior Distribution for Tyrannosaurus Sexual Dimorphism") +
    custom_theme()
ggsave('../images/tyrannosaurDimorphism.png', plot = p, width = 15, height = 8)







###### Comparison
dimorphism_df <- rbind(
    data.frame(
        x = maiasaura_samples[['relative_difference']] - 1,
        species = "Maiasaura peeblesorum"
    ),
    data.frame(
        x = psittacosaurus_samples[['relative_difference']] - 1,
        species = "Psittacosaurus lujiatunensis"
    ),
    data.frame(
        x = tyrannosaurus_samples[['relative_difference']] - 1,
        species = "Tyrannosaurus rex"
    )
)

# values from recreation of saitta et al
maiasaura_dimorphism <- 1.474308346902
psittacosaurus_dimorphism <- 1.20719588632773
tyrannosaur_dimorphism <- 1.12623031006368

saitta_dimorphism_df <- data.frame(
    species = c("Maiasaura peeblesorum", "Psittacosaurus lujiatunensis", "Tyrannosaurus rex"),
    dimorphism = c(maiasaura_dimorphism, psittacosaurus_dimorphism, tyrannosaur_dimorphism) - 1
)

p <- ggplot(dimorphism_df, aes(x = x, colour = species)) +
    custom_theme() +
    scale_colour_grey() +
    scale_fill_grey() +
    geom_density(aes(y = after_stat(density), fill = species), alpha = 0.2, linewidth = 2) +
    geom_vline(data = saitta_dimorphism_df, mapping = aes(xintercept = dimorphism, colour = species, linetype = "Saitta et al (2020)"), linewidth = 2) +
    scale_linetype_manual(values = c("Saitta et al (2020)" = 'dashed')) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Percent Dimorphism", y = "Density", title = "A Comparison of Different Distributions of Dimorphism", colour = "Species", fill = "Species", linetype = NULL)
ggsave("../images/combinedDimorphism.png", plot = p, width = 15, height = 8)
