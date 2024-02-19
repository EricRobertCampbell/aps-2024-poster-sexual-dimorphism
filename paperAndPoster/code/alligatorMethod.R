source("functions.R")

# generate the sample
num_samples <- 1e3
combined_data <- generate_combined_alligator_sample(generate_male_sample, generate_female_sample, num_samples)

population_fit <- generate_vb_curve_fit(combined_data)
population_params <- coef(population_fit)

curve_ages <- seq(1, ALLIGATOR_MAX_LIFESPAN, length.out = 100)
curve_lengths <- generalized_von_bertalanffy(curve_ages, L = population_params[['L']], A = population_params[['A']], K = population_params[['K']])
curve_data <- data.frame(age = curve_ages, length = curve_lengths)

p <- ggplot(combined_data, aes(age, length)) +
    geom_point(aes(shape = sex, colour = sex)) +
    geom_line(data = curve_data, mapping = aes(age, length), linewidth = 2) +
    annotate('label', label = "Predicted Male", x = 25, y = 3.7, size = unit(8, 'pt'), fill = 'white', label.size = NA) +
    annotate('label', label = "Predicted Female", x = 25, y = 1.2, size = unit(8, 'pt'), fill = 'white', label.size = NA) +
    annotate('label', label = "Population Curve", x = 35, y = 2.8, size = unit(8, 'pt'), fill = 'white', label.size = NA) +
    labs(x = "Age (years)", y = "Length (m)", title = "Alligator Population and Curve", colour = "Actual Sex", shape = "Actual Sex")

paper_image <- paper_colours(p) + 
    paper_theme(legend.position = c(0.85, 0.25), legend.background = element_rect(fill = 'white'))
poster_image <- poster_colours(p) + 
    poster_theme(legend.position = c(0.85, 0.25), legend.background = element_rect(fill = 'white'))
ggsave_with_defaults(generate_filename("alligatorMethod.png"), plot = paper_image)
ggsave_with_defaults(generate_filename("alligatorMethod-poster.png"), plot = poster_image)
