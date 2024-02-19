source("functions.R")

set.seed(2024)

prior_h <- c(160, 10)
actual_h <- c(176, 7)

sample <- rnorm(10, actual_h[1], actual_h[2])

heights <- seq(140, 200, length.out = 1000)
prior_likelihood <- dnorm(heights, prior_h[1], prior_h[2])
likelihood_likelihood <- sapply(heights, function(h) sum(dnorm(sample, h, actual_h[2])))
posterior <- prior_likelihood * likelihood_likelihood

prior_likelihood <- prior_likelihood / sum(prior_likelihood)
likelihood_likelihood <- likelihood_likelihood / sum(likelihood_likelihood)
posterior <- posterior / sum(posterior)

df <- data.frame(
    height = rep(heights, 3),
    density = c(prior_likelihood, likelihood_likelihood, posterior),
    type = rep(c("Prior", "Likelihood", "Posterior"), each = length(heights))
)

p <- ggplot(df, aes(height, density, group = type, colour = type)) +
    geom_line(linewidth = 2) +
    labs(x = "Height (cm)", y = "Density", colour = "Type") +
    scale_y_continuous(labels = NULL) +
    # scale_colour_grey() +
    theme()

paper_image <- paper_colours(p) +
    paper_theme(legend.position = c(0.8, 0.8), legend.background = element_rect(fill = 'white'))
poster_image <- poster_colours(p) + 
    poster_theme(legend.position = c(0.8, 0.8), legend.background = element_rect(fill = 'white'))

ggsave_with_defaults(generate_filename("bayesianExample.png"), plot = paper_image)
ggsave_with_defaults(generate_filename("bayesianExample-poster.png"), plot = poster_image)