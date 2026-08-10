# eQTL violin plot with beta and p-value annotation
#   Simulated_test_data_geneA_eQTL.csv
# with columns:
#   subject_id, genotype_effect_allele_dosage, age_years, sex, batch,
#   timepoint, hours_post_infection, infection_context, geneA_expression_log2

library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(interactions)
library(rlang)
library(scales)

# Read data
df <- read.csv("~/Simulated_test_data_geneA_eQTL.csv")

# Make genotype a factor for the violin plot
df <- df %>%
  mutate(
    genotype = factor(genotype_effect_allele_dosage, levels = c(0, 1, 2)),
    exp = geneA_expression_log2
  )

# Fit mixed-effects eQTL model
fit <- lmer(exp ~ genotype_effect_allele_dosage + age_years + sex + batch + (1 | subject_id), data = df)
print(summary(fit))


# Remove covariates but NOT genotype
fit0 <- lmer(exp ~ age_years + sex + batch + (1 | subject_id), data = df)
# Residualised expression
df$corrected_exp <- resid(fit0) + mean(df$exp, na.rm = TRUE)


# Extract beta and p-value for genotype
fit_sum <- summary(fit)
beta <- coef(fit_sum)["genotype_effect_allele_dosage","Estimate"]
P <- coef(fit_sum)["genotype_effect_allele_dosage", "Pr(>|t|)"]

# Add group sizes for labels
data4 <- df %>%
  group_by(genotype) %>%
  mutate(n = n()) %>%
  ungroup()

data4$name <- paste0(data4$genotype, "\n(n=", data4$n, ")")

# Plot
ggplot(data4, aes(x = name, y = corrected_exp)) +
  ggtitle(bquote(atop("Gene A eQTL",
        italic(P) == .(signif(P, 3)) * ";" ~ beta == .(signif(beta, 3)))))+
  ylab("Gene expression") +
  xlab("Genotype") +
  geom_violin(trim = FALSE, fill = "darkorange", color = "grey", alpha = 0.5) +
  geom_jitter(shape = 16, position = position_jitter(width = 0.1), alpha = 0.5, size = 1) +
  geom_boxplot(outlier.size = 0, fill = "white", width = 0.4, alpha = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

# Save plot
#ggsave("geneA_violin_plot_with_beta_pvalue.png", width = 8, height = 5, dpi = 300)


#------------------------------------------------------------------------------#
## does the eQTL effect interact with infection time points
#------------------------------------------------------------------------------#

df$timepoint <- factor(df$timepoint,
                       levels = c("baseline", "6h_post_infection", "24h_post_infection"))
#
# Main-effects model
fit_main <- lmer(
  corrected_exp ~ genotype_effect_allele_dosage + timepoint + (1 | subject_id),
  data = df,
  REML = FALSE)

# Interaction model
fit_int <- lmer(
  corrected_exp ~ genotype_effect_allele_dosage * timepoint + (1 | subject_id),
  data = df,
  REML = FALSE)

# Overall interaction test
lrt <- anova(fit_main, fit_int)
P.int <- lrt$`Pr(>Chisq)`[2]

# Standard interaction plot
int.plot <- interact_plot(
  fit_int,
  pred = genotype_effect_allele_dosage,
  modx = !!sym("timepoint"),
  plot.points = TRUE,
  jitter = c(0.1, 0),
  colors = c("grey", "darkorange", "red"),
  centered = "none",
  interval = TRUE,
  partial.residuals = FALSE)

#
int.plot +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = ",")) +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  ggtitle(
    bquote(atop(
      "Gene A eQTL by timepoint",
      italic(P)[interaction] ~ "=" ~ .(format.pval(P.int, digits = 2, eps = 1e-3)))))



