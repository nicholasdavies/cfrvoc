# Fig. 1
# Requires plots from: survival.R, kaplan_meier.R, death_rates.R.

# pl_fig1_old = plot_grid(
#     plot_grid(
#         plot_grid(
#             plot_grid(plot_samples, plot_deaths, ncol = 1, labels = letters, label_size = 10, align = "v", axis = "bottom"),
#             plAA60_both, 
#             plBB60_both, 
#             plCC60_both, 
#             nrow = 2, labels = c("", letters[3:5]), label_size = 10)
#     ),
#     pl_deathrates_sgtf,
#     nrow = 1,
#     rel_widths = c(1, 0.9)
# )

pl_fig1 = plot_grid(
    plot_grid(
        plot_grid(
            plot_samples, 
            plot_deaths, 
            ncol = 1, labels = letters, label_size = 10, align = "v", axis = "bottom"
        ),
        plAA60_both, 
        ncol = 1, rel_heights = c(0.8, 1), labels = c("", "c"), label_size = 10
    ),
    pl_deathrates_sgtf,
    nrow = 1,
    rel_widths = c(0.5, 1)
)

ggsave("./output/fig1.png", pl_fig1, width = 40, height = 20, units = "cm")
ggsave("./output/fig1.pdf", pl_fig1, width = 40, height = 20, units = "cm", useDingbats = FALSE)



# Fig. 2
# Requires plots from survival.R

pl_fig2 = plot_grid(
    plot_grid(hp_sgtf_cc, hp_voc_cc, hp_sgtf_ipw, hp_voc_ipw, ncol = 1, labels = letters[1:4], label_size = 10, align = "hv"),
    pl_effects, labels = c("", letters[5]), label_size = 10, nrow = 1, rel_widths = c(10, 20)
)

#ggsave("./output/fig2.pdf", pl_fig2, width = 30, height = 25, units = "cm", useDingbats = FALSE)
ggsave("./output/fig2.png", pl_fig2, width = 30, height = 21, units = "cm")



# Viral load

theme_set(theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

plO = ggplot(dataS[specimen_date >= "2021-01-01"]) + 
    geom_density(aes(x = ctORF1ab, colour = factor(ifelse(sgtf == 1, "SGTF", "Other"), levels = c("SGTF", "Other")))) + 
    facet_wrap(~NHSER_name) +
    labs(x = "Ct ORF1ab", y = "Density", colour = NULL) +
    theme(legend.position = c(0.4, 0.2))

plN = ggplot(dataS[specimen_date >= "2021-01-01"]) + 
    geom_density(aes(x = ctN, colour = factor(ifelse(sgtf == 1, "SGTF", "Other"), levels = c("SGTF", "Other")))) + 
    facet_wrap(~NHSER_name) +
    labs(x = "Ct N", y = "Density", colour = NULL) +
    theme(legend.position = c(0.4, 0.2))

plot_grid(plO, plN, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/ct.pdf", width = 25, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/ct.png", width = 25, height = 12, units = "cm")

dataS[specimen_date >= "2021-01-01", mean(ctORF1ab), by = .(sgtf, NHSER_name)]
dataS[specimen_date >= "2021-01-01", mean(ctN), by = .(sgtf, NHSER_name)]
