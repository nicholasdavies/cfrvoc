# Fig. 1
# Requires plots from: survival.R, kaplan_meier.R, death_rates.R.

pl_fig1 = plot_grid(
    plot_grid(
        plot_grid(
            plot_grid(plot_samples, plot_deaths, ncol = 1, labels = letters, label_size = 10, align = "v", axis = "bottom"),
            plAA60_both, 
            plBB60_both, 
            plCC60_both, 
            nrow = 2, labels = c("", letters[3:5]), label_size = 10)
    ),
    pl_deathrates_sgtf,
    nrow = 1,
    rel_widths = c(1, 1)
)

ggsave("./output/fig1.pdf", pl_fig1, width = 50, height = 25, units = "cm", useDingbats = FALSE)
ggsave("./output/fig1.png", pl_fig1, width = 50, height = 25, units = "cm")

# Fig. 2
# Requires plots from survival.R

pl_fig2 = plot_grid(
    plot_grid(hp_sgtf, hp_voc, ncol = 1, labels = letters[1:2], label_size = 10),
    pl_effects, labels = c("", letters[3]), label_size = 10, nrow = 1
)

ggsave("./output/fig2.pdf", pl_fig2, width = 30, height = 25, units = "cm", useDingbats = FALSE)
ggsave("./output/fig2.png", pl_fig2, width = 30, height = 25, units = "cm")
