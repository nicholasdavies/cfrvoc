library(writexl)

# Fig. 1
# Requires plots from: survival.R, kaplan_meier.R, death_rates.R.

pl_fig1 = plot_grid(
    plot_grid(
        plot_grid(
            plot_samples, 
            plot_deaths, 
            ncol = 1, labels = letters, label_size = 11, align = "v", axis = "bottom"
        ),
        plAA60_both, 
        ncol = 1, rel_heights = c(1, 1), labels = c("", "c"), label_size = 11
    ),
    pl_deathrates_sgtf,
    nrow = 1,
    rel_widths = c(0.5, 1)
)

ggsave("./output/fig1.pdf", pl_fig1, width = 30, height = 18, units = "cm", useDingbats = FALSE)
ggsave("./output/fig1.png", pl_fig1, width = 30, height = 18, units = "cm")

# Write source data
write_xlsx(
    list(
        `Fig 1a` = plot_samples$data,
        `Fig 1b` = plot_deaths$data,
        `Fig 1c` = plAA60_inset$data,
        `Fig 1d` = rbind(pla$layers[[1]]$data, pla$layers[[3]]$data, fill = TRUE),
        `Fig 1e` = rbind(plb$layers[[1]]$data, plb$layers[[3]]$data, fill = TRUE),
        `Fig 1f` = rbind(plc$layers[[1]]$data, plc$layers[[3]]$data, fill = TRUE),
        `Fig 1g` = rbind(pld$layers[[1]]$data, pld$layers[[3]]$data, fill = TRUE),
        `Fig 1h` = rbind(ple$layers[[1]]$data, ple$layers[[3]]$data, fill = TRUE),
        `Fig 1i` = rbind(plf$layers[[1]]$data, plf$layers[[3]]$data, fill = TRUE)
    ), "./manuscript/sd1.xlsx")



# Fig. 2
# Requires plots from survival.R

library(Cairo)

pl_fig2 = plot_grid(
    plot_grid(
        hp_sgtf_cc + annotate("text", x = 28, y = 0.25, hjust = 1, vjust = 0.5, 
            label = "list(SGTF,complete~cases)", parse = TRUE, size = 9/ggplot2:::.pt) + labs(x = NULL) + ylim(c(0, 3.1)), 
        hp_sgtf_ipw + annotate("text", x = 28, y = 0.25, hjust = 1, vjust = 0.5, 
            label = "list(SGTF,IPW)", parse = TRUE, size = 9/ggplot2:::.pt) + labs(x = NULL) + ylim(c(0, 3.1)), 
        hp_voc_cc + annotate("text", x = 28, y = 0.25, hjust = 1, vjust = 0.5, 
            label = "list(p[VOC],complete~cases)", parse = TRUE, size = 9/ggplot2:::.pt) + labs(x = NULL) + ylim(c(0, 3.1)), 
        hp_voc_ipw + annotate("text", x = 28, y = 0.25, hjust = 1, vjust = 0.5, 
            label = "list(p[VOC],IPW)", parse = TRUE, size = 9/ggplot2:::.pt) + labs(x = NULL) + ylim(c(0, 3.1)), 
        nrow = 1, labels = letters[1:4], label_size = 11, align = "h", rel_widths = c(1.09, 1, 1, 1), label_x = c(0, -0.05, -0.05, -0.05)),
    ggdraw() + annotate("text", x = 0.52, y = 1, label = "Days since positive test", size = 11/ggplot2:::.pt),
    pl_sensitivity, labels = c("", "", letters[5]), label_size = 11, nrow = 3, rel_heights = c(3.8, 0.2, 25)
)

ggsave("./output/fig2.pdf", pl_fig2, width = 20, height = 29, units = "cm", device = cairo_pdf)
ggsave("./output/fig2.png", pl_fig2, width = 20, height = 29, units = "cm")

# Write source data
data2e = fread("./output/table_sensitivity_effects.csv")
write_xlsx(
    list(
        `Fig 2a` = hp_sgtf_cc$data,
        `Fig 2b` = hp_sgtf_ipw$data,
        `Fig 2c` = hp_voc_cc$data,
        `Fig 2d` = hp_voc_ipw$data,
        `Fig 2e` = data2e
    ), "./manuscript/sd2.xlsx")


# Viral load

theme_set(theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

plO = ggplot(dataS[specimen_date >= "2021-01-01"]) + 
    geom_density(aes(x = ctORF1ab, colour = factor(ifelse(sgtf == 1, "SGTF", "Non-SGTF"), levels = c("SGTF", "Non-SGTF")))) + 
    facet_wrap(~NHSER_name) +
    labs(x = "Ct ORF1ab", y = "Density", colour = NULL) +
    theme(legend.position = c(0.4, 0.2))

plN = ggplot(dataS[specimen_date >= "2021-01-01"]) + 
    geom_density(aes(x = ctN, colour = factor(ifelse(sgtf == 1, "SGTF", "Non-SGTF"), levels = c("SGTF", "Non-SGTF")))) + 
    facet_wrap(~NHSER_name) +
    labs(x = "Ct N", y = "Density", colour = NULL) +
    theme(legend.position = c(0.4, 0.2))

plot_grid(plO, plN, nrow = 1, labels = letters, label_size = 10)
ggsave("./output/ct.pdf", width = 25, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/ct.png", width = 25, height = 12, units = "cm")

# source data
counts = dataS[specimen_date >= "2020-01-01"]
counts[, ct_ORF1ab := cut(ctORF1ab, 0:30)]
counts[, ct_N := cut(ctN, 0:30)]

library(writexl)

write_xlsx(list(
    ORF1ab = counts[, .N, keyby = .(sgtf, NHSER_name, ct_ORF1ab)],
    N = counts[, .N, keyby = .(sgtf, NHSER_name, ct_N)]
), "./manuscript/sdE_ct.xlsx")

dataS[specimen_date >= "2021-01-01", mean(ctORF1ab), by = .(sgtf, NHSER_name)]
dataS[specimen_date >= "2021-01-01", mean(ctN), by = .(sgtf, NHSER_name)]
