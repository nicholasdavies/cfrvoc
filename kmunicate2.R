library(KMunicate)

# This is just the original KMunicate function, adding .legend_position, .xlim, .ylim, .margin, .title and .align arguments.
KMunicate2 = function (fit, time_scale, .risk_table = "KMunicate", .reverse = FALSE, 
    .theme = NULL, .color_scale = NULL, .fill_scale = NULL, .linetype_scale = NULL, 
    .annotate = NULL, .xlab = "Time", .ylab = ifelse(.reverse, 
        "Estimated (1 - survival)", "Estimated survival"), .alpha = 0.25, 
    .rel_heights = NULL, .ff = NULL, .risk_table_base_size = 11, 
    .size = NULL, .legend_position = c(1, 1), .xlim = NULL, .ylim = c(0, 1), .margin = 0, .title = NULL, .align = "hv") 
{
    arg_checks <- checkmate::makeAssertCollection()
    checkmate::assert_class(x = fit, classes = "survfit", add = arg_checks)
    checkmate::assert_numeric(x = time_scale, add = arg_checks)
    checkmate::assert_logical(x = .reverse, len = 1, add = arg_checks)
    checkmate::assert_string(x = .risk_table, null.ok = TRUE, 
        add = arg_checks)
    if (!is.null(.risk_table)) {
        .risk_table <- match.arg(.risk_table, choices = c("KMunicate", 
            "survfit", NULL))
        checkmate::assert_true(x = .risk_table %in% c("KMunicate", 
            "survfit"), add = arg_checks)
    }
    checkmate::assert_class(x = .theme, classes = c("theme", 
        "gg"), null.ok = TRUE, add = arg_checks)
    if (!is.null(.theme)) 
        checkmate::assert_true(x = ggplot2::is.theme(.theme), 
            add = arg_checks)
    if (!is.null(.color_scale)) 
        checkmate::assert_true(x = ggplot2::is.ggproto(.color_scale), 
            add = arg_checks)
    if (!is.null(.fill_scale)) 
        checkmate::assert_true(x = ggplot2::is.ggproto(.fill_scale), 
            add = arg_checks)
    if (!is.null(.linetype_scale)) 
        checkmate::assert_true(x = ggplot2::is.ggproto(.linetype_scale), 
            add = arg_checks)
    if (!is.null(.annotate)) 
        checkmate::assert_true(x = ggplot2::is.ggproto(.annotate), 
            add = arg_checks)
    checkmate::assert_string(x = .xlab, add = arg_checks)
    checkmate::assert_string(x = .ylab, add = arg_checks)
    checkmate::assert_number(x = .alpha, add = arg_checks)
    checkmate::assert_numeric(x = .rel_heights, null.ok = TRUE, 
        add = arg_checks)
    checkmate::assert_string(x = .ff, null.ok = TRUE, add = arg_checks)
    checkmate::assert_number(x = .risk_table_base_size, add = arg_checks)
    checkmate::assert_number(x = .size, null.ok = TRUE, add = arg_checks)
    checkmate::assert_true(x = (is.numeric(.legend_position) & 
        length(.legend_position) == 2) | (is.character(.legend_position) & 
        length(.legend_position) == 1), add = arg_checks)
    if (!arg_checks$isEmpty()) 
        checkmate::reportAssertions(arg_checks)
    data <- KMunicate:::.fortify(fit)
    if (.reverse) {
        data$surv <- 1 - data$surv
        data$lower <- 1 - data$lower
        data$upper <- 1 - data$upper
    }
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = time, y = surv, 
        ymin = lower, ymax = upper))
    if ("strata" %in% names(fit)) {
        plot <- plot + pammtools::geom_stepribbon(ggplot2::aes(fill = strata), 
            alpha = .alpha)
        if (!is.null(.size)) {
            plot <- plot + ggplot2::geom_step(ggplot2::aes(color = strata, 
                linetype = strata), size = .size)
        }
        else {
            plot <- plot + ggplot2::geom_step(ggplot2::aes(color = strata, 
                linetype = strata))
        }
    }
    else {
        if (!is.null(.size)) {
            plot <- plot + pammtools::geom_stepribbon(alpha = .alpha, 
                size = .size)
        }
        else {
            plot <- plot + pammtools::geom_stepribbon(alpha = .alpha)
        }
        plot <- plot + ggplot2::geom_step()
    }
    plot <- plot + ggplot2::scale_x_continuous(breaks = time_scale) + 
        ggplot2::coord_cartesian(ylim = .ylim, xlim = if (is.null(.xlim)) range(time_scale) else .xlim) + 
        ggplot2::labs(color = "", fill = "", linetype = "", x = .xlab, 
            y = .ylab, title = .title)
    if (!is.null(.theme)) {
        plot <- plot + .theme
    }
    else if (!is.null(.ff)) {
        plot <- plot + ggplot2::theme_gray(base_family = .ff)
    }
    plot <- plot + ggplot2::theme(legend.position = .legend_position, 
        legend.background = ggplot2::element_blank(), legend.key = ggplot2::element_blank(), 
        plot.margin = ggplot2::unit(c(.margin, .margin, .margin, .margin), "cm"))
    if (length(.legend_position) > 1) {
        plot <- plot + ggplot2::theme(legend.justification = .legend_position)
    }
    if (!is.null(.color_scale)) {
        plot <- plot + .color_scale
    }
    if (!is.null(.fill_scale)) {
        plot <- plot + .fill_scale
    }
    if (!is.null(.linetype_scale)) {
        plot <- plot + .linetype_scale
    }
    if (!is.null(.annotate)) {
        plot <- plot + .annotate
    }
    if (!is.null(.risk_table)) {
        table_data <- KMunicate:::.fortify_summary(fit = fit, time_scale = time_scale, 
            risk_table = .risk_table)
        table_data <- tidyr::pivot_longer(data = table_data, 
            cols = c("n.risk", "n.event", "n.censor"))
        table_data$name <- factor(table_data$name, levels = c("n.event", 
            "n.censor", "n.risk"), labels = c("Events", "Censored", 
            "At risk"))
        if (!("strata" %in% names(fit))) {
            table_data$strata <- "Overall"
        }
        tds <- split(table_data, f = table_data$strata)
        tds <- lapply(seq_along(tds), function(i) {
            p <- ggplot2::ggplot(tds[[i]], ggplot2::aes(x = time, 
                y = name, label = value))
            if (is.null(.ff)) {
                p <- p + ggplot2::geom_text(size = .risk_table_base_size/3)
            }
            else {
                p <- p + ggplot2::geom_text(mapping = ggplot2::aes(family = .ff), 
                  size = .risk_table_base_size/3)
            }
            p <- p + ggplot2::scale_x_continuous(breaks = time_scale) + 
                ggplot2::coord_cartesian(xlim = if (is.null(.xlim)) range(time_scale) else .xlim) + 
                ggplot2::theme_void(base_size = .risk_table_base_size) + 
                ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"), 
                    plot.title = element_text(size = .risk_table_base_size))
            if (is.null(.ff)) {
                p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
            }
            else {
                p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic", 
                  family = .ff), plot.title = ggplot2::element_text(family = .ff))
            }
            p <- p + ggplot2::labs(title = names(tds)[i])
        })
        if (is.null(.rel_heights)) {
            .rel_heights <- c(3, rep(1, length(tds)))
        }
        KMunicate_plot <- cowplot::plot_grid(plotlist = c(list(plot), 
            tds), align = .align, axis = "tlbr", ncol = 1, rel_heights = .rel_heights)
    }
    else {
        KMunicate_plot <- plot
    }
    return(KMunicate_plot)
}
