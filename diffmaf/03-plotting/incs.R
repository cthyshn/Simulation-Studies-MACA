library(data.table)
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)
library(tidyr)
library(cowplot)

cvsall <- read.table("../cvsall2.txt") 
cvsall <- na.omit(cvsall)
plots <- list()
vec <- 1:12

for (i in vec) {
  
  cv <- cvsall$id[i]
  results5050 <- fread(paste0("../results/50:50/", cv, "/results.tsv"))
  results8020 <- fread(paste0("../results/80:20/", cv, "/results.tsv"))
  results5050$setting <- "50:50"
  results8020$setting <- "80:20"
  
  results <- rbind(results5050, results8020)
  
  
  prop_incs <- results %>% group_by(setting) %>%
    summarise(
      prop_cs_incs = sum(cs_incs, na.rm = TRUE)/100,
      prop_em_incs = sum(em_incs, na.rm = TRUE)/100,
      prop_es_incs = sum(es_incs, na.rm = TRUE)/100,
      prop_cm_incs = sum(cm_incs, na.rm = TRUE)/100
    )
  prop_incs_long <- pivot_longer(
    prop_incs, 
    cols = c(prop_cs_incs, prop_em_incs, prop_es_incs, prop_cm_incs),
    names_to = "Method", 
    values_to = "prop_incs"
  )
  
  prop_incs_long$Method <- factor(
    prop_incs_long$Method,
    levels = c("prop_cs_incs", "prop_cm_incs", "prop_es_incs", "prop_em_incs"),
    labels = c("coloc_SuSiEx", "coloc_MsCAVIAR", "eCAVIAR_SuSiEx", "eMsCAVIAR")
  )
  
  prop_incs_long$se <- sqrt(prop_incs_long$prop_incs*(1-prop_incs_long$prop_incs)/100)
  
  prop_incs_long$barmin <- pmax(0, prop_incs_long$prop_incs - 1.96*prop_incs_long$se)
  prop_incs_long$barmax <- pmin(prop_incs_long$prop_incs + 1.96*prop_incs_long$se, 1)
  
  
  incs_plot <- ggplot(prop_incs_long, aes(x = (factor(setting)), y = prop_incs,
                                          fill = Method, shape = Method)) +
    geom_point(position = position_dodge(width = 0.7),
               size = 3, stroke = 0.5) + 
    geom_hline(yintercept=0.95, linetype = "longdash", color = "red", linewidth = 0.25) +
    geom_errorbar(aes(ymin=barmin, ymax=barmax), width=.2,
                  position=position_dodge(.7)) +
    theme_classic() +
    labs(y = "Coverage", fill = "Method", shape = "Method", x = "EUR:AFR") +
    
    scale_fill_manual(values = c("coloc_SuSiEx" = "#f6a6b2",
                                 "eMsCAVIAR" = "#ffecb8",
                                 "eCAVIAR_SuSiEx" = "#b7ded2",
                                 "coloc_MsCAVIAR" = "#90d2d8")) +
    
    scale_shape_manual(values = c("coloc_SuSiEx" = 21,   # Circle
                                  "eMsCAVIAR" = 22,      # Square
                                  "eCAVIAR_SuSiEx" = 23, # Diamond
                                  "coloc_MsCAVIAR" = 24  # Triangle
    )) +
    
    theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
          legend.position = "none") +   scale_y_continuous(expand = expansion(mult = c(0.2, 0.05)),
                                                           labels = scales::number_format(accuracy = 0.01))
  
  
  plots[[i]] <- incs_plot
  
}

legend <- get_legend(incs_plot + theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)))
legend_row <- plot_grid(legend, nrow = 1)


for (j in seq_along(plots)) {
  plots[[j]] <- wrap_elements(plots[[j]]) + theme(legend.position = "none") + 
    plot_annotation(title = bquote(bold(.(LETTERS[j])* ".") ~ .(cvsall$id[vec][j])),
                    theme = theme(plot.title = element_text(size = 10)))
  plots[[j]] <- plot_grid(plots[[j]], ncol = 1, rel_heights = c(0.1, 1))
  
}

allplots <- wrap_plots(plots, ncol = 4)

final_plot <- plot_grid(allplots, legend_row, ncol = 1, rel_heights = c(5, 0.5))  
print(final_plot)



