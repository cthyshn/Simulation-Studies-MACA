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
  
  
  variantlevel_long <- melt(results, id.vars = "setting", measure.vars = c("cs_vl", "cm_vl", "es_vl", "em_vl"), 
                            variable.name = "Method", value.name = "variantlevel")
  variantlevel_long$Method <- factor(variantlevel_long$Method, levels = c("cs_vl", "cm_vl", "es_vl", "em_vl"),
                                     labels = c("coloc_SuSiEx", "coloc_MsCAVIAR", "eCAVIAR_SuSiEx", "eMsCAVIAR"))
  

  variantlevel_plot <- ggplot(na.omit(variantlevel_long), aes(x=setting, y=variantlevel, fill = Method)) + geom_boxplot(outlier.size = 0.75) +
    labs(y="cvCLPP", fill = "Method", x = "EUR:AFR") + theme_classic() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9, color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 9),
          # plot.title = element_text(size = 10, face = "bold"),
          plot.margin = margin(t=1,r=1,b=1,l=1, unit = "pt")) + ylim(c(0,1)) + 
    scale_fill_manual(values = c("coloc_SuSiEx" = "#f6a6b2",
                                 "eMsCAVIAR" = "#ffecb8",
                                 "eCAVIAR_SuSiEx" = "#b7ded2",
                                 "coloc_MsCAVIAR" = "#90d2d8")) + theme(legend.position = "none")
  

  plots[[i]] <- variantlevel_plot
  
}

legend <- get_legend(variantlevel_plot + theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)))
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



