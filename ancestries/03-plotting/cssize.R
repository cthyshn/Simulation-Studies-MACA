### updated version that gives me the plots that I will use
library(data.table)
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)
library(tidyr)
library(cowplot)

cvsall <- read.table("../cvsall.txt") 
cvsall <- na.omit(cvsall)
plots <- list()
vec <- 1:24

for (i in vec) {
  
  cv <- cvsall$id[i]
  resultseven <- fread(paste0("../results/EUR:EAS:AFR/1:1:1/", cv, "/results.tsv"))
  resultsuneven <- fread(paste0("../results/EUR:EAS:AFR/75:15:10/", cv, "/results.tsv"))
  resultseven$setting <- "33:33:33"
  resultsuneven$setting <- "75:15:10"
  
  results <- rbind(resultseven, resultsuneven)
  
  
  cssize_long <- melt(results, id.vars = "setting", measure.vars = c("cs_size", "cm_size", "es_size", "em_size"), 
                      variable.name = "Method", value.name = "cssize")
  cssize_long$Method <- factor(cssize_long$Method, levels = c("cs_size", "cm_size", "es_size", "em_size"),
                               labels = c("coloc_SuSiEx", "coloc_MsCAVIAR", "eCAVIAR_SuSiEx", "eMsCAVIAR"))
  
  
  cssize_plot <- ggplot(na.omit(cssize_long), aes(x=setting, y=cssize, fill = Method)) + geom_boxplot(outlier.size = 0.75) + 
    theme_classic() + labs(y="Credible Set Size", fill = "Method", x = "EUR:AFR") +
    scale_fill_manual(values = c("coloc_SuSiEx" = "#f6a6b2",
                                 "eMsCAVIAR" = "#ffecb8",
                                 "eCAVIAR_SuSiEx" = "#b7ded2",
                                 "coloc_MsCAVIAR" = "#90d2d8")) +
    #   theme(legend.position = ifelse(i == ind[(length(ind))], "bottom", "none")) + 
    theme(legend.position = "none")   + 
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9, color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 9),
          # plot.title = element_text(size = 10, face = "bold"),
          plot.margin = margin(t=1,r=1,b=1,l=1, unit = "pt")) 
  
  
  plots[[i]] <- cssize_plot
  
}

legend <- get_legend(cssize_plot + theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)))
legend_row <- plot_grid(legend, nrow = 1)


for (j in seq_along(plots)) {
  plots[[j]] <- wrap_elements(plots[[j]]) + theme(legend.position = "none") + 
    plot_annotation(title = bquote(bold(.(LETTERS[j])* ".") ~ .(cvsall$id[vec][j])),
                    theme = theme(plot.title = element_text(size = 10)))
  plots[[j]] <- plot_grid(plots[[j]], ncol = 1, rel_heights = c(0.1, 1))
  
}

allplots <- wrap_plots(plots, ncol = 4)

final_plot <- plot_grid(allplots, legend_row, ncol = 1, rel_heights = c(5, 0.1))  
print(final_plot)

ggsave("cssize.pdf",
       plot = final_plot,
       width = 8.5, height = 16, units = "in")

ggsave("cssize.png",
       plot = final_plot,
       width = 8.5, height = 16, units = "in")

