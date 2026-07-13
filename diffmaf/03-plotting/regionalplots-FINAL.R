### generates figure 4 in the main text
library(data.table)
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggbreak)
library(ggforce)
library(cowplot)

i=1

cvsall <- read.table("../cvsall2.txt")
plots <- list()

cv <- cvsall$id[i]
results5050 <- fread(paste0("../results/50:50/", cv, "/results.tsv"))
results8020 <- fread(paste0("../results/80:20/", cv, "/results.tsv"))
results5050$setting <- "50:50"
results8020$setting <- "80:20"

results <- rbind(results5050, results8020)

variantlevel_long <- melt(results, id.vars = "setting", measure.vars = c("cs_vl", "em_vl", "es_vl", "cm_vl"), variable.name = "Method", value.name = "variantlevel")


variantlevel_long$Method <- factor(
  variantlevel_long$Method,
  levels = c("cs_vl", "cm_vl", "es_vl", "em_vl"),
  labels = c("coloc_SuSiEx", "coloc_MsCAVIAR", "eCAVIAR_SuSiEx", "eMsCAVIAR")
)

locilevel_long <- melt(results, id.vars = "setting", measure.vars = c("cs_ll", "em_ll", "es_ll", "cm_ll"), variable.name = "Method", value.name = "locilevel")

locilevel_long$Method <- factor(
  locilevel_long$Method,
  levels = c("cs_ll", "cm_ll", "es_ll", "em_ll"),
  labels = c("coloc_SuSiEx", "coloc_MsCAVIAR", "eCAVIAR_SuSiEx", "eMsCAVIAR")
)


cssize_long <- melt(results, id.vars = "setting", measure.vars = c("cs_size", "em_size", "es_size", "cm_size"), variable.name = "Method", value.name = "cssize")


cssize_long$Method <- factor(
  cssize_long$Method,
  levels = c("cs_size", "cm_size", "es_size", "em_size"),
  labels = c("coloc_SuSiEx", "coloc_MsCAVIAR", "eCAVIAR_SuSiEx", "eMsCAVIAR")
)

sum_incs <- results %>% group_by(setting) %>%
  summarise(
    sum_cs_incs = sum(cs_incs, na.rm = TRUE)/100,
    sum_cm_incs = sum(cm_incs, na.rm = TRUE)/100,
    sum_es_incs = sum(es_incs, na.rm = TRUE)/100,
    sum_em_incs = sum(em_incs, na.rm = TRUE)/100
  )

sum_incs_long <- pivot_longer(
  sum_incs, 
  cols = c(sum_cs_incs, sum_cm_incs, sum_es_incs, sum_em_incs),
  names_to = "Method", 
  values_to = "sum_incs"
)

sum_incs_long$Method <- factor(sum_incs_long$Method, levels = c("sum_cs_incs", "sum_cm_incs", "sum_es_incs", "sum_em_incs"),
                               labels = c("coloc_SuSiEx", "coloc_MsCAVIAR", "eCAVIAR_SuSiEx", "eMsCAVIAR"))

sum_incs_long$se <- sqrt(sum_incs_long$sum_incs*(1-sum_incs_long$sum_incs)/100)
# sum_incs_long$se <- sqrt(100*sum_incs_long$sum_incs/100*(1-(sum_incs_long$sum_incs/100)))

sum_incs_long$barmin <- pmax(0, sum_incs_long$sum_incs - 1.96*sum_incs_long$se)
sum_incs_long$barmax <- pmin(sum_incs_long$sum_incs + 1.96*sum_incs_long$se, 1)



cssize_plot <- ggplot(na.omit(cssize_long), aes(x=setting, y=cssize, fill = Method)) + geom_boxplot(outlier.size = 0.75) + 
  theme_classic() + labs(y="95% Credible Set Size", fill = "Method", x = "EUR:AFR") +
  scale_fill_manual(values = c("coloc_SuSiEx" = "#f6a6b2",
                               "eMsCAVIAR" = "#ffecb8",
                               "eCAVIAR_SuSiEx" = "#b7ded2",
                               "coloc_MsCAVIAR" = "#90d2d8")) +
  #   theme(legend.position = ifelse(i == ind[(length(ind))], "bottom", "none")) + 
  theme(legend.position = "none")   + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 10, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold"),
        plot.margin = margin(t=1,r=1,b=1,l=1, unit = "pt")) 

incs_plot <- ggplot(sum_incs_long, aes(x = setting, y = sum_incs,
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
        legend.position = "none") + 
  
  coord_cartesian(clip = "off", ylim = c(0.50, 1.0), xlim = c(1, 2)) +
  scale_y_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0.2, 0.05)),
                     labels = scales::number_format(accuracy = 0.01),
                     breaks = c(0,0.50, 0.60, 0.70, 0.80, 0.90, 1)) +
  # scale_y_break(c(0.06, 0.84)) +
  geom_segment(aes(x=0.37, xend = 0.45, y= 0.455, yend = 0.459)) +
  geom_segment(aes(x=0.37, xend = 0.45, y= 0.470, yend = 0.474)) +
  annotate("text",
           x = 0.39,                     
           y = 0.43,                    
           label = "0.00",
           hjust = 1.25, vjust = 1,      
           size = 3) + 
  annotate("segment",
           x = 0.37, xend = 0.40 ,
           y = 0.420, yend = 0.420,
           color = "black")







variantlevel_plot <- ggplot(na.omit(variantlevel_long), aes(x=setting, y=variantlevel, fill = Method)) + geom_boxplot(outlier.size = 0.75) +
  labs(y="Conditional Variant Level PP", fill = "Method", x = "EUR:AFR") + theme_classic() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 10, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold"),
        plot.margin = margin(t=1,r=1,b=1,l=1, unit = "pt")) + ylim(c(0,1)) + 
  scale_fill_manual(values = c("coloc_SuSiEx" = "#f6a6b2",
                               "eMsCAVIAR" = "#ffecb8",
                               "eCAVIAR_SuSiEx" = "#b7ded2",
                               "coloc_MsCAVIAR" = "#90d2d8")) + theme(legend.position = "none") 

locilevel_plot <- ggplot(na.omit(locilevel_long), aes(x=setting, y=locilevel, fill = Method)) + geom_boxplot(outlier.size = 0.75) +
  labs(y="Loci Level PP", fill = "Method", x = "EUR:AFR") + theme_classic() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 10, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold"),
        plot.margin = margin(t=1,r=1,b=1,l=1, unit = "pt")) + ylim(c(0,1)) +
  scale_fill_manual(values = c("coloc_SuSiEx" = "#f6a6b2",
                               "eMsCAVIAR" = "#ffecb8",
                               "eCAVIAR_SuSiEx" = "#b7ded2",
                               "coloc_MsCAVIAR" = "#90d2d8")) + theme(legend.position = "none") 



plots <- cssize_plot + incs_plot + variantlevel_plot + locilevel_plot +
  plot_layout(ncol = 2, nrow = 2, heights = c(1.2, 1), widths = c(1, 1.2)) +
  plot_annotation(tag_levels = 'A')  





dotplot_legend <- get_legend(incs_plot + theme(legend.position = "bottom"))
boxplot_legend <- get_legend(cssize_plot + theme(legend.position = "bottom"))
legend_row <- plot_grid(boxplot_legend, dotplot_legend, nrow = 2)


final_plot <- plot_grid(plots, legend_row, ncol = 1, rel_heights = c(2, 0.3))
print(final_plot)

tiff("Figure4.tiff", units="in", width=6, height=6, res=300)
print(final_plot)
dev.off()
