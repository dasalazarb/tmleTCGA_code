library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)

load("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/os_lgg_data_var50/tmle_comparison_var50.rda")
var50 <- tmle_comparison
load("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/os_lgg_data_var100/tmle_comparison_var100.rda")
var100 <- tmle_comparison
load("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/os_lgg_data_var500/tmle_comparison_var500.rda")
var500 <- tmle_comparison
load("G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/os_lgg_data_var1000/tmle_comparison_var1000.rda")
var1000 <- tmle_comparison; rm(tmle_comparison)


## ------------------------------------------------------------ ##
### ... Graficos y m?s graficos de tmle y del superlearner ... ###
## ------------------------------------------------------------ ##

### ATE plot
giveMeplot <- function(tmle_comparison) {
  psi.ate <- ggplot(tmle_comparison, aes(x=reorder(learner_strategy, -ATE.psi), y=ATE.psi)) + 
    geom_line() + 
    geom_point() + 
    geom_hline(yintercept=0, linetype="dashed", color = "black") + 
    geom_text_repel(label = round(tmle_comparison$ATE.psi, digits = 4), size=3) + 
    geom_errorbar(aes(ymin=ATE.CI.LO,ymax=ATE.CI.UP), width=.2,
                  position=position_dodge(0.05)) + 
    theme(axis.text.x = element_text(angle = 0, hjust=.75,vjust=0.95), 
          text = element_text(size = 13)) + xlab("Method")
  return(psi.ate)
}

plot50 <- giveMeplot(var50)
plot100 <- giveMeplot(var100)
plot500 <- giveMeplot(var500)
plot1000 <- giveMeplot(var1000)

all.ate.plots <- ggarrange(plot50, plot100, plot500, plot1000, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

###first set of plots

jpeg(filename = "G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/all.ate.plots.jpg", width = 15, height = 8, units = "in", res=300)
all.ate.plots
dev.off()

# tiff(filename = paste0(ruta, "os_lgg_data_var", nvar, "/fig_ATE_RR_OR_var", nvar, ".tiff"), width = 13, height = 8, units = "in", res=300)
# ggarrange(psi, or.rr, ncol = 2)
# dev.off()

pdf(file = "G:/Mi unidad/Tutorial_2018-10/03_Resultados/DataTCGA/tmleTCGA/all.ate.plots.pdf", width = 15, height = 8)
all.ate.plots
dev.off()
