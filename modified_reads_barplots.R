library(tidyverse)
library(reshape2)
library(scales)
library(ggpubr)

data1 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV678_3days_Cas9_/CRISPResso_quantification_of_editing_frequency.txt")

data2 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV678_3days_Cas9-/CRISPResso_quantification_of_editing_frequency.txt")

d <- rbind(data1,data2) %>% as.data.frame()

d <- cbind(c("Cas9+ 3days", "Cas9- 3days"),d[,4]) %>% cbind(d[,6:8])
colnames(d)[1:2] <- c("Sample", "Total_reads")


d2 <- melt(d, id.vars = "Sample")
d3 <- filter(d2, Sample == "Cas9+ 3days") %>% mutate(pct = (as.integer(value)/1070876)*100)
d4 <- filter(d2, Sample == "Cas9- 3days") %>% mutate(pct = (as.integer(value)/1103296)*100) %>% rbind(d3)
colnames(d4)[2] <- "Type"

p <- ggplot(d4, aes(x = Sample, y = as.integer(value))) +
  geom_bar( aes(fill = as.factor(Type)), stat="identity",position=position_dodge(), colour = "black") +
  geom_text(mapping = aes(label = format(round(pct, 1), nsmall = 1)), data = d4, vjust=-0.5, color="black",position = position_dodge(0.87), size=5.5)+
  scale_fill_hue(c = 40)


p + scale_y_continuous(name = "Reads", limits = c(0,1200000), breaks = c(200000, 400000, 600000, 800000,1000000,1200000) ) +
  theme(legend.title=element_blank(),
        axis.ticks = element_line(size = 1),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.box.background = element_rect(),
        legend.text = element_text(face="bold", color="black", size = 12),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size=16,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=14),
        axis.text.y = element_text(face="bold", color="black", size=14),
        panel.grid.major.y = element_line(color = "darkgray", linetype = "dashed"),
        axis.line = element_line(colour = "black", size = 0.8),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),)
#p
########################################################################## LV641 3days #####################################################################################

data1 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_3days_Cas9_/CRISPResso_quantification_of_editing_frequency.txt")

data2 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_3days_Cas9-/CRISPResso_quantification_of_editing_frequency.txt")

d <- rbind(data1,data2) %>% as.data.frame()

d <- cbind(c("Cas9+ 3days", "Cas9- 3days"),d[,4]) %>% cbind(d[,6:8])
colnames(d)[1:2] <- c("Sample", "Total_reads")


d2 <- melt(d, id.vars = "Sample")
d3 <- filter(d2, Sample == "Cas9+ 3days") %>% mutate(pct = (as.integer(value)/879992)*100)
d4 <- filter(d2, Sample == "Cas9- 3days") %>% mutate(pct = (as.integer(value)/1035369)*100) %>% rbind(d3)
colnames(d4)[2] <- "Type"



p2 <- ggplot(d4, aes(x = Sample, y = as.integer(value))) +
  geom_bar( aes(fill = as.factor(Type)), stat="identity",position=position_dodge(), colour = "black") +
  geom_text(mapping = aes(label = format(round(pct, 1), nsmall = 1)), data = d4, vjust=-0.5, color="black",position = position_dodge(0.87), size=5.5)+
  scale_fill_hue(c = 40)


p2 + scale_y_continuous(name = "Reads", limits = c(0,1200000), breaks = c(200000, 400000, 600000, 800000,1000000,1200000) ) +
  theme(legend.title=element_blank(),
        axis.ticks = element_line(size = 1),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.box.background = element_rect(),
        legend.text = element_text(face="bold", color="black", size = 12),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size=16,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=14),
        axis.text.y = element_text(face="bold", color="black", size=14),
        panel.grid.major.y = element_line(color = "darkgray", linetype = "dashed"),
        axis.line = element_line(colour = "black", size = 0.8),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),)



#################################################################### LV641 15days #################################################################

data1 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_15days_Cas9_/CRISPResso_quantification_of_editing_frequency.txt")

data2 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_15days_Cas9-/CRISPResso_quantification_of_editing_frequency.txt")

d <- rbind(data1,data2) %>% as.data.frame()

d <- cbind(c("Cas9+ 15days", "Cas9- 15days"),d[,4]) %>% cbind(d[,6:8])
colnames(d)[1:2] <- c("Sample", "Total_reads")


d2 <- melt(d, id.vars = "Sample")
d3 <- filter(d2, Sample == "Cas9+ 15days") %>% mutate(pct = (as.integer(value)/2217740)*100)
d4 <- filter(d2, Sample == "Cas9- 15days") %>% mutate(pct = (as.integer(value)/2089539)*100) %>% rbind(d3)
colnames(d4)[2] <- "Type"


p3 <- ggplot(d4, aes(x = Sample, y = as.integer(value))) +
  geom_bar( aes(fill = as.factor(Type)), stat="identity",position=position_dodge(), colour = "black") +
  geom_text(mapping = aes(label = format(round(pct, 1), nsmall = 1)), data = d4, vjust=-0.5, color="black",position = position_dodge(0.87), size=5.5)+
  scale_fill_hue(c = 40)


p3 + scale_y_continuous(name = "Reads", limits = c(0,2400000), breaks = c(400000, 800000, 1200000, 1600000,2000000,2400000)) +
  theme(legend.title=element_blank(),
        axis.ticks = element_line(size = 1),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0.3, "cm"),
        legend.box.background = element_rect(),
        legend.text = element_text(face="bold", color="black", size = 12),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size=16,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=14),
        axis.text.y = element_text(face="bold", color="black", size=14),
        panel.grid.major.y = element_line(color = "darkgray", linetype = "dashed"),
        axis.line = element_line(colour = "black", size = 0.8),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),)


# p3 <- ggplot(d4, aes(x = Sample, y = as.integer(value), fill = Type)) +
#   geom_bar(stat="identity",position=position_dodge()) +
#   geom_text(mapping = aes(label = format(round(pct, 1), nsmall = 1)), data = d4, vjust=-0.5, color="black",position = position_dodge(0.87), size=5.5)+
#   theme_minimal()
# p3 <- p3 + scale_fill_brewer(palette="Reds") +
#   scale_y_continuous(name = "Reads", limits = c(0,2400000), breaks = c(400000, 800000, 1200000, 1600000,2000000,2400000) ) +
#   ggtitle("LV641 15days, read distribution") +
#   theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title=element_text(size=14)) 

#######################################################################################################################################3


figure <- ggarrange(p, p2, p3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure
