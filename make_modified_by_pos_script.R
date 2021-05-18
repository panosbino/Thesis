library(tidyverse)
library(reshape2)
library(scales)



data <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_15days_Cas9_/Effect_vector_combined.txt")
data2 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_15days_Cas9-/Effect_vector_combined.txt")

comb_data <- cbind(data, data2$effect)



colnames(comb_data) <- c("Position", "Cas9+", "Cas9-")

d <- melt(comb_data, id.vars="Position")
colnames(d) <- c("Position", "Sample", "Percentage_of_modified_alleles")
p <- ggplot(d, aes(x = Position, y= Percentage_of_modified_alleles/100, col = Sample )) +
  geom_rect(aes(xmin=48, xmax=71, ymin=0, ymax=0.6), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=75, xmax=98, ymin=0, ymax=0.6), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=102, xmax=125, ymin=0, ymax=0.6), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=129, xmax=152, ymin=0, ymax=0.6), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=156, xmax=179, ymin=0, ymax=0.6), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  # geom_vline(xintercept = 48, linetype = "dashed") +
  # geom_vline(xintercept = 71, linetype = "dashed") +
  # geom_vline(xintercept = 75, linetype = "dashed") +
  # geom_vline(xintercept = 98, linetype = "dashed") +
  # geom_vline(xintercept = 102, linetype = "dashed") +
  # geom_vline(xintercept = 125, linetype = "dashed") +
  # geom_vline(xintercept = 129, linetype = "dashed") +
  # geom_vline(xintercept = 152, linetype = "dashed") +
  # geom_vline(xintercept = 156, linetype = "dashed") +
  # geom_vline(xintercept = 179, linetype = "dashed") +
  geom_line(lwd = 1)


#breaks = c(48,71,75,98,102,125,129,152,156,179)
p + scale_x_continuous(name="Position", limits=c(40, 200), breaks = c() ) +
  scale_y_continuous(name="Percentage of modified alleles (%)", limits=c(0, 0.6), labels = percent, breaks = c(0,0.2,0.4,0.6)) + 
  theme(plot.title = element_text(face = "bold", size = 18, color = "black", hjust = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.8),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.text = element_blank(),
        #axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.title.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.y = element_text(face = "bold", color = "black", size = 14),
        legend.background = element_rect(),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.text = element_text(face = "bold", color = "black", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "black", linetype = "dashed", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.ontop = TRUE)



############################################################ 3_days ######################################################################3


data <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_3days_Cas9_/Effect_vector_combined.txt")
data2 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV641_3days_Cas9-/Effect_vector_combined.txt")

comb_data <- cbind(data, data2$effect)



colnames(comb_data) <- c("Position", "Cas9+", "Cas9-")

d <- melt(comb_data, id.vars="Position")
colnames(d) <- c("Position", "Sample", "Percentage_of_modified_alleles")

p2 <- ggplot(d, aes(x = Position, y= Percentage_of_modified_alleles/100, col = Sample )) +
  geom_rect(aes(xmin=48, xmax=71, ymin=0, ymax=0.15), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=75, xmax=98, ymin=0, ymax=0.15), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=102, xmax=125, ymin=0, ymax=0.15), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=129, xmax=152, ymin=0, ymax=0.15), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  geom_rect(aes(xmin=156, xmax=179, ymin=0, ymax=0.15), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  # geom_vline(xintercept = 48, linetype = "dashed") +
  # geom_vline(xintercept = 71, linetype = "dashed") +
  # geom_vline(xintercept = 75, linetype = "dashed") +
  # geom_vline(xintercept = 98, linetype = "dashed") +
  # geom_vline(xintercept = 102, linetype = "dashed") +
  # geom_vline(xintercept = 125, linetype = "dashed") +
  # geom_vline(xintercept = 129, linetype = "dashed") +
  # geom_vline(xintercept = 152, linetype = "dashed") +
  # geom_vline(xintercept = 156, linetype = "dashed") +
  # geom_vline(xintercept = 179, linetype = "dashed") +
  geom_line(lwd = 1)


#breaks = c(48,71,75,98,102,125,129,152,156,179)
p2 + scale_x_continuous(name="Position", limits=c(40, 200), breaks = c() ) +
  scale_y_continuous(name="Percentage of modified alleles (%)", limits=c(0, 0.15), labels = percent) +
  theme(plot.title = element_text(face = "bold", size = 18, color = "black", hjust = 0.5),
        axis.text = element_blank(),
        #axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        axis.line.y = element_line(colour = "black", size = 0.8),
        axis.ticks.length.y = unit(0.3, "cm"),
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.title.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.y = element_text(face = "bold", color = "black", size = 14),
        legend.background = element_rect(),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.text = element_text(face = "bold", color = "black", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "black", linetype = "dashed", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.ontop = TRUE)
# p <- ggplot(d, aes(x = Position, y= Percentage_of_modified_alleles/100, col = Sample )) +
#   geom_rect(aes(xmin=48, xmax=71, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_rect(aes(xmin=75, xmax=98, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_rect(aes(xmin=102, xmax=125, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_rect(aes(xmin=129, xmax=152, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_rect(aes(xmin=156, xmax=179, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_vline(xintercept = 48, linetype = "dashed") +
#   geom_vline(xintercept = 71, linetype = "dashed") +
#   geom_vline(xintercept = 75, linetype = "dashed") +
#   geom_vline(xintercept = 98, linetype = "dashed") +
#   geom_vline(xintercept = 102, linetype = "dashed") +
#   geom_vline(xintercept = 125, linetype = "dashed") +
#   geom_vline(xintercept = 129, linetype = "dashed") +
#   geom_vline(xintercept = 152, linetype = "dashed") +
#   geom_vline(xintercept = 156, linetype = "dashed") +
#   geom_vline(xintercept = 179, linetype = "dashed") +
#   geom_line()
# 
# p + scale_x_continuous(name="Position", limits=c(0, 200), breaks = c(48,71,75,98,102,125,129,152,156,179)) +
#   scale_y_continuous(name="Percentage of modified alleles (%)", limits=c(0, 0.15), labels = percent) + 
#   ggtitle("Percentage of modified alleles by position for LV641 \n3 days post transduction") +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(size = 16, hjust = 0.5))
# 

############################################################ LV678 #######################################################################

data <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV678_3days_Cas9_/Effect_vector_combined.txt")
data2 <- read_tsv("/home/panosbino/Desktop/Thesis/Results/In_vitro_amp130-148/pooled_crispresso/CRISPResso_on_LV678_3days_Cas9-/Effect_vector_combined.txt")

comb_data <- cbind(data, data2$effect)



colnames(comb_data) <- c("Position", "Cas9+", "Cas9-")

d <- melt(comb_data, id.vars="Position")
colnames(d) <- c("Position", "Sample", "Percentage_of_modified_alleles")

p3 <- ggplot(d, aes(x = Position, y= Percentage_of_modified_alleles/100, col = Sample )) +
  geom_rect(aes(xmin=49, xmax=72, ymin=0, ymax=0.05), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
  geom_rect(aes(xmin=95, xmax=119, ymin=0, ymax=0.05), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
  geom_rect(aes(xmin=142, xmax=165, ymin=0, ymax=0.05), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
  # geom_rect(aes(xmin=48, xmax=71, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  # geom_rect(aes(xmin=75, xmax=98, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  # geom_rect(aes(xmin=102, xmax=125, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  # geom_rect(aes(xmin=129, xmax=152, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  # geom_rect(aes(xmin=156, xmax=179, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.5) +
  # # geom_vline(xintercept = 48, linetype = "dashed") +
  # geom_vline(xintercept = 71, linetype = "dashed") +
  # geom_vline(xintercept = 75, linetype = "dashed") +
  # geom_vline(xintercept = 98, linetype = "dashed") +
  # geom_vline(xintercept = 102, linetype = "dashed") +
  # geom_vline(xintercept = 125, linetype = "dashed") +
  # geom_vline(xintercept = 129, linetype = "dashed") +
  # geom_vline(xintercept = 152, linetype = "dashed") +
  # geom_vline(xintercept = 156, linetype = "dashed") +
  # geom_vline(xintercept = 179, linetype = "dashed") +
  geom_line(lwd = 1)


#breaks = c(48,71,75,98,102,125,129,152,156,179)
p3 + scale_x_continuous(name="Position", limits=c(40, 200), breaks = c() ) +
  scale_y_continuous(name="Percentage of modified alleles (%)", limits=c(0, 0.05), labels = percent) + 
  theme(plot.title = element_text(face = "bold", size = 18, color = "black", hjust = 0.5),
        axis.text = element_blank(),
        axis.line.y = element_line(colour = "black", size = 0.8),
        axis.ticks.length.y = unit(0.3, "cm"),
        #axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        axis.title.x = element_text(face = "bold", color = "black", size = 14),
        axis.title.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.y = element_text(face = "bold", color = "black", size = 14),
        legend.background = element_rect(),
        legend.title = element_blank(),
        legend.box.background = element_rect(),
        legend.text = element_text(face = "bold", color = "black", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "black", linetype = "dashed", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.ontop = TRUE)


# p <- ggplot(d, aes(x = Position, y= Percentage_of_modified_alleles/100, col = Sample )) +
#   geom_rect(aes(xmin=49, xmax=72, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_rect(aes(xmin=95, xmax=119, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_rect(aes(xmin=142, xmax=165, ymin=0, ymax=Inf), color= "lightgrey", fill = "lightgrey", alpha = 0.9) +
#   geom_vline(xintercept = 49, linetype = "dashed") +
#   geom_vline(xintercept = 72, linetype = "dashed") +
#   geom_vline(xintercept = 95, linetype = "dashed") +
#   geom_vline(xintercept = 119, linetype = "dashed") +
#   geom_vline(xintercept = 142, linetype = "dashed") +
#   geom_vline(xintercept = 165, linetype = "dashed") +
#   geom_line()
# 
# p + scale_x_continuous(name="Position", limits=c(0, 200), breaks = c(49,72,95,119,142,165)) +
#   scale_y_continuous(name="Percentage of modified alleles (%)", limits=c(0, 0.05), labels = percent) + 
#   ggtitle("Percentage of modified alleles by position for LV678 \n3 days post transduction") +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(size = 16, hjust = 0.5))


