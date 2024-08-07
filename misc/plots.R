setwd("/home1/04179/mshabb/GMIP/results/")
library(data.table)
library(ggplot2)
library(ggrepel) # For better text label placement
library(RColorBrewer)
library(reshape2)

#Read in the resulst files
res_dt=data.table::fread("/home1/04179/mshabb/GMIP/results/2024_04_20_01_GMIP_noFeatureClustering.final.full.report.tsv.gz")
#Clean up and get methods
res_dt$Method=gsub(".*_method_", "", res_dt$Trait)
res_dt$Method=gsub(".gmip_withoutFC.*", ".gmip_withoutFC", res_dt$Method)
dim(res_dt)
res_dt=res_dt[res_dt$Method!="PLSRegressionCV.gmip_withoutFC", ]
dim(res_dt)
res_dt=res_dt[res_dt$Method!="PLSRegressionCV.pops", ]
dim(res_dt)
res_dt=res_dt[res_dt$Method!="PC_RidgeCV.gmip_withoutFC", ]
dim(res_dt)
res_dt=res_dt[res_dt$Method!="PC_RidgeCV.pops", ]
dim(res_dt)
res_dt=res_dt[res_dt$Method!="PC_RidgeCV.gmip_withFC_PC_RidgeCV", ]
dim(res_dt)

table(res_dt$Method)

res_dt$Method=gsub("PLSRegression_nc1.gmip_withoutFC", "GMIP_PLSR_nc1", res_dt$Method)
res_dt$Method=gsub("PLSRegression_nc2.gmip_withoutFC", "GMIP_PLSR_nc2", res_dt$Method)
res_dt$Method=gsub("PLSRegression_nc3.gmip_withoutFC", "GMIP_PLSR_nc3", res_dt$Method)
res_dt$Method=gsub("PLSRegression_nc5.gmip_withoutFC", "GMIP_PLSR_nc5", res_dt$Method)
res_dt$Method=gsub("PLSRegression_nc10.gmip_withoutFC", "GMIP_PLSR_nc10", res_dt$Method)
res_dt$Method=gsub("pops", "PoPS", res_dt$Method)

table(res_dt$Method)
#Cleanup and get trait names
res_dt$Trait2=gsub("_collapsing_.*", "", res_dt$Trait)
res_dt=res_dt[res_dt$Trait2!="Total_bilirubin__umol.L__Code30840", ]
res_dt=res_dt[res_dt$Trait2!="Direct_bilirubin__umol.L__Code30660", ]
res_dt=res_dt[res_dt$Trait2!="Rheumatoid_factor__IU.ml__Code30820", ]
dim(res_dt)
#Cleanup and get collapsing info
res_dt$collapsing=gsub("_method_.*", "", gsub(".*_collapsing_", "", res_dt$Trait))
table(res_dt$collapsing)
res_dt=res_dt[res_dt$collapsing!="true", ]
dim(res_dt)
#Cleanup the Label
res_dt$Label=gsub("pops_features_top", "", res_dt$Label)
#Correct all pops names
res_dt[grepl("PoPS", res_dt$Method),]$Method="PoPS"
table(res_dt$Method)
#Create a ID col
res_dt$ID=paste0(res_dt$Trait2, "__", res_dt$Method, "__", res_dt$Label)
#Subset to have a minimal info data.table
res_dt_minimal=res_dt[,c("Trait2", "Method", "Label", "tau_normalized", "tau_standard_error_normalized", "Enrichment_p", "Prop._h2")]
dim(res_dt_minimal)
#Remove repeating pops duplicates
res_dt_minimal=unique(res_dt_minimal)
dim(res_dt_minimal)
table(res_dt_minimal$Method)
length(table(res_dt$Trait2))
dim(res_dt_minimal)
head(res_dt_minimal)
table(res_dt_minimal$Method)

fwrite(data.table(GWAS_list=unique(res_dt_minimal$Trait2)), "GWAS_list_used.txt")
#Calculate p-adj using bonferroni correction
res_dt_minimal$padj=p.adjust(res_dt_minimal$Enrichment_p, method="bonferroni")
dim(res_dt_minimal)
head(res_dt_minimal)

res_dt_minimal_filt=res_dt_minimal[res_dt_minimal$padj<=0.01 & res_dt_minimal$tau_normalized >=0,]
dim(res_dt_minimal_filt)
length(table(res_dt_minimal_filt$Trait2))

# Collapse data to retain the row with the largest tau_normalized for each Trait2
collapsed_res_dt_minimal_filt=res_dt_minimal_filt[, .SD[which.max(tau_normalized)], by = Trait2]
dim(collapsed_res_dt_minimal_filt)
table(collapsed_res_dt_minimal_filt$Label)
table(collapsed_res_dt_minimal_filt$Method)
data.table::fwrite(collapsed_res_dt_minimal_filt,"filtered_collapsed.2024_04_20_01_GMIP_noFeatureClustering.final.full.report.tsv.gz")
#collapsed_res_dt_minimal_filt2=res_dt_minimal_filt[, .SD[which.max(tau_normalized)], by = .(Trait2, Method)]
########################################################################################################################
# Convert Label to a factor for color representation
collapsed_res_dt_minimal_filt[, Label := as.factor(Label)]

# Add a small constant to tau_normalized values and then round to 2 decimal places
collapsed_res_dt_minimal_filt[, tau_normalized := round(tau_normalized + 1e-5, 2)]

# Order Trait2 by tau_normalized
collapsed_res_dt_minimal_filt[, Trait2 := reorder(Trait2, Prop._h2, median)]

# Get Set1 palette and add a 10th color (black in this case)
set1_colors <- brewer.pal(n = 9, name = "Set1")
extended_colors <- c(set1_colors, "black")

c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "black")

method_order <- rev(c("PoPS", "GMIP_PLSR_nc1", "GMIP_PLSR_nc2", "GMIP_PLSR_nc3", "GMIP_PLSR_nc5", "GMIP_PLSR_nc10"))
collapsed_res_dt_minimal_filt$Method <- factor(collapsed_res_dt_minimal_filt$Method, levels = method_order)

# Assuming `collapsed_res_dt_minimal_filt` is already your data.table object
g1 <- ggplot(collapsed_res_dt_minimal_filt, aes(x = Trait2, y = Method, size = tau_normalized, color = Label)) +
  geom_point() +
  geom_point(data = subset(collapsed_res_dt_minimal_filt, Label == "0500"), shape = 21, color = "black", stroke = 1.5) +
  scale_size_continuous(
    name = "Tau Normalized",
    range = c(3, 10),
    breaks = seq(min(collapsed_res_dt_minimal_filt$tau_normalized),
                 max(collapsed_res_dt_minimal_filt$tau_normalized),
                 length.out = 5)) +
  scale_color_manual(
    values = extended_colors,
    guide = guide_legend(override.aes = list(size = 6))  # Increase size of colored dots in legend
  ) +
  theme_minimal() +
  labs(
    title = "Dot Plot of Traits by Method",
    x = "Trait",
    y = "Method",
    color = "Label"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, face = "bold")  # Increase legend title size
  )

# Print the plot
print(g1)

ggsave(plot = g1, "dotplot_figure2.pdf", width = 15, height = 7)
ggsave(plot = g1, "dotplot_figure2.png", width = 15, height = 7)

g2=g1 + coord_flip() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, face = "bold")  # Increase legend title size
  )

ggsave(plot = g2, "dotplot_figure3_vertical.pdf", width = 12, height = 12)
ggsave(plot = g2, "dotplot_figure3_vertical.png", width = 12, height = 12)

method_order <- c("PoPS", "GMIP_PLSR_nc1", "GMIP_PLSR_nc2", "GMIP_PLSR_nc3", "GMIP_PLSR_nc5", "GMIP_PLSR_nc10")
collapsed_res_dt_minimal_filt$Method <- factor(collapsed_res_dt_minimal_filt$Method, levels = method_order)

# Count the actual counts for each label
label_counts <- collapsed_res_dt_minimal_filt[, .N, by = Label]

# Ensure the colors match the previous plot
names(extended_colors) <- levels(collapsed_res_dt_minimal_filt$Label)

# Calculate the maximum y-axis limit
max_y <- ceiling(max(label_counts$N))

# Create the horizontal bar plot with a border for the "0500" bar
g2=ggplot(label_counts, aes(x = Label, y = N, fill = Label)) +
  geom_bar(stat = "identity", aes(color = ifelse(Label == "0500", "darkred", NA)), size = 1.5) +
  geom_text(aes(label = N), hjust = -0.5, size = 12) + # Add values on top of the bars and increase text size
  scale_fill_manual(values = extended_colors) + # Use the extended color palette
  scale_color_identity() + # Use the color identity scale to apply the border color
  theme_minimal() +
  labs(x = "Number of top N Genes",
       y = "Count",
       fill = "topN") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 32, face = "bold"),
    axis.title.x = element_text(size = 32, face = "bold"),
    axis.title.y = element_text(size = 32, face = "bold"),
    axis.text.x = element_text(size = 32, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 32),  # Increase legend text size
    legend.title = element_text(size = 32, face = "bold")  # Increase legend title size
  ) + ylim(0, 16) +
  coord_flip() # Flip the coordinates for a horizontal bar plot

print(g2)
ggsave(plot=g2, "barplot1.pdf", width = 13, height = 8)
ggsave(plot=g2, "barplot1.png", width = 13, height = 8)


g2b <- ggplot(label_counts, aes(x = Label, y = N, fill = Label)) +
  geom_bar(stat = "identity", aes(color = ifelse(Label == "0500", "darkred", NA)), size = 1.5) +
  geom_text(aes(label = N), vjust = -0.5, size = 5) + # Adjust the text position for vertical bars
  scale_fill_manual(values = extended_colors) + # Use the extended color palette
  scale_color_identity() + # Use the color identity scale to apply the border color
  theme_minimal() +
  labs(
    title = "Count of Traits with Best Tau when Top N Genes Considered",
    x = "Number of Top N Genes Considered",
    y = "Count",
    fill = "Label"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, face = "bold")  # Increase legend title size
  )

# Print the plot
print(g2b)
ggsave(plot=g2b, "barplot1b.pdf", width = 14, height = 9)
ggsave(plot=g2b, "barplot1b.png", width = 14, height = 9)


# Convert Method to a factor for color representation (if not already done)
collapsed_res_dt_minimal_filt[, Method := as.factor(Method)]

# Count the actual counts for each Method
method_counts <- collapsed_res_dt_minimal_filt[, .N, by = Method]

# Create the horizontal bar plot with increased text sizes
g3=ggplot(method_counts, aes(x = Method, y = N, fill = Method)) +
  geom_bar(stat = "identity") + # Add border to all bars
  geom_text(aes(label = N), hjust = -0.5, size = 12) + # Add values on top of the bars and increase text size
  theme_minimal() +
  labs(
       x = "Method",
       y = "Count",
       fill = "Method") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 32, face = "bold"),
    axis.title.x = element_text(size = 32, face = "bold"),
    axis.title.y = element_text(size = 32, face = "bold"),
    axis.text.x = element_text(size = 32, angle = 90),
    axis.text.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 32),  # Increase legend text size
    legend.title = element_text(size = 32, face = "bold"),  # Increase legend title size,
    legend.position = "None"
  ) + ylim(0, 27) +
  coord_flip() # Flip the coordinates for a horizontal bar plot
print(g3)
ggsave(plot=g3, "barplot2.pdf", width = 14, height = 8)
ggsave(plot=g3, "barplot2.png", width = 14, height = 8)

library(ggplot2)

# Assuming `method_counts` is already your data.frame or data.table object
g3b <- ggplot(method_counts, aes(x = Method, y = N, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = N), vjust = -0.5, size = 5) + # Add values on top of the bars and increase text size
  theme_minimal() +
  labs(
    title = "Count of Methods",
    x = "Method",
    y = "Count",
    fill = "Method"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, face = "bold")  # Increase legend title size
  )

# Print the plot
print(g3b)
ggsave(plot=g3b, "barplot2b.pdf", width = 14, height = 9)
ggsave(plot=g3b, "barplot2b.png", width = 14, height = 9)

########################################################################################################################
df=res_dt_minimal[res_dt_minimal$Label=="0500",]
table(df$Method)

# Subset the data for the two methods
pops_data=df[df$Method == "PoPS"]
gmip_nc3_data=df[df$Method == "GMIP_PLSR_nc3"]

# Merge the two subsets on Trait2 and Label
merged_data <- merge(pops_data, gmip_nc3_data, by = "Trait2", suffixes = c("_pops", "_gmip_nc3"))

# Assuming `merged_data` is already your data.frame or data.table object
# Find the common limit for both axes
max_limit <- max(c(merged_data$tau_normalized_pops, merged_data$tau_normalized_gmip_nc3))

h1=ggplot(merged_data, aes(x = tau_normalized_pops, y = tau_normalized_gmip_nc3)) +
  geom_point(size = 3, color = 'blue', alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") + # Add diagonal line
  labs(
    title = "Comparison of Tau Values Between Methods",
    x = "Tau (pops)",
    y = "Tau (GMIP_PLSR_nc3)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, face = "bold")  # Increase legend title size
  )
print(h1)

library(ggplot2)

# Assuming `merged_data` is already your data.frame or data.table object
# Find the common limit for both axes
max_limit <- max(c(merged_data$tau_normalized_pops, merged_data$tau_normalized_gmip_nc3))

h2 <- ggplot(merged_data, aes(x = tau_normalized_pops, y = tau_normalized_gmip_nc3, size = -log10(padj_gmip_nc3))) +
  geom_point(color = "#984EA3", alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red") + # Add diagonal line
  labs(
    title = "Comparison of Tau Values Between Methods",
    x = "Tau normalized (pops)",
    y = "Tau normalized (GMIP_PLSR_nc3)",
    size = expression(-log[10](padj) ~ "(nc3)")  # Label for size legend
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, face = "bold")  # Increase legend title size
  ) +
  coord_cartesian(xlim = c(0, max_limit), ylim = c(0, max_limit))
# Print the plot
print(h2)
ggsave(plot=h2, "gmip_pops_scatter_top500.pdf", width = 12, height = 10)
ggsave(plot=h2, "gmip_pops_scatter_top500.png", width = 12, height = 10)
########################################################################################################################


