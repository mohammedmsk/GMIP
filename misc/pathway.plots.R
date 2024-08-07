data=read.csv('~/GMIP/results/all_in_one.pathways.txt', sep='\t')
dim(data)
data=data[grepl("lipid|liver", data$description, ignore.case=T), ]
dim(data)
table(data$List)
data$List=gsub("GWAS_WKP", "GWAS", data$List)
data$List=gsub("popsFeatures_GMIP_nc3_DIS", "popsFeatures_GMIP_nc3", data$List)
data$List=gsub("popsFeatures_GMIP_nc3_WKP", "popsFeatures_GMIP_nc3", data$List)
data$List=gsub("GSE166504_GMIP_nc3_DIS", "GSE166504_GMIP_nc3", data$List)
table(data$List)

table(data$description)
data$description=gsub("Familial hyperlipidemia type 1", "Familial hyperlipidemia", data$description)
data$description=gsub("Familial hyperlipidemia type 2", "Familial hyperlipidemia", data$description)
data$description=gsub("Familial hyperlipidemia type 3", "Familial hyperlipidemia", data$description)
data$description=gsub("Familial hyperlipidemia type 4", "Familial hyperlipidemia", data$description)
data$description=gsub("Familial hyperlipidemia type 5", "Familial hyperlipidemia", data$description)
table(data$description)
data$description=gsub("Fatty Liver", "Non-alcoholic fatty liver disease", data$description)
data$description=gsub("Nonalcoholic fatty liver disease", "Non-alcoholic fatty liver disease", data$description)
data$description=gsub("Non-alcoholic Non-alcoholic fatty liver disease Disease", "Non-alcoholic fatty liver disease", data$description)
table(data$description)
data=data[!grepl("Sphingolipid", data$description, ignore.case=T), ]
table(data$description)
data$description=gsub("Liver neoplasms", "Liver carcinoma", data$description)
data$description=gsub("Liver Neoplasms, Experimental", "Liver carcinoma", data$description)
data$description=gsub("Liver Cirrhosis, Alcoholic", "Liver Cirrhosis", data$description)
data$description=gsub("Liver Cirrhosis, Experimental", "Liver Cirrhosis", data$description)
data$description=gsub("Lipid metabolism in senescent cells", "Lipid metabolism & pathway related", data$description)
data$description=gsub("Nuclear receptors in lipid metabolism and toxicity", "Lipid metabolism & pathway related", data$description)
data$description=gsub("SREBF and miR33 in cholesterol and lipid homeostasis", "Lipid metabolism & pathway related", data$description)
data$description=gsub("Lipid metabolism pathway", "Lipid metabolism & pathway related", data$description)
table(data$description)
<<<<<<< HEAD
table(data$List)
=======
data.table::fwrite(data, "/gstore/scratch/u/kanchwam/gp/05577fee5848687884crrerfefe845/results_2/pathway_analysis/all_in_one.pathways.filtered_for_fig.txt", sep="\t")
>>>>>>> 59b15dc92768b2c2a54f6de5b1d27cdf0b2cb28b
#######################################################################################################################
library(dplyr)
library(tidyr)

# Ensure the userId column is treated as strings and split it into individual genes
data <- data %>% mutate(userId = strsplit(as.character(userId), ";"))

# Prepare the links data for List to Description
links1 <- data %>%
  select(source = List, target = description) %>%
  group_by(source, target) %>%
  summarize(value = n()) %>%
  ungroup()

# Prepare the links data for Description to Genes (userId)
links2 <- data %>%
  unnest(userId) %>%
  select(source = description, target = userId) %>%
  group_by(source, target) %>%
  summarize(value = n()) %>%
  ungroup()

# Combine both link datasets
links <- bind_rows(links1, links2)

# Function to create flow string
create_flow_string <- function(source, target, value, color = NULL) {
  if (is.null(color)) {
    return(paste(source, "[", value, "]", target))
  } else {
    return(paste(source, "[", value, "]", target, color))
  }
}

# Create flow strings
flow_strings <- links %>%
  rowwise() %>%
  mutate(flow_string = create_flow_string(source, target, value)) %>%
  pull(flow_string)

# Print flow strings
sink("flow_strings.txt")
cat(flow_strings, sep = "\n")
sink()

# Define node colors
node_colors <- data.frame(
  node = c("List", "Description", "Genes"),
  color = c("#1f77b4", "#ff7f0e", "#2ca02c")
)

# Create node color strings
node_color_strings <- node_colors %>%
  rowwise() %>%
  mutate(color_string = paste0(":", node, " ", color)) %>%
  pull(color_string)

# Print node color strings
cat(node_color_strings, sep = "\n")

#######################################################################################################################

# Ensure the userId column is treated as strings and split it into individual genes
data <- data %>% mutate(userId = strsplit(as.character(userId), ";"))

# Prepare the nodes data
nodes <- data.frame(name = unique(c(data$List, data$description, unlist(data$userId))))

# Add a column to classify the type of node
nodes$type <- ifelse(nodes$name %in% data$List, "List",
                     ifelse(nodes$name %in% data$description, "Description", "Genes"))

# Prepare the links data for List to Description
links1 <- data %>%
  select(source = List, target = description) %>%
  group_by(source, target) %>%
  summarize(value = n()) %>%
  ungroup()

# Prepare the links data for Description to Genes (userId)
links2 <- data %>%
  unnest(userId) %>%
  select(source = description, target = userId) %>%
  group_by(source, target) %>%
  summarize(value = n()) %>%
  ungroup()

# Combine both link datasets
links <- bind_rows(links1, links2)

# Map source and target to their node indexes
links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1

# Convert links to a plain data frame
links <- as.data.frame(links)

# Custom colors for nodes
my_color <- 'd3.scaleOrdinal()
                .domain(["List", "Description", "Genes"])
                .range(["#1f77b4", "#ff7f0e", "#2ca02c"])'

# Create the Sankey diagram with custom colors
sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                        Value = "value", NodeID = "name", units = "Genes", fontSize = 12,
                        nodeWidth = 30, colourScale = my_color, NodeGroup = "type")

# Display the Sankey diagram
sankey

