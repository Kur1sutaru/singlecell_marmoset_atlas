library(tidyverse)

# Load data
mouse_prop <- read.csv("mouse_cell_count_proportion.csv")
mouse_count <- read.csv("mouse_cell_count_sex.csv")
marmoset_prop <- read.csv("marmoset_cell_proportion.csv")
marmoset_count <- read.csv("marmoset_cells_count.csv")

# ----------------------------
# 1. Mouse Female vs Male - Proportions
# ----------------------------
mouse_prop_long <- mouse_prop %>%
  pivot_longer(cols = c(FC, FC2, FC3, MC, MC2, MC3),
               names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sex = ifelse(grepl("^F", Sample), "Female", "Male"))

ggplot(mouse_prop_long, aes(x = id, y = Proportion, fill = Sex)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Mouse Female vs Male - Cell Proportions",
       x = "Cell Type", y = "Proportion")

# ----------------------------
# 2. Mouse Female vs Male - Counts
# ----------------------------
mouse_count_long <- mouse_count %>%
  pivot_longer(cols = c(Female_count, Male_count),
               names_to = "Sex", values_to = "Count") %>%
  mutate(Sex = ifelse(Sex == "Female_count", "Female", "Male"))

ggplot(mouse_count_long, aes(x = id, y = Count, fill = Sex)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Mouse Female vs Male - Cell Counts",
       x = "Cell Type", y = "Count")

# ----------------------------
# 3. Mouse vs Marmoset - Mean Proportions
# ----------------------------
avg_mouse <- mouse_prop %>%
  mutate(Mouse_mean = rowMeans(select(., FC, FC2, FC3, MC, MC2, MC3))) %>%
  select(id, Mouse_mean) %>%
  mutate(Species = "Mouse", Proportion = Mouse_mean) %>%
  select(id, Species, Proportion)

avg_marmoset <- marmoset_prop %>%
  mutate(Marmoset_mean = rowMeans(select(., marmoset, marmoset2))) %>%
  select(id, Marmoset_mean) %>%
  mutate(Species = "Marmoset", Proportion = Marmoset_mean) %>%
  select(id, Species, Proportion)

combined_prop <- bind_rows(avg_mouse, avg_marmoset)

ggplot(combined_prop, aes(x = id, y = Proportion, fill = Species)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Mouse vs Marmoset - Mean Cell Proportions",
       x = "Cell Type", y = "Proportion")

# ----------------------------
# 4. Example statistical tests (t-test for proportions)
# ----------------------------
stats_mouse_sex <- mouse_prop_long %>%
  group_by(id) %>%
  summarise(p_value = t.test(Proportion ~ Sex)$p.value) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

stats_species <- combined_prop %>%
  group_by(id) %>%
  summarise(p_value = t.test(Proportion ~ Species)$p.value) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

########### optional #########
library(ggpubr)

# Mouse Female vs Male - Proportions
mouse_prop_long <- mouse_prop %>%
  pivot_longer(cols = c(FC, FC2, FC3, MC, MC2, MC3),
               names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sex = ifelse(grepl("^F", Sample), "Female", "Male"))

ggplot(mouse_prop_long, aes(x = Sex, y = Proportion, fill = Sex)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ id, scales = "free_y") +
  theme_bw() +
  stat_compare_means(method = "t.test", label = "p.signif") +
  labs(title = "Mouse Female vs Male - Cell Proportions",
       x = "Sex", y = "Proportion")

# Mouse vs Marmoset - Proportions
avg_mouse <- mouse_prop %>%
  mutate(Proportion = rowMeans(select(., FC, FC2, FC3, MC, MC2, MC3)),
         Species = "Mouse") %>%
  select(id, Species, Proportion)

avg_marmoset <- marmoset_prop %>%
  mutate(Proportion = rowMeans(select(., marmoset, marmoset2)),
         Species = "Marmoset") %>%
  select(id, Species, Proportion)

combined_prop <- bind_rows(avg_mouse, avg_marmoset)

# Filter out ids without both species
combined_prop_filtered <- combined_prop %>%
  group_by(id) %>%
  filter(n_distinct(Species) == 2)

ggplot(combined_prop_filtered, aes(x = Species, y = Proportion, fill = Species)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ id, scales = "free_y") +
  theme_bw() +
  stat_compare_means(method = "t.test", label = "p.signif") +
  labs(title = "Mouse vs Marmoset - Mean Cell Proportions",
       x = "Species", y = "Proportion")

write.csv(marmoset_count, "marmosetcount.csv")
write.csv(mouse_count, "mousecounts.csv")

# Save stats tables
write.csv(stats_mouse_sex, "mouse_female_vs_male_stats.csv", row.names = FALSE)
write.csv(stats_species, "mouse_vs_marmoset_stats.csv", row.names = FALSE)

## test 1
library(tidyverse)

# Step 1: Prepare mouse replicate proportions
mouse_long <- mouse_prop %>%
  pivot_longer(cols = c(FC, FC2, FC3, MC, MC2, MC3),
               names_to = "Sample", values_to = "Proportion") %>%
  mutate(Species = "Mouse")

# Step 2: Prepare marmoset replicate proportions
marmoset_long <- marmoset_prop %>%
  pivot_longer(cols = c(marmoset, marmoset2),
               names_to = "Sample", values_to = "Proportion") %>%
  mutate(Species = "Marmoset")

# Step 3: Combine and filter to shared cell types
shared_celltypes <- c(
  "Basal keratinocyte",
  "Lymphatic endothelial",
  "Macrophages",
  "Myelinating Schwann cells",
  "Pericytes",
  "Unmyelinating Schwann cells"
)

all_long <- bind_rows(mouse_long, marmoset_long) %>%
  filter(id %in% shared_celltypes)

# Step 4: Aggregate proportion per replicate (Sample)
bar_data <- all_long %>%
  group_by(Species, Sample, id) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop")

# Step 5: Stacked barplot
ggplot(bar_data, aes(x = Sample, y = Proportion, fill = id)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Stacked Barplot of Cell Type Proportions",
    x = "Sample", y = "Proportion",
    fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###### Mean ± SD Bar Plot per Cell Type and Species

library(tidyverse)

# Prepare long-format data
mouse_long <- mouse_prop %>%
  pivot_longer(cols = c(FC, FC2, FC3, MC, MC2, MC3),
               names_to = "Sample", values_to = "Proportion") %>%
  mutate(Species = "Mouse")

marmoset_long <- marmoset_prop %>%
  pivot_longer(cols = c(marmoset, marmoset2),
               names_to = "Sample", values_to = "Proportion") %>%
  mutate(Species = "Marmoset")

shared_celltypes <- c(
  "Basal keratinocyte",
  "Lymphatic endothelial",
  "Macrophages",
  "Myelinating Schwann cells",
  "Pericytes",
  "Unmyelinating Schwann cells"
)

# Combine and filter
all_long <- bind_rows(mouse_long, marmoset_long) %>%
  filter(id %in% shared_celltypes)

# Calculate mean and SD per cell type and species
summary_df <- all_long %>%
  group_by(Species, id) %>%
  summarise(
    Mean = mean(Proportion, na.rm = TRUE),
    SD = sd(Proportion, na.rm = TRUE),
    .groups = "drop"
  )

# Plot mean ± SD
ggplot(summary_df, aes(x = id, y = Mean, fill = Species)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mean ± SD of Cell Type Proportions",
    x = "Cell Type", y = "Proportion", fill = "Species"
  )

#### mouse female x male
library(tidyverse)

# Step 1: Long format for all mouse replicates
mouse_long <- mouse_prop %>%
  pivot_longer(cols = c(FC, FC2, FC3, MC, MC2, MC3),
               names_to = "Sample", values_to = "Proportion") %>%
  mutate(
    Sex = case_when(
      str_starts(Sample, "F") ~ "Female",
      str_starts(Sample, "M") ~ "Male",
      TRUE ~ NA_character_
    )
  )

# Step 2: Summarize proportions for each sample and cell type
bar_data_all <- mouse_long %>%
  group_by(Sex, Sample, id) %>%
  summarise(Proportion = sum(Proportion, na.rm = TRUE), .groups = "drop")

# Step 3: Stacked barplot of all cell types per sample
ggplot(bar_data_all, aes(x = Sample, y = Proportion, fill = id)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Sex, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Stacked Barplot of All Cell Type Proportions - Mouse Sex Comparison",
    x = "Sample", y = "Proportion",
    fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### trying to add number of cells for each celltype
# Total per sample
library(tidyverse)

# Step 1: Prepare data
mouse_bar <- mouse_count %>%
  pivot_longer(cols = c(Female_count, Male_count),
               names_to = "Sex", values_to = "CellCount") %>%
  mutate(Sex = ifelse(Sex == "Female_count", "Female", "Male"))

# Step 2: Sort cell types by total abundance
celltype_order <- mouse_bar %>%
  group_by(id) %>%
  summarise(Total = sum(CellCount), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(id)

mouse_bar <- mouse_bar %>%
  mutate(id = factor(id, levels = celltype_order))  # reorder by abundance

# Step 3: Compute proportions
mouse_bar <- mouse_bar %>%
  group_by(Sex) %>%
  mutate(Proportion = CellCount / sum(CellCount)) %>%
  ungroup()

# Step 4: Plot with cell count stacks and proportion labels
ggplot(mouse_bar, aes(x = Sex, y = CellCount, fill = id)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Proportion * 100), "%")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  theme_bw() +
  labs(
    title = "Mouse Male vs Female – Stacked Barplot by Cell Type",
    x = "Sex", y = "Number of Cells",
    fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

### trying a sankey plot
library(tidyverse)
library(networkD3)

# Step 1: Long format
df_sankey <- mouse_count %>%
  pivot_longer(cols = c(Female_count, Male_count),
               names_to = "Sex", values_to = "Count") %>%
  mutate(Sex = ifelse(Sex == "Female_count", "Female", "Male"))

# Step 2: Create node list
nodes <- data.frame(name = c(unique(df_sankey$Sex), unique(df_sankey$id)))

# Step 3: Create link list
df_sankey_links <- df_sankey %>%
  mutate(
    source = match(Sex, nodes$name) - 1,
    target = match(id, nodes$name) - 1,
    value = Count
  ) %>%
  select(source, target, value)

# Step 4: Plot Sankey
sankeyNetwork(Links = df_sankey_links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize = 12, nodeWidth = 30)

### plan B
library(tidyverse)
library(networkD3)   # for sankey plot
library(ggplot2)
# Assume mouse_count has: id, Female_count, Male_count
df_sankey <- mouse_count %>%
  pivot_longer(cols = c(Female_count, Male_count),
               names_to = "Sex", values_to = "Count") %>%
  mutate(Sex = ifelse(Sex == "Female_count", "Female", "Male"))

# Create node list
nodes <- data.frame(name = c(unique(df_sankey$Sex), unique(df_sankey$id)))

# Create link list (indexing starts at 0)
df_links <- df_sankey %>%
  mutate(
    source = match(Sex, nodes$name) - 1,
    target = match(id, nodes$name) - 1,
    value = Count
  ) %>%
  select(source, target, value)

# Plot sankey
sankeyNetwork(Links = df_links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize = 13, nodeWidth = 30)
# Prepare data
df_bar <- mouse_count %>%
  pivot_longer(cols = c(Female_count, Male_count),
               names_to = "Sex", values_to = "Count") %>%
  mutate(Sex = ifelse(Sex == "Female_count", "Female", "Male"))

# Sort by total abundance
cell_order <- df_bar %>%
  group_by(id) %>%
  summarise(total = sum(Count)) %>%
  arrange(desc(total)) %>%
  pull(id)

df_bar <- df_bar %>%
  mutate(id = factor(id, levels = cell_order))

# Plot
ggplot(df_bar, aes(x = id, y = Count, fill = Sex)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(title = "Male vs Female – Cell Counts per Type",
       x = "Cell Type", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
