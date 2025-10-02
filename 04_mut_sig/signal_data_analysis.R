# Load required libraries
library(tidyverse)
library(fs)
library(ggplot2)

# Set the base directory path
base_dir <- "C:/Stuff/uni_work/project/mibc_landscape/04_mut_sig/signal_download/MIBC_landscape_summer"

# Function to read contributions for a single sample (signature set 1 only)
read_sample_contributions <- function(sample_dir) {
  # Only read from signature set 1
  contrib_file <- file.path(sample_dir, "sbs_signatures", "1", "contributions.csv")
  
  if (file_exists(contrib_file)) {
    # Read the contributions file
    contrib_data <- read_csv(contrib_file, show_col_types = FALSE) %>%
      mutate(sample_id = basename(sample_dir))
    return(contrib_data)
  } else {
    warning(paste("No contributions.csv found for sample:", basename(sample_dir)))
    return(NULL)
  }
}

# Get all sample directories
sample_dirs <- dir_ls(base_dir, type = "directory", regexp = "TCGA-")

all_contributions <- map_dfr(sample_dirs, read_sample_contributions, .progress = TRUE)

# Reshape the data for easier analysis
contributions_wide <- all_contributions %>%
  select(sample_id, Signature, Contribution) %>%
  pivot_wider(
    names_from = Signature,
    values_from = Contribution,
    values_fill = 0
  )

# Plot signature contributions heatmap
sig_data <- all_contributions %>%
  filter(Signature != "Unassigned") %>%
  select(sample_id, Signature, Contribution) %>%
  pivot_wider(names_from = Signature, values_from = Contribution, values_fill = 0)

non_zero_cols <- sig_data %>%
  select(-sample_id) %>%
  summarise(across(everything(), ~ sum(.x) > 0)) %>%
  unlist()

sig_data <- sig_data %>%
  select(sample_id, all_of(names(non_zero_cols)[non_zero_cols]))


# plot --------------------------------------------------------------------

sig_data_long_all <- sig_data %>%
  pivot_longer(cols = -sample_id, 
               names_to = "Signature", 
               values_to = "Contribution")

ggplot(sig_data_long_all, aes(x = Signature, y = Contribution + 1)) +
    geom_boxplot(aes(fill = Signature), alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.35, alpha = 0.25, size = 0.5) +
    scale_y_log10() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = "Mutational Signature Contributions (Log Scale)",
      x = "Signature",
      y = "Contribution + 1 (log scale)"
    ) +
    scale_fill_viridis_d()

# stats -------------------------------------------------------------------

sig_stat <- sig_data |> summarise(across(starts_with("SBS"),
                                         list(median = median, sd = sd),
                                         .names = "{.col}_{.fn}")) |> 
  pivot_longer(
    everything(),
    names_to = c("signature", "statistic"),
    names_sep = "_",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = "signature",
    values_from = "value"
  ) |> column_to_rownames("statistic")


