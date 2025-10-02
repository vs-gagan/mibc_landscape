## VAF Analysis -----------------------------------------------------------

# Calculate VAF for both datasets
wxs_trimmed[, VAF := t_alt_count / (t_ref_count + t_alt_count)]
wgs_trimmed[, VAF := t_alt_count / (t_ref_count + t_alt_count)]

# Create VAF comparison dataset
vaf_data <- rbind(
  wxs_trimmed[mut_id %in% wxs_only_ids, .(Category = "WXS_Unique", VAF, Hugo_Symbol, mut_id)],
  wgs_trimmed[mut_id %in% wgs_only_ids, .(Category = "WGS_Unique", VAF, Hugo_Symbol, mut_id)],
  wxs_trimmed[mut_id %in% shared_ids, .(Category = "Shared_WXS", VAF, Hugo_Symbol, mut_id)],
  wgs_trimmed[mut_id %in% shared_ids, .(Category = "Shared_WGS", VAF, Hugo_Symbol, mut_id)]
)

# Remove any rows with missing VAF values
vaf_data <- vaf_data[!is.na(VAF) & !is.infinite(VAF)]

# Summary statistics by category
vaf_summary <- vaf_data[, .(
  Count = .N,
  Mean_VAF = round(mean(VAF), 3),
  Median_VAF = round(median(VAF), 3),
  Q25_VAF = round(quantile(VAF, 0.25), 3),
  Q75_VAF = round(quantile(VAF, 0.75), 3),
  Low_VAF_percent = round(sum(VAF < 0.1) / .N * 100, 1),
  High_VAF_percent = round(sum(VAF > 0.4) / .N * 100, 1)
), by = Category]

