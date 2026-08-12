#######
# PD PROGRESSION ANALYSIS — R SCRIPT
# Cox-MCP Penalized Survival Analysis | PPMI Cohort
#
# Title:   A Transcriptomic Signature of Parkinson's Disease Progression
# Authors: Amevor R, McLain AC, Zhang J, Saha E, Ghosal R
#
# NOTE: PPMI data are access-controlled and are NOT included in this repository.
#       Register and download from https://www.ppmi-info.org to reproduce.
#
# ── USER CONFIG: set these paths to your local data before running 
data_dir    <- "data"                 # folder containing the PPMI CSV files
gmt_dir     <- "data/gmt"             # folder containing the MSigDB .gmt files
output_dir  <- "output"              # folder where plots/tables will be saved

tpm_file    <- file.path(data_dir, "counts_matrix_TPM.csv")
meta_file   <- file.path(data_dir, "PPMI_Curated.csv")
kegg_path     <- file.path(gmt_dir, "c2.cp.kegg.v2022.1.Hs.symbols.gmt")
reactome_path <- file.path(gmt_dir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt.txt")
gobp_path     <- file.path(gmt_dir, "c5.go.bp.v2026.1.Hs.symbols.gmt.txt")


#  SECTION 0: PACKAGES 
library(readr);     library(dplyr);      library(stringr)
library(lubridate); library(survival);   library(survminer)
library(purrr);     library(broom);      library(tibble)
library(glmnet);    library(ncvreg);     library(caret)
library(org.Hs.eg.db); library(AnnotationDbi)
library(ggplot2);   library(tableone);   library(gridExtra)
library(openxlsx);  library(pec)

#  PUBLICATION THEME & COLORS 
PUB_COLORS <- list(
  black      = "#000000",
  dark_gray  = "#404040",
  mid_gray   = "#808080",
  light_gray = "#BDBDBD",
  white      = "#FFFFFF"
)

KM_2_GROUPS <- c("#000000", "#808080")
KM_3_GROUPS <- c("#000000", "#606060", "#BDBDBD")

theme_pub <- function() {
  theme_classic() +
    theme(
      plot.title      = element_text(size = 13, face = "bold",   hjust = 0),
      plot.subtitle   = element_text(size = 11, color = "#404040", hjust = 0),
      axis.title      = element_text(size = 11, face = "bold"),
      axis.text       = element_text(size = 10, color = "black"),
      legend.title    = element_text(size = 10, face = "bold"),
      legend.text     = element_text(size = 10),
      legend.position = "bottom",
      panel.grid      = element_blank(),
      strip.text      = element_text(size = 10, face = "bold")
    )
}

#  OUTPUT FOLDER 
plot_dir <- output_dir
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
cat("Plots will be saved to:", plot_dir, "\n\n")

save_plot <- function(plot_obj, filename, width = 7, height = 6) {
  full_path <- file.path(plot_dir, filename)
  ggsave(paste0(full_path, ".pdf"), plot = plot_obj,
         width = width, height = height, device = "pdf")
  ggsave(paste0(full_path, ".png"), plot = plot_obj,
         width = width, height = height, dpi = 300, device = "png")
  cat("Saved:", filename, "\n")
}

save_km_plot <- function(km_obj, filename, width = 8, height = 7) {
  full_path <- file.path(plot_dir, filename)
  pdf(paste0(full_path, ".pdf"), width = width, height = height)
  print(km_obj)
  dev.off()
  png(paste0(full_path, ".png"), width = width * 300, height = height * 300,
      res = 300)
  print(km_obj)
  dev.off()
  cat("Saved:", filename, "\n")
}


# SECTION 1: DATA LOADING AND PREPROCESSING

tpm <- read_csv(
  tpm_file,
  show_col_types = FALSE
)
gene_ids <- as.character(tpm[[1]])
tpm_mat  <- as.matrix(tpm[, -1, drop = FALSE])
rownames(tpm_mat) <- gene_ids
storage.mode(tpm_mat) <- "numeric"
rm(tpm)
cat("Original TPM:", nrow(tpm_mat), "genes x", ncol(tpm_mat), "samples\n")

prop_required <- 0.10
min_samples   <- ceiling(prop_required * ncol(tpm_mat))
keep_genes    <- rowSums(tpm_mat > 1, na.rm = TRUE) >= min_samples
tpm_filt      <- tpm_mat[keep_genes, , drop = FALSE]
rm(tpm_mat)

cat("Filtered:", nrow(tpm_filt), "genes x", ncol(tpm_filt), "samples\n")
logtpm      <- log2(tpm_filt + 1)
genes_all   <- rownames(logtpm)
samples_all <- colnames(logtpm)
rm(tpm_filt)

sample_info <- tibble(sample = samples_all) %>%
  mutate(
    sample_base = str_remove(sample, "_\\d+$"),
    patno       = as.character(str_extract(sample_base, "^\\d{3,6}")),
    event       = toupper(as.character(str_extract(sample_base, "(?<=_)\\w+")))
  ) %>%
  filter(!is.na(patno))

cat("Total samples:", nrow(sample_info), "\n")
cat("Unique patients:", n_distinct(sample_info$patno), "\n")

### Load Metadata 

meta <- read_csv(
  meta_file,
  show_col_types = FALSE
) %>%
  rename_with(tolower)

meta2 <- meta %>%
  mutate(
    patno             = as.character(patno),
    event_id          = toupper(as.character(event_id)),
    hy                = as.numeric(hy),
    cohort            = as.numeric(cohort),
    age_at_visit      = as.numeric(age_at_visit),
    sex               = sex,
    visit_date_parsed = as.Date(
      parse_date_time(visit_date,
                      orders = c("m/Y","my","mdy","dmy","ymd","Y-m-d"),
                      exact  = FALSE)
    )
  )

cat("Total records:", nrow(meta2), "\n")


# SECTION 2: COHORT DEFINITION

pd_baseline <- meta2 %>%
  filter(cohort == 1, event_id == "BL", hy %in% c(1, 2)) %>%
  dplyr::select(patno, bl_date = visit_date_parsed,
                age_at_visit, sex, bl_hy = hy) %>%
  filter(!is.na(bl_date))

cat("PD patients HY 1-2 at baseline:", nrow(pd_baseline), "\n")
cat("  HY 1:", sum(pd_baseline$bl_hy == 1), "\n")
cat("  HY 2:", sum(pd_baseline$bl_hy == 2), "\n")

visits_all <- meta2 %>%
  dplyr::select(patno, visit_date_parsed, hy) %>%
  filter(!is.na(visit_date_parsed))

conv_all <- visits_all %>%
  inner_join(pd_baseline %>% dplyr::select(patno, bl_date), by = "patno") %>%
  filter(visit_date_parsed > bl_date, !is.na(hy), hy >= 3) %>%
  group_by(patno) %>%
  summarise(conv_date = min(visit_date_parsed), .groups = "drop")

last_vis_all <- visits_all %>%
  group_by(patno) %>%
  summarise(last_date = max(visit_date_parsed, na.rm = TRUE), .groups = "drop")

surv_all <- pd_baseline %>%
  left_join(conv_all,     by = "patno") %>%
  left_join(last_vis_all, by = "patno") %>%
  mutate(
    event       = as.integer(!is.na(conv_date)),
    end_date    = if_else(event == 1L, conv_date, last_date),
    time_months = as.numeric(difftime(end_date, bl_date, units = "days")) / 30.44
  ) %>%
  filter(time_months > 0) %>%
  dplyr::select(patno, bl_hy, event, time_months)

cat("Total with follow-up:", nrow(surv_all), "\n")
cat("Transitioned to HY >= 3:", sum(surv_all$event), "\n")
cat("Censored:", sum(surv_all$event == 0), "\n")

rna_baseline_samples <- sample_info %>%
  filter(event == "BL") %>%
  inner_join(pd_baseline, by = "patno")

final_patnos <- unique(rna_baseline_samples$patno)
cat("Patients with RNA-seq:", length(final_patnos), "\n")

visits <- visits_all %>% filter(patno %in% final_patnos)

conv <- visits %>%
  inner_join(rna_baseline_samples %>% dplyr::select(patno, bl_date),
             by = "patno") %>%
  filter(visit_date_parsed > bl_date, !is.na(hy), hy >= 3) %>%
  group_by(patno) %>%
  summarise(conv_date = min(visit_date_parsed), .groups = "drop")

last_vis <- visits %>%
  group_by(patno) %>%
  summarise(last_date = max(visit_date_parsed, na.rm = TRUE), .groups = "drop")

surv_tbl <- rna_baseline_samples %>%
  dplyr::select(patno, bl_sample = sample, bl_date,
                age_at_visit, sex, bl_hy) %>%
  left_join(conv,     by = "patno") %>%
  left_join(last_vis, by = "patno") %>%
  mutate(
    event       = as.integer(!is.na(conv_date)),
    end_date    = if_else(event == 1L, conv_date, last_date),
    time_months = as.numeric(difftime(end_date, bl_date, units = "days")) / 30.44,
    sex         = factor(sex)
  ) %>%
  filter(time_months > 0)

cat("Final cohort:", nrow(surv_tbl), "patients\n")
cat("  HY 1:", sum(surv_tbl$bl_hy == 1), "\n")
cat("  HY 2:", sum(surv_tbl$bl_hy == 2), "\n")
cat("Transitioned to HY >= 3:", sum(surv_tbl$event), "\n")
cat("Censored:", sum(surv_tbl$event == 0), "\n")
cat("Median follow-up:", round(median(surv_tbl$time_months), 1), "months\n")

meta_bl <- meta %>%
  filter(event_id == "BL") %>%
  mutate(patno = as.character(patno)) %>%
  dplyr::select(patno, updrs3_score, td_pigd, moca, rem)

surv_tbl_ext <- surv_tbl %>%
  left_join(meta_bl, by = "patno")

cat("Clinical variables merged. Complete cases:",
    sum(complete.cases(surv_tbl_ext[,
                                    c("updrs3_score","td_pigd","moca","rem")])), "\n")

# SECTION 3: DESCRIPTIVE STATISTICS & TABLE 1

desc_data <- surv_tbl_ext %>%
  mutate(
    progression = factor(event, levels = c(0,1),
                         labels = c("No Progression", "Progressed to HY \u2265 3")),
    sex_label   = factor(sex, levels = levels(sex),
                         labels = c("Female", "Male")),
    bl_hy_label = factor(bl_hy, levels = c(1, 2),
                         labels = c("HY Stage 1", "HY Stage 2"))
  )

vars     <- c("age_at_visit", "sex_label", "bl_hy_label",
              "time_months", "updrs3_score", "td_pigd", "moca", "rem")
factVars <- c("sex_label", "bl_hy_label")

tab1 <- CreateTableOne(
  vars = vars, strata = "progression",
  data = desc_data, factorVars = factVars, addOverall = TRUE
)

print(tab1, varLabels = TRUE, showAllLevels = TRUE,
      printToggle = TRUE, quote = FALSE, noSpaces = TRUE, test = TRUE)





#  KM PLOTS 
cat("\n--- Figure 1A: KM Overall ---\n")
km_fit  <- survfit(Surv(time_months, event) ~ 1, data = surv_tbl)
km_plot <- ggsurvplot(
  km_fit, data = surv_tbl,
  risk.table = TRUE, conf.int = TRUE,
  conf.int.fill = PUB_COLORS$light_gray,
  palette = PUB_COLORS$black,
  title = "Time to HY \u2265 3 Transition",
  subtitle = "Overall cohort (n = 417, events = 133)",
  xlab = "Time (Months)", ylab = "Progression-Free Probability",
  risk.table.title = "Number at risk", risk.table.y.text = FALSE,
  tables.height = 0.25, ggtheme = theme_pub(), fontsize = 4
)
print(km_plot)
save_km_plot(km_plot, "Figure_1A_KM_Overall", width = 8, height = 7)

cat("--- Figure 1B: KM by HY Stage ---\n")
km_hy <- survfit(Surv(time_months, event) ~ bl_hy, data = surv_tbl)
km_hy_plot <- ggsurvplot(
  km_hy, data = surv_tbl,
  risk.table = TRUE, conf.int = TRUE, pval = TRUE, pval.method = TRUE,
  palette = KM_2_GROUPS,
  legend.labs = c("HY Stage 1", "HY Stage 2"), legend.title = "Baseline HY Stage",
  title = "Time to HY \u2265 3 by Baseline HY Stage",
  xlab = "Time (Months)", ylab = "Progression-Free Probability",
  risk.table.y.text = FALSE, tables.height = 0.25,
  ggtheme = theme_pub(), fontsize = 4
)
print(km_hy_plot)
save_km_plot(km_hy_plot, "Figure_1B_KM_by_HY_Stage", width = 8, height = 7)

cat("--- Figure 1C: KM by Sex ---\n")
km_sex <- survfit(Surv(time_months, event) ~ sex, data = surv_tbl)
km_sex_plot <- ggsurvplot(
  km_sex, data = surv_tbl,
  risk.table = TRUE, conf.int = TRUE, pval = TRUE, pval.method = TRUE,
  palette = KM_2_GROUPS,
  legend.labs = c("Female", "Male"), legend.title = "Sex",
  title = "Time to HY \u2265 3 by Sex",
  xlab = "Time (Months)", ylab = "Progression-Free Probability",
  risk.table.y.text = FALSE, tables.height = 0.25,
  ggtheme = theme_pub(), fontsize = 4
)
print(km_sex_plot)
save_km_plot(km_sex_plot, "Figure_1C_KM_by_Sex", width = 8, height = 7)


# SECTION 4: GENE SCREENING


genes_to_test <- rownames(logtpm)
cat("Total genes:", length(genes_to_test), "\n")

get_baseline_expr <- function(gene) {
  as.numeric(logtpm[gene, surv_tbl$bl_sample])
}

cat("\n--- Univariate Cox ---\n")
fit_cox_univar <- function(gene) {
  expr_vec <- get_baseline_expr(gene)
  df   <- surv_tbl %>% mutate(gexpr = as.numeric(scale(expr_vec)))
  fit  <- coxph(Surv(time_months, event) ~ gexpr, data = df)
  raw  <- tidy(fit, exponentiate = FALSE)
  expd <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tibble(
    gene    = gene,
    coef    = raw$estimate[raw$term == "gexpr"],
    se      = raw$std.error[raw$term == "gexpr"],
    p_value = raw$p.value[raw$term == "gexpr"],
    HR      = expd$estimate[expd$term == "gexpr"],
    HR_low  = expd$conf.low[expd$term == "gexpr"],
    HR_high = expd$conf.high[expd$term == "gexpr"]
  )
}

cox_res_all_univar <- map_dfr(genes_to_test, fit_cox_univar) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

cat("Univariate: p < 0.05:", sum(cox_res_all_univar$p_value < 0.05),
    "| FDR < 0.05:", sum(cox_res_all_univar$p_adj < 0.05), "\n")

cat("\n--- Age/Sex Adjusted Cox ---\n")
fit_cox_adj <- function(gene) {
  expr_vec <- get_baseline_expr(gene)
  df   <- surv_tbl %>% mutate(gexpr = as.numeric(scale(expr_vec)))
  fit  <- coxph(Surv(time_months, event) ~ gexpr + age_at_visit + sex, data = df)
  raw  <- tidy(fit, exponentiate = FALSE)
  expd <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tibble(
    gene    = gene,
    coef    = raw$estimate[raw$term == "gexpr"],
    se      = raw$std.error[raw$term == "gexpr"],
    p_value = raw$p.value[raw$term == "gexpr"],
    HR      = expd$estimate[expd$term == "gexpr"],
    HR_low  = expd$conf.low[expd$term == "gexpr"],
    HR_high = expd$conf.high[expd$term == "gexpr"]
  )
}

cox_res_all_adj <- map_dfr(genes_to_test, fit_cox_adj) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

cat("Adjusted: p < 0.05:", sum(cox_res_all_adj$p_value < 0.05),
    "| FDR < 0.05:", sum(cox_res_all_adj$p_adj < 0.05), "\n")

cox_res_all_adj <- cox_res_all_adj %>%
  mutate(gene_clean = sub("\\..*", "", gene))

gene_info <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(cox_res_all_adj$gene_clean),
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENSEMBL"
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

cox_res_all_adj <- cox_res_all_adj %>%
  left_join(gene_info, by = c("gene_clean" = "ENSEMBL"))

# ALL 62 SIGNIFICANT GENES 
sig_62 <- cox_res_all_adj %>%
  filter(p_value < 0.05) %>%
  dplyr::select(gene, SYMBOL, GENENAME, coef, HR, HR_low, HR_high,
                p_value, p_adj) %>%
  mutate(
    coef      = round(coef, 3),
    HR        = round(HR, 3),
    HR_low    = round(HR_low, 3),
    HR_high   = round(HR_high, 3),
    p_value   = signif(p_value, 3),
    p_adj     = signif(p_adj, 3),
    direction = ifelse(coef > 0, "Risk", "Protective")
  ) %>%
  arrange(p_value)


cat("★ ALL", nrow(sig_62), "SIGNIFICANT GENES (Adjusted Cox, p < 0.05) ★\n")

print(as.data.frame(sig_62), row.names = FALSE)


# SECTION 5: PENALIZED GENE SELECTION


sig_genes <- sig_62$gene
X_sig     <- t(logtpm[sig_genes, surv_tbl$bl_sample])
Y_surv    <- Surv(surv_tbl$time_months, surv_tbl$event)
set.seed(2025)
fold      <- createFolds(cut(surv_tbl$time_months, breaks = 5), k = 5, list = FALSE)

# Z-SCORE GENE EXPRESSION 
X_sig <- scale(X_sig)   # mean = 0, sd = 1 per gene
# 

X_all          <- cbind(X_sig, age = surv_tbl$age_at_visit,
                        sex = as.numeric(surv_tbl$sex))
penalty_factor <- c(rep(1, ncol(X_sig)), 0, 0)

# LASSO
lasso_fit   <- cv.glmnet(X_all, Y_surv, family = "cox", alpha = 1,
                         penalty.factor = penalty_factor)
lasso_coefs <- coef(lasso_fit, s = "lambda.min")
lasso_genes <- rownames(lasso_coefs)[which(lasso_coefs != 0 &
                                             !rownames(lasso_coefs) %in% c("age", "sex"))]

# MCP
mcp_fit   <- cv.ncvsurv(X_all, Y_surv, penalty = "MCP",
                        penalty.factor = penalty_factor)
mcp_coefs <- coef(mcp_fit, lambda = mcp_fit$lambda.min)
mcp_genes <- names(mcp_coefs)[which(mcp_coefs != 0 &
                                      !names(mcp_coefs) %in% c("age", "sex"))]

# SCAD
scad_fit   <- cv.ncvsurv(X_all, Y_surv, penalty = "SCAD",
                         penalty.factor = penalty_factor)
scad_coefs <- coef(scad_fit, lambda = scad_fit$lambda.min)
scad_genes <- names(scad_coefs)[which(scad_coefs != 0 &
                                        !names(scad_coefs) %in% c("age", "sex"))]

cat("LASSO:", length(lasso_genes), "genes\n")
cat("MCP:  ", length(mcp_genes),   "genes\n")
cat("SCAD: ", length(scad_genes),  "genes\n")
cat("\nMCP genes:\n"); print(mcp_genes)

core_genes <- Reduce(intersect, list(lasso_genes, mcp_genes, scad_genes))
cat("\nCore (all 3 methods):", length(core_genes), "genes\n")
print(core_genes)

# FIGURE 2: Venn Diagram 

if (!requireNamespace("ggVennDiagram", quietly = TRUE))
  install.packages("ggVennDiagram")
library(ggVennDiagram)

gene_sets <- list(LASSO = lasso_genes, MCP = mcp_genes, SCAD = scad_genes)

venn_plot <- ggVennDiagram(
  gene_sets, label_alpha = 0,
  category.names = c(
    paste0("LASSO\n(", length(lasso_genes), " genes)"),
    paste0("MCP\n(",   length(mcp_genes),   " genes)"),
    paste0("SCAD\n(",  length(scad_genes),  " genes)")
  )
) +
  scale_fill_gradient(low = "white", high = PUB_COLORS$light_gray) +
  scale_color_manual(values = rep(PUB_COLORS$dark_gray, 3)) +
  labs(
    title    = "Gene Overlap Across Penalized Methods",
    subtitle = paste0(length(core_genes),
                      " genes consistently selected by all three methods")
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.line  = element_blank(), axis.text  = element_blank(),
    axis.ticks = element_blank(), axis.title = element_blank()
  )

print(venn_plot)
save_plot(venn_plot, "Figure_2_Venn_Diagram", width = 7, height = 6)


# SECTION 6: CROSS-VALIDATION C-INDICES


cal_c <- function(marker, Stime, status) {
  utimes <- sort(unique(Stime[status == 1]))
  num <- denom <- 0
  for (ti in utimes) {
    inx_i <- which(Stime == ti & status == 1)
    inx_j <- which(Stime > ti)
    for (i in seq_along(inx_i)) {
      num <- num + sum(marker[inx_j] > marker[inx_i[i]]) +
        0.5 * sum(marker[inx_j] == marker[inx_i[i]])
    }
    denom <- denom + length(inx_i) * length(inx_j)
  }
  1 - num / denom
}

# Model 1: Age + Sex — full-data concordance (CV unstable with 2 predictors)
fit_m1  <- coxph(Surv(time_months, event) ~ age_at_visit + sex, data = surv_tbl)
c1_full <- summary(fit_m1)$concordance[1]
cat("Model 1 (Age+Sex): Full-data C =", round(c1_full, 3), "\n")

# Complete cases for full models 
complete_idx    <- complete.cases(surv_tbl_ext[, c("age_at_visit","sex",
                                                   "updrs3_score","td_pigd","moca","rem")])
surv_complete   <- surv_tbl_ext[complete_idx, ]
fold_complete   <- fold[complete_idx]
X_sig_complete  <- X_sig[complete_idx, ]   # inherits z-scoring from Section 5
y_surv_complete <- Surv(surv_complete$time_months, surv_complete$event)

cat("Complete cases:", sum(complete_idx), "| Events:", sum(surv_complete$event), "\n")

# Model 2: Clinical only
cv2 <- numeric(5)
for (k in 1:5) {
  tr <- fold_complete != k; te <- fold_complete == k
  dat_tr <- surv_complete[tr, c("age_at_visit","sex","updrs3_score","td_pigd","moca","rem")]
  dat_te <- surv_complete[te, c("age_at_visit","sex","updrs3_score","td_pigd","moca","rem")]
  fit  <- coxph(y_surv_complete[tr] ~ ., data = dat_tr)
  pred <- predict(fit, newdata = dat_te, type = "lp")
  cv2[k] <- cal_c(pred, surv_complete$time_months[te], surv_complete$event[te])
}
cat("Model 2 (Clinical): CV C =", round(mean(cv2), 3), "\n")

# Model 3: Cox-MCP full (primary)
data_mcp_full <- data.frame(
  X_sig_complete[, mcp_genes, drop = FALSE],
  age     = surv_complete$age_at_visit,
  sex     = surv_complete$sex,
  updrs3  = surv_complete$updrs3_score,
  td_pigd = surv_complete$td_pigd,
  moca    = surv_complete$moca,
  rem     = surv_complete$rem
)

# x = TRUE required for pec (Brier score)
cox_mcp_full    <- coxph(y_surv_complete ~ ., data = data_mcp_full, x = TRUE)
pred_mcp_full   <- predict(cox_mcp_full, type = "lp")
c_full_mcp_full <- cal_c(pred_mcp_full, surv_complete$time_months,
                         surv_complete$event)

cv_mcp_full <- numeric(5)
for (k in 1:5) {
  tr <- fold_complete != k; te <- fold_complete == k
  dat_tr <- data.frame(X_sig_complete[tr, mcp_genes, drop=FALSE],
                       age=surv_complete$age_at_visit[tr], sex=surv_complete$sex[tr],
                       updrs3=surv_complete$updrs3_score[tr], td_pigd=surv_complete$td_pigd[tr],
                       moca=surv_complete$moca[tr], rem=surv_complete$rem[tr])
  dat_te <- data.frame(X_sig_complete[te, mcp_genes, drop=FALSE],
                       age=surv_complete$age_at_visit[te], sex=surv_complete$sex[te],
                       updrs3=surv_complete$updrs3_score[te], td_pigd=surv_complete$td_pigd[te],
                       moca=surv_complete$moca[te], rem=surv_complete$rem[te])
  fit  <- coxph(y_surv_complete[tr] ~ ., data = dat_tr)
  pred <- predict(fit, newdata = dat_te, type = "lp")
  cv_mcp_full[k] <- cal_c(pred, surv_complete$time_months[te],
                          surv_complete$event[te])
}
cat("Model 3 (Cox-MCP Full): CV C =", round(mean(cv_mcp_full), 3),
    "| Full C =", round(c_full_mcp_full, 3), "\n")

cat("\n--- Summary: Cox-MCP Full Model ---\n")
print(summary(cox_mcp_full))

cat("\n--- Proportional Hazards Test ---\n")
print(cox.zph(cox_mcp_full))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6b: LEAKAGE-FREE NESTED CROSS-VALIDATION
# MCP selection is repeated INSIDE each training fold, so held-out folds never
# contribute to gene selection.
# ══════════════════════════════════════════════════════════════════════════════


X_all_c <- X_all[complete_idx, ]

set.seed(2025)
cv_nested <- numeric(5); n_genes_perfold <- numeric(5)
for (k in 1:5) {
  tr <- fold_complete != k
  te <- fold_complete == k
  fit_k  <- cv.ncvsurv(X_all_c[tr, ], y_surv_complete[tr],
                       penalty = "MCP", penalty.factor = penalty_factor)
  beta_k <- as.numeric(coef(fit_k, lambda = fit_k$lambda.min))
  n_genes_perfold[k] <- sum(beta_k[seq_len(ncol(X_sig))] != 0)
  lp_te  <- as.numeric(X_all_c[te, ] %*% beta_k)
  cv_nested[k] <- cal_c(lp_te, surv_complete$time_months[te], surv_complete$event[te])
}
cat("Nested-CV MCP (62 genes + age/sex): C =", round(mean(cv_nested), 3), "\n")
cat("  Genes selected per fold:", paste(n_genes_perfold, collapse = ", "), "\n")

cat("\n--- Nested-CV comparison ---\n")
cat("Age + Sex + Clinical    :", round(mean(cv2), 3), "\n")
cat("Nested MCP + Age/Sex    :", round(mean(cv_nested), 3), "\n")
cat("(Fixed-signature MCP Model 3 =", round(mean(cv_mcp_full), 3), ")\n")


# SECTION 7: FOREST PLOT — MCP GENE EFFECTS


mcp_results <- tidy(cox_mcp_full, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(!term %in% c("age", "sex", "updrs3", "td_pigd", "moca", "rem")) %>%
  mutate(
    gene_clean = sub("\\..*", "", term),
    p_val_star = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  left_join(gene_info, by = c("gene_clean" = "ENSEMBL")) %>%
  mutate(
    label     = ifelse(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, term),
    direction = ifelse(estimate > 1, "Risk (HR > 1)", "Protective (HR < 1)")
  ) %>%
  filter(p.value < 0.05) %>%
  arrange(estimate)

forest_plot <- ggplot(mcp_results,
                      aes(x = estimate, y = reorder(label, estimate),
                          xmin = conf.low, xmax = conf.high, color = direction)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.25, linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = PUB_COLORS$dark_gray, linewidth = 0.6) +
  geom_text(aes(label = p_val_star, x = conf.high + 0.05),
            hjust = 0, size = 3.5, color = "black") +
  scale_color_manual(values = c(
    "Risk (HR > 1)"       = PUB_COLORS$dark_gray,
    "Protective (HR < 1)" = PUB_COLORS$mid_gray)) +
  scale_x_log10() +
  labs(
    title    = "Forest Plot: Cox-MCP Gene Effects",
    subtitle = "Primary model (n = 411, 130 events). Genes reaching p < 0.05 shown.",
    x        = "Hazard Ratio (95% CI) \u2014 log scale",
    y        = "Gene", color = "Direction"
  ) +
  theme_pub() +
  theme(legend.position = "right", axis.text.y = element_text(size = 9))

print(forest_plot)
save_plot(forest_plot, "Figure_3_Forest_Plot_MCP_Genes", width = 9, height = 7)

# SECTION 8: RISK GROUP STRATIFICATION + KM PLOT


pred_risk  <- predict(cox_mcp_full, type = "lp")
tertiles   <- quantile(pred_risk, probs = c(1/3, 2/3))
risk_group <- cut(pred_risk,
                  breaks = c(-Inf, tertiles[1], tertiles[2], Inf),
                  labels = c("Low Risk", "Mid Risk", "High Risk"))

surv_complete_risk <- surv_complete %>% mutate(risk_group = risk_group)

cat("Risk group distribution:\n")
print(table(risk_group))

cat("\nEvent rates by risk group:\n")
surv_complete_risk %>%
  group_by(risk_group) %>%
  summarise(n = n(), events = sum(event),
            pct = round(100 * mean(event), 1), .groups = "drop") %>%
  print()

lr_test <- survdiff(Surv(time_months, event) ~ risk_group,
                    data = surv_complete_risk)
cat("\nLog-rank p =", round(1 - pchisq(lr_test$chisq, df = 2), 5), "\n")

km_risk <- survfit(Surv(time_months, event) ~ risk_group,
                   data = surv_complete_risk)
km_risk_plot <- ggsurvplot(
  km_risk, data = surv_complete_risk,
  risk.table = TRUE, conf.int = TRUE, pval = TRUE, pval.method = TRUE,
  palette = KM_3_GROUPS, legend.title = "Risk Group",
  legend.labs = c("Low Risk", "Mid Risk", "High Risk"),
  title    = "Kaplan\u2013Meier Curves by Predicted Risk Group",
  subtitle = "Final Cox-MCP Model: 25 genes + age + sex + clinical variables",
  xlab = "Time (Months)", ylab = "Progression-Free Probability",
  risk.table.title = "Number at risk", risk.table.y.text = FALSE,
  tables.height = 0.28, ggtheme = theme_pub(), fontsize = 4,
  font.main = c(13, "bold"), font.subtitle = c(10, "plain", "#404040")
)
print(km_risk_plot)
save_km_plot(km_risk_plot, "Figure_4_KM_Risk_Groups", width = 9, height = 8)


# SECTION 9: BRIER SCORE (CALIBRATION)


data_for_pec             <- data_mcp_full
data_for_pec$time_months <- surv_complete$time_months
data_for_pec$event       <- surv_complete$event

brier_result <- pec(
  object  = list("Cox-MCP" = cox_mcp_full),
  formula = Surv(time_months, event) ~ 1,
  data    = data_for_pec,
  times   = seq(12, 120, by = 12)
)

print(brier_result)
cat("\n--- Integrated Brier Score ---\n")
crps(brier_result)

# Save Brier score plot — PDF
pdf(file.path(plot_dir, "Figure_5_Brier_Score.pdf"), width = 8, height = 6)
plot(brier_result,
     main = "Brier Score Over Time \u2014 Cox-MCP Model",
     xlab = "Time (Months)", ylab = "Brier Score",
     col  = c("black", "gray50"), lwd  = 2)
legend("topleft",
       legend = c("Reference (Kaplan-Meier)", "Cox-MCP Model"),
       col    = c("black", "gray50"),   # ← swap these
       lwd    = 2)
dev.off()

# Save Brier score plot — PNG
png(file.path(plot_dir, "Figure_5_Brier_Score.png"),
    width = 8*300, height = 6*300, res = 300)
plot(brier_result,
     main = "Brier Score Over Time \u2014 Cox-MCP Model",
     xlab = "Time (Months)", ylab = "Brier Score",
     col  = c("black", "gray50"), lwd  = 2)
legend("topleft",
       legend = c("Reference (Kaplan-Meier)", "Cox-MCP Model"),
       col    = c("black", "gray50"),   # ← swap these
       lwd    = 2)

dev.off()

cat("Saved: Figure_5_Brier_Score\n")


####################
# GENE SET ENRICHMENT ANALYSIS (GSEA) — PD Progression
# Uses: KEGG, Reactome, GO Biological Process (MSigDB)
# Input: Univariate Cox HR for all 17,168 genes
################

library(fgsea)
library(dplyr)
library(tibble)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

#  Output folder 
plot_dir <- output_dir
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# GMT file paths 
# kegg_path defined in USER CONFIG above
# reactome_path defined in USER CONFIG above
# gobp_path defined in USER CONFIG above

##############################
# STEP 1: PREPARE RANKED GENE LIST
# Rank all 17,168 genes by coefficient (log HR) from univariate Cox
# Positive = risk (faster progression), Negative = protective
# ###################

# using cox_res_all_univar from main analysis

# Map ENSG IDs to gene symbols
gene_clean <- sub("\\..*", "", cox_res_all_univar$gene)

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(gene_clean),
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# Build ranked vector: name = gene symbol, value = coefficient (log HR)
ranked_df <- cox_res_all_univar %>%
  mutate(gene_clean = sub("\\..*", "", gene)) %>%
  left_join(symbol_map, by = c("gene_clean" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  # Handle duplicate symbols: keep highest absolute coefficient
  group_by(SYMBOL) %>%
  slice_max(abs(coef), n = 1, with_ties = FALSE) %>%
  ungroup()

# Create named vector ranked by coefficient
ranked_gene_list        <- ranked_df$coef
names(ranked_gene_list) <- ranked_df$SYMBOL
ranked_gene_list        <- sort(ranked_gene_list, decreasing = TRUE)

cat("Total genes with symbols:", length(ranked_gene_list), "\n")
cat("Top 5 risk genes (highest coef):\n")
print(head(ranked_gene_list, 5))
cat("Top 5 protective genes (lowest coef):\n")
print(tail(ranked_gene_list, 5))

# STEP 2: LOAD PATHWAYS


pathways_kegg     <- gmtPathways(kegg_path)
pathways_reactome <- gmtPathways(reactome_path)
pathways_gobp     <- gmtPathways(gobp_path)

cat("KEGG pathways:     ", length(pathways_kegg),     "\n")
cat("Reactome pathways: ", length(pathways_reactome), "\n")
cat("GO BP pathways:    ", length(pathways_gobp),     "\n")


# STEP 3: RUN FGSEA — ALL THREE DATABASES

set.seed(2025)

run_fgsea <- function(pathways, ranked_list, db_name) {
  cat("\n--- Running fgsea:", db_name, "---\n")
  res <- fgsea(
    pathways  = pathways,
    stats     = ranked_list,
    minSize   = 15,
    maxSize   = 500
  ) %>%
    arrange(padj) %>%
    mutate(database = db_name)
  cat("Total pathways tested:", nrow(res), "\n")
  cat("Significant (padj < 0.05):", sum(res$padj < 0.05, na.rm=TRUE), "\n")
  cat("Significant (padj < 0.25):", sum(res$padj < 0.25, na.rm=TRUE), "\n")
  res
}

res_kegg     <- run_fgsea(pathways_kegg,     ranked_gene_list, "KEGG")
res_reactome <- run_fgsea(pathways_reactome, ranked_gene_list, "Reactome")
res_gobp     <- run_fgsea(pathways_gobp,     ranked_gene_list, "GO BP")


# STEP 4: RESULTS SUMMARY

print_top10 <- function(res, db_name) {
  cat("\n=== TOP 10 SIGNIFICANT PATHWAYS:", db_name, "===\n")
  top10 <- res %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(10) %>%
    dplyr::select(pathway, NES, pval, padj, size) %>%
    mutate(
      NES       = round(NES,  3),
      pval      = signif(pval, 3),
      padj      = signif(padj, 3),
      direction = ifelse(NES > 0, "Risk", "Protective")
    )
  print(as.data.frame(top10))
  top10
}

# the summaries
top10_kegg     <- print_top10(res_kegg,     "KEGG")
top10_reactome <- print_top10(res_reactome, "Reactome")
top10_gobp     <- print_top10(res_gobp,     "GO BP")

# STEP 5: PLOTS

cat("\n========== STEP 5: PLOTS ==========\n")

plot_top_pathways <- function(res, db_name, n = 15) {
  top_pos <- res %>% filter(NES > 0, padj < 0.05) %>% arrange(padj) %>% head(ceiling(n/2))
  top_neg <- res %>% filter(NES < 0, padj < 0.05) %>% arrange(padj) %>% head(floor(n/2))
  top_plot <- bind_rows(top_pos, top_neg) %>%
    mutate(
      pathway_clean = str_remove(pathway, paste0("^", db_name, "_")) %>%
        str_replace_all("_", " ") %>%
        str_to_title() %>%
        str_trunc(55),
      direction = ifelse(NES > 0, "Risk (NES > 0)", "Protective (NES < 0)")
    ) %>%
    arrange(NES)
  
  if (nrow(top_plot) == 0) {
    cat("No significant pathways for", db_name, "- skipping plot\n")
    return(NULL)
  }
  
  p <- ggplot(top_plot,
              aes(x = NES, y = reorder(pathway_clean, NES), fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_fill_manual(
      values = c("Risk (NES > 0)"       = "#404040",
                 "Protective (NES < 0)" = "#BDBDBD")
    ) +
    labs(
      title    = paste0("GSEA \u2014 ", db_name, " Pathways"),
      subtitle = paste0("Pathways associated with PD HY \u2265 3 progression (padj < 0.05)"),
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL,
      fill     = "Direction"
    ) +
    theme_classic() +
    theme(
      plot.title      = element_text(size = 12, face = "bold"),
      plot.subtitle   = element_text(size = 10, color = "#404040"),
      axis.text.y     = element_text(size = 9),
      axis.text.x     = element_text(size = 9),
      legend.position = "bottom"
    )
  print(p)
  
  fname <- paste0("Figure_GSEA_", db_name)
  ggsave(file.path(plot_dir, paste0(fname, ".pdf")), p,
         width = 10, height = 7, device = "pdf")
  ggsave(file.path(plot_dir, paste0(fname, ".png")), p,
         width = 10, height = 7, dpi = 300, device = "png")
  cat("Saved:", fname, "\n")
  p
}

# the plots
library(stringr)
p_kegg     <- plot_top_pathways(res_kegg,     "KEGG")
p_reactome <- plot_top_pathways(res_reactome, "Reactome")
p_gobp     <- plot_top_pathways(res_gobp,     "GOBP")



# STEP 6: SAVE RESULTS TO EXCEL
library(openxlsx)

# Convert leadingEdge list to string for Excel
clean_for_excel <- function(res) {
  res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse=", ")) %>%
    dplyr::select(pathway, NES, pval, padj, size, leadingEdge, database)
}

wb <- createWorkbook()
addWorksheet(wb, "KEGG")
addWorksheet(wb, "Reactome")
addWorksheet(wb, "GO_BP")
addWorksheet(wb, "Top10_Summary")

writeData(wb, "KEGG",     clean_for_excel(res_kegg))
writeData(wb, "Reactome", clean_for_excel(res_reactome))
writeData(wb, "GO_BP",    clean_for_excel(res_gobp))

# Summary sheet — top 10 from each
summary_all <- bind_rows(
  top10_kegg %>% mutate(database = "KEGG"),
  top10_reactome %>% mutate(database = "Reactome"),
  top10_gobp %>% mutate(database = "GO BP")
) %>% dplyr::select(database, pathway, NES, pval, padj, size, direction)

writeData(wb, "Top10_Summary", summary_all)
saveWorkbook(wb, file.path(plot_dir, "GSEA_Results_All_Databases.xlsx"),
             overwrite = TRUE)
cat("\nSaved: GSEA_Results_All_Databases.xlsx\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10: DIFFERENTIAL EXPRESSION — PD vs HEALTHY CONTROLS (limma-trend)
# Complementary analysis testing whether a disease-presence signature (PD vs HC)
# predicts progression. Produces Supplementary Table S4.
# Uses limma-trend on log2-TPM (eBayes trend = TRUE, robust = TRUE).
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== SECTION 10: DIFFERENTIAL EXPRESSION (PD vs Control) ==========\n")

library(limma)

# PD patients at baseline (cohort 1)
pd_baseline_de <- meta2 %>%
  filter(cohort == 1, event_id == "BL") %>%
  dplyr::select(patno, age_at_visit, sex) %>%
  filter(!is.na(age_at_visit), !is.na(sex))
pd_samples <- sample_info %>%
  filter(event == "BL") %>%
  inner_join(pd_baseline_de, by = "patno") %>%
  mutate(group = "PD") %>%
  dplyr::select(sample, patno, age = age_at_visit, sex, group)

# Healthy controls at baseline (cohort 2)
controls_baseline <- meta2 %>%
  filter(cohort == 2, event_id == "BL") %>%
  dplyr::select(patno, age_at_visit, sex) %>%
  filter(!is.na(age_at_visit), !is.na(sex))
control_samples <- sample_info %>%
  filter(event == "BL") %>%
  inner_join(controls_baseline, by = "patno") %>%
  mutate(group = "Control") %>%
  dplyr::select(sample, patno, age = age_at_visit, sex, group)

# Combine (Control = reference; positive logFC = higher in PD)
limma_samples <- bind_rows(pd_samples, control_samples) %>%
  mutate(group = factor(group, levels = c("Control", "PD")), sex = factor(sex))
cat("PD:", sum(limma_samples$group == "PD"),
    "| Control:", sum(limma_samples$group == "Control"),
    "| Total:", nrow(limma_samples), "\n")

# Expression matrix (log2-TPM), design, fit
expr_limma <- logtpm[, limma_samples$sample, drop = FALSE]
design <- model.matrix(~ age + sex + group, data = limma_samples)
fit    <- lmFit(expr_limma, design)
fit    <- eBayes(fit, trend = TRUE, robust = TRUE)

de_results <- topTable(fit, coef = "groupPD", number = Inf, adjust.method = "BH")
de_results$gene <- rownames(de_results)
cat("Genes tested:", nrow(de_results),
    "| Significant (adj.P < 0.05):", sum(de_results$adj.P.Val < 0.05), "\n")

# Top 50 for Supplementary Table S4
de_top50 <- head(de_results[order(de_results$P.Value), ], 50)
cat("\n--- Top 50 DE genes (Table S4) ---\n")
print(de_top50[, c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")],
      row.names = FALSE)


