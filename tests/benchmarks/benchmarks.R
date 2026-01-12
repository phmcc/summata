################################################################################
##                    summata Performance Benchmarks
## 
## Comparing summata functions against popular alternatives:
##   - Descriptive tables: desctable() vs tableby, tbl_summary, CreateTableOne,
##                         summary_factorlist
##   - Regression output: fit() vs tbl_regression, broom::tidy, glmmulti+fit2df
##   - Univariable screening: uniscreen() vs tbl_uvregression, modelsum, 
##                            glmuni+fit2df
##   - Complete workflows: fullfit() vs finalfit(), gtsummary multi-step
##   - Forest plots: glmforest/coxforest() vs ggforest, manual ggplot2
##   - Mixed-effects models: fit() with lmer/glmer/coxme vs alternatives
##
## Author: Paul H. McClelland
## Date: 2025-12-20
################################################################################

## =============================================================================
## SETUP
## =============================================================================

required_packages <- c(
    "summata", "finalfit", "arsenal", "gtsummary", "tableone",
    "broom", "broom.mixed", "survminer", "data.table", "microbenchmark",
    "ggplot2", "survival", "lme4", "coxme"
)

install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing package: ", pkg)
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}
invisible(lapply(required_packages, install_if_missing))

library(summata)
library(finalfit)
library(arsenal)
library(gtsummary)
library(tableone)
library(broom)
library(broom.mixed)
library(survminer)
library(data.table)
library(microbenchmark)
library(ggplot2)
library(survival)
library(lme4)
library(coxme)
library(ragg)

set.seed(71)
setDTthreads(0)

## =============================================================================
## GENERATE BENCHMARK DATASETS
## =============================================================================

generate_benchmark_data <- function(n, seed = 71) {
    set.seed(seed)
    
    data.table(
        id = 1:n,
        site = factor(sample(paste0("Site_", 1:10), n, replace = TRUE)),
        age = round(rnorm(n, mean = 55, sd = 12)),
        sex = factor(sample(c("Male", "Female"), n, replace = TRUE)),
        race = factor(sample(c("White", "Black", "Asian", "Other"), n, 
                             replace = TRUE, prob = c(0.6, 0.2, 0.15, 0.05)),
                      levels = c("White", "Black", "Asian", "Other")),
        bmi = round(rnorm(n, mean = 27, sd = 5), 1),
        smoking = factor(sample(c("Never", "Former", "Current"), n, 
                                replace = TRUE, prob = c(0.4, 0.35, 0.25)),
                         levels = c("Never", "Former", "Current")),
        hypertension = factor(sample(0:1, n, replace = TRUE, prob = c(0.7, 0.3))),
        diabetes = factor(sample(0:1, n, replace = TRUE, prob = c(0.85, 0.15))),
        treatment = factor(sample(c("Placebo", "Drug A", "Drug B"), n, replace = TRUE),
                           levels = c("Placebo", "Drug A", "Drug B")),
        stage = factor(sample(c("I", "II", "III", "IV"), n, 
                              replace = TRUE, prob = c(0.2, 0.3, 0.35, 0.15)),
                       levels = c("I", "II", "III", "IV")),
        grade = factor(sample(1:3, n, replace = TRUE), 
                       levels = 1:3, labels = c("Low", "Moderate", "High")),
        biomarker1 = round(rlnorm(n, meanlog = 2, sdlog = 0.8), 2),
        biomarker2 = round(rnorm(n, mean = 100, sd = 25), 1),
        creatinine = round(rnorm(n, mean = 1.1, sd = 0.3), 2),
        hemoglobin = round(rnorm(n, mean = 13.5, sd = 1.8), 1),
        los_days = round(rexp(n, rate = 0.1) + 1, 1),
        event_count = rpois(n, lambda = 2),
        response = rbinom(n, 1, prob = 0.3),
        os_time = round(rexp(n, rate = 0.02) + 1, 1),
        os_status = rbinom(n, 1, prob = 0.4)
    )
}

cat("Generating benchmark datasets...\n")
data_50k <- generate_benchmark_data(50000)

data_500 <- data_50k[1:500]
data_1k <- data_50k[1:1000]
data_5k <- data_50k[1:5000]
data_10k <- data_50k[1:10000]

df_500 <- as.data.frame(data_500)
df_1k <- as.data.frame(data_1k)
df_5k <- as.data.frame(data_5k)
df_10k <- as.data.frame(data_10k)
df_50k <- as.data.frame(data_50k)

vars_continuous <- c("age", "bmi", "biomarker1", "biomarker2", "creatinine", "hemoglobin")
vars_categorical <- c("sex", "race", "smoking", "hypertension", "diabetes", "treatment", "stage", "grade")
vars_all <- c(vars_continuous, vars_categorical)

cat("Datasets ready: 500, 1k, 5k, 10k, 50k observations\n\n")

## =============================================================================
## COMMON PLOTTING THEME AND COLOR PALETTE
## =============================================================================

## Define consistent colors for all packages
pkg_colors <- c(
    "summata" = "#FC8D62",
    "summata_minimal" = "#FCB362",
    "finalfit" = "#5597DF",
    "gtsummary" = "#66C2A5",
    "arsenal" = "#A6D854",
    "tableone" = "#8DA0CB",
    "broom" = "#FFD92F",
    "broom_loop" = "#E5C494",
    "manual" = "#E78AC3",
    "survminer" = "#FB9A99"
)

## Helper function to reorder factor levels by median time (fastest to slowest)
order_by_speed <- function(dt, time_col = "time_ms", group_col = "expr") {
    ## Calculate median time for each package
    medians <- dt[, .(med_time = median(get(time_col))), by = group_col]
    setorder(medians, med_time)
    
    ## Reorder factor levels
    dt[, (group_col) := factor(get(group_col), levels = medians[[group_col]])]
    return(dt)
}

## Font setup for consistent rendering across devices
if (requireNamespace("ragg", quietly = TRUE)) {
    ## Use ragg for high-quality PNG output with proper font rendering
    knitr::opts_chunk$set(dev = "ragg_png") # For knitr/rmarkdown
    
    ## For ggsave, we'll use the device argument
    png_device <- ragg::agg_png
    cat("Using ragg for PNG rendering (recommended)\n")
} else if (capabilities("cairo")) {
    png_device <- function(...) grDevices::png(..., type = "cairo")
    cat("Using cairo for PNG rendering\n")
} else {
    png_device <- grDevices::png
    cat("Using default PNG device (fonts may differ from PDF)\n")
    message("For consistent fonts, install the 'ragg' package: install.packages('ragg')")
}

## Custom ggsave wrapper that uses appropriate device for each format
save_plot <- function(filename, plot, width = 10, height = 5, dpi = 150, ...) {
    ext <- tools::file_ext(filename)
    
    if (ext == "png" && exists("png_device")) {
        ggsave(filename, plot, width = width, height = height, dpi = dpi,
               device = png_device, ...)
    } else {
        ggsave(filename, plot, width = width, height = height, dpi = dpi, ...)
    }
    cat("Saved:", filename, "\n")
}

theme_benchmark <- function() {
    theme_minimal(base_family = "Helvetica") +
        theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10, color = "gray40"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "gray95", color = NA),
            strip.text = element_text(face = "bold")
        )
}

## =============================================================================
## BENCHMARK 1: DESCRIPTIVE STATISTICS TABLES
## =============================================================================

cat(strrep("=", 70), "\n")
cat("BENCHMARK 1: DESCRIPTIVE STATISTICS TABLES\n")
cat(strrep("=", 70), "\n\n")

benchmark_desctable <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    
    vars_for_finalfit <- setdiff(vars_all, "treatment")
    
    microbenchmark(
        summata = desctable(data, by = "treatment", variables = vars_all,
                            stats_continuous = c("mean_sd", "median_iqr"), test = TRUE),
        finalfit = summary_factorlist(df, dependent = "treatment",
                                      explanatory = vars_for_finalfit, p = TRUE, cont = "mean"),
        arsenal = {
            tab <- tableby(treatment ~ age + bmi + biomarker1 + biomarker2 + creatinine + 
                               hemoglobin + sex + race + smoking + hypertension + diabetes + stage + grade,
                           data = df, test = TRUE, numeric.stats = c("meansd", "medianq1q3"))
            as.data.frame(summary(tab, text = TRUE))
        },
        gtsummary = add_p(tbl_summary(df, by = treatment, include = all_of(vars_all),
                                      statistic = list(all_continuous() ~ "{mean} ({sd}), {median} ({p25}, {p75})",
                                                       all_categorical() ~ "{n} ({p}%)"), missing = "no")),
        tableone = CreateTableOne(vars = vars_all, strata = "treatment", data = df, test = TRUE, smd = FALSE),
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Descriptive Table Benchmarks ---\n\n")
bench_desc_1k <- benchmark_desctable(data_1k, df_1k, n_times = 15)
bench_desc_5k <- benchmark_desctable(data_5k, df_5k, n_times = 10)
bench_desc_10k <- benchmark_desctable(data_10k, df_10k, n_times = 5)

bench_desc_all <- rbindlist(list(
    as.data.table(bench_desc_1k)[, n := 1000],
    as.data.table(bench_desc_5k)[, n := 5000],
    as.data.table(bench_desc_10k)[, n := 10000]
))
bench_desc_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest (based on overall median)
bench_desc_all <- order_by_speed(bench_desc_all)

p_desc_box <- ggplot(bench_desc_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(1000, 5000, 10000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Descriptive Tables: desctable() vs Alternatives",
         subtitle = "Execution time in milliseconds (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_desctable.png", p_desc_box, width = 10, height = 5)

desc_scaling <- bench_desc_all[, .(median_time = median(time_ms)), by = .(expr, n)]
p_desc_scale <- ggplot(desc_scaling, aes(x = n, y = median_time, color = expr)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    scale_x_log10(labels = scales::comma) + scale_y_log10() +
    scale_color_manual(values = pkg_colors) +
    labs(title = "Descriptive Tables: Scaling Behavior",
         subtitle = "Median execution time vs dataset size (log-log scale)",
         x = "Number of Observations", y = "Median Time (ms)", color = "Package") +
    theme_benchmark()
save_plot("benchmark_desctable_scaling.png", p_desc_scale, width = 8, height = 5)
cat("Saved: benchmark_desctable.png, benchmark_desctable_scaling.png\n")

## =============================================================================
## BENCHMARK 1B: SURVIVAL SUMMARY TABLES
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 1B: SURVIVAL SUMMARY TABLES\n")
cat(strrep("=", 70), "\n\n")

benchmark_survtable <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    
    microbenchmark(
        summata = survtable(data, outcome = "Surv(os_time, os_status)", by = "treatment",
                            times = c(12, 24, 36), probs = 0.5),
        ## Note: add_p() removed due to gtsummary 2.5.0 compatibility issue
        ## where it cannot extract data from survfit objects created within microbenchmark
        gtsummary = {
            fit <- survfit(Surv(os_time, os_status) ~ treatment, data = df)
            tbl_survfit(fit, times = c(12, 24, 36), label_header = "**{time}**")
        },
        manual = {
            fit <- survfit(Surv(os_time, os_status) ~ treatment, data = df)
            summ <- summary(fit, times = c(12, 24, 36))
            quant <- quantile(fit, probs = 0.5)
            list(survival = summ, median = quant)
        },
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Survival Table Benchmarks ---\n\n")
bench_surv_1k <- benchmark_survtable(data_1k, df_1k, n_times = 15)
bench_surv_5k <- benchmark_survtable(data_5k, df_5k, n_times = 10)
bench_surv_10k <- benchmark_survtable(data_10k, df_10k, n_times = 5)

bench_surv_all <- rbindlist(list(
    as.data.table(bench_surv_1k)[, n := 1000],
    as.data.table(bench_surv_5k)[, n := 5000],
    as.data.table(bench_surv_10k)[, n := 10000]
))
bench_surv_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_surv_all <- order_by_speed(bench_surv_all)

p_surv_box <- ggplot(bench_surv_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(1000, 5000, 10000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Survival Tables: survtable() vs Alternatives",
         subtitle = "Execution time in milliseconds (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_survtable.png", p_surv_box, width = 10, height = 5)
cat("Saved: benchmark_survtable.png\n")

## =============================================================================
## BENCHMARK 2: LOGISTIC REGRESSION OUTPUT
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 2: LOGISTIC REGRESSION OUTPUT\n")
cat(strrep("=", 70), "\n\n")

benchmark_logistic <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    predictors <- c("age", "sex", "bmi", "smoking", "hypertension", "diabetes", "treatment", "stage")
    
    microbenchmark(
        summata = fit(data, outcome = "response", predictors = predictors,
                      model_type = "glm", family = "binomial"),
        summata_minimal = fit(data, outcome = "response", predictors = predictors,
                              model_type = "glm", family = "binomial",
                              show_n = FALSE, show_events = FALSE, reference_rows = FALSE),
        finalfit = fit2df(glmmulti(df, "response", predictors), estimate_suffix = " (multivariable)"),
        gtsummary = tbl_regression(glm(response ~ age + sex + bmi + smoking + hypertension + 
                                           diabetes + treatment + stage, data = df, family = binomial), exponentiate = TRUE),
        broom = tidy(glm(response ~ age + sex + bmi + smoking + hypertension + 
                             diabetes + treatment + stage, data = df, family = binomial), conf.int = TRUE, exponentiate = TRUE),
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Logistic Regression Benchmarks ---\n\n")
bench_log_500 <- benchmark_logistic(data_500, df_500, n_times = 20)
bench_log_1k <- benchmark_logistic(data_1k, df_1k, n_times = 15)
bench_log_5k <- benchmark_logistic(data_5k, df_5k, n_times = 10)
bench_log_10k <- benchmark_logistic(data_10k, df_10k, n_times = 5)

bench_log_all <- rbindlist(list(
    as.data.table(bench_log_500)[, n := 500], as.data.table(bench_log_1k)[, n := 1000],
    as.data.table(bench_log_5k)[, n := 5000], as.data.table(bench_log_10k)[, n := 10000]
))
bench_log_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_log_all <- order_by_speed(bench_log_all)

p_log_box <- ggplot(bench_log_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000, 10000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Logistic Regression: fit() vs Alternatives",
         subtitle = "Execution time in milliseconds (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_logistic.png", p_log_box, width = 10, height = 5)
cat("Saved: benchmark_logistic.png\n")

## =============================================================================
## BENCHMARK 3: LINEAR REGRESSION OUTPUT
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 3: LINEAR REGRESSION OUTPUT\n")
cat(strrep("=", 70), "\n\n")

benchmark_linear <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    predictors <- c("age", "sex", "bmi", "smoking", "hypertension", "diabetes", "treatment", "stage")
    
    microbenchmark(
        summata = fit(data, outcome = "los_days", predictors = predictors, model_type = "lm"),
        summata_minimal = fit(data, outcome = "los_days", predictors = predictors, model_type = "lm",
                              show_n = FALSE, reference_rows = FALSE),
        finalfit = fit2df(lmmulti(df, "los_days", predictors)),
        gtsummary = tbl_regression(lm(los_days ~ age + sex + bmi + smoking + hypertension + 
                                          diabetes + treatment + stage, data = df)),
        broom = tidy(lm(los_days ~ age + sex + bmi + smoking + hypertension + 
                            diabetes + treatment + stage, data = df), conf.int = TRUE),
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Linear Regression Benchmarks ---\n\n")
bench_lm_500 <- benchmark_linear(data_500, df_500, n_times = 20)
bench_lm_1k <- benchmark_linear(data_1k, df_1k, n_times = 15)
bench_lm_5k <- benchmark_linear(data_5k, df_5k, n_times = 10)
bench_lm_10k <- benchmark_linear(data_10k, df_10k, n_times = 5)

bench_lm_all <- rbindlist(list(
    as.data.table(bench_lm_500)[, n := 500], as.data.table(bench_lm_1k)[, n := 1000],
    as.data.table(bench_lm_5k)[, n := 5000], as.data.table(bench_lm_10k)[, n := 10000]
))
bench_lm_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_lm_all <- order_by_speed(bench_lm_all)

p_lm_box <- ggplot(bench_lm_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000, 10000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Linear Regression: fit() vs Alternatives",
         subtitle = "Execution time in milliseconds (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_linear.png", p_lm_box, width = 10, height = 5)
cat("Saved: benchmark_linear.png\n")

## =============================================================================
## BENCHMARK 4: POISSON REGRESSION OUTPUT
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 4: POISSON REGRESSION OUTPUT\n")
cat(strrep("=", 70), "\n\n")

benchmark_poisson <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    predictors <- c("age", "sex", "bmi", "smoking", "treatment", "stage")
    
    microbenchmark(
        summata = fit(data, outcome = "event_count", predictors = predictors,
                      model_type = "glm", family = "poisson"),
        summata_minimal = fit(data, outcome = "event_count", predictors = predictors,
                              model_type = "glm", family = "poisson",
                              show_n = FALSE, show_events = FALSE, reference_rows = FALSE),
        finalfit = fit2df(glm(event_count ~ age + sex + bmi + smoking + treatment + stage,
                              data = df, family = poisson)),
        gtsummary = tbl_regression(glm(event_count ~ age + sex + bmi + smoking + treatment + stage,
                                       data = df, family = poisson), exponentiate = TRUE),
        broom = tidy(glm(event_count ~ age + sex + bmi + smoking + treatment + stage,
                         data = df, family = poisson), conf.int = TRUE, exponentiate = TRUE),
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Poisson Regression Benchmarks ---\n\n")
bench_pois_500 <- benchmark_poisson(data_500, df_500, n_times = 20)
bench_pois_1k <- benchmark_poisson(data_1k, df_1k, n_times = 15)
bench_pois_5k <- benchmark_poisson(data_5k, df_5k, n_times = 10)
bench_pois_10k <- benchmark_poisson(data_10k, df_10k, n_times = 5)

bench_pois_all <- rbindlist(list(
    as.data.table(bench_pois_500)[, n := 500], as.data.table(bench_pois_1k)[, n := 1000],
    as.data.table(bench_pois_5k)[, n := 5000], as.data.table(bench_pois_10k)[, n := 10000]
))
bench_pois_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_pois_all <- order_by_speed(bench_pois_all)

p_pois_box <- ggplot(bench_pois_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000, 10000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Poisson Regression: fit() vs Alternatives",
         subtitle = "Execution time in milliseconds (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_poisson.png", p_pois_box, width = 10, height = 5)
cat("Saved: benchmark_poisson.png\n")

## =============================================================================
## BENCHMARK 5: COX REGRESSION
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 5: COX PROPORTIONAL HAZARDS REGRESSION\n")
cat(strrep("=", 70), "\n\n")

benchmark_cox <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    predictors <- c("age", "sex", "bmi", "treatment", "stage")
    
    microbenchmark(
        summata = fit(data, outcome = "Surv(os_time, os_status)", predictors = predictors, model_type = "coxph"),
        summata_minimal = fit(data, outcome = "Surv(os_time, os_status)", predictors = predictors, model_type = "coxph",
                              show_n = FALSE, show_events = FALSE, reference_rows = FALSE),
        finalfit = fit2df(coxphmulti(df, "Surv(os_time, os_status)", predictors)),
        gtsummary = tbl_regression(coxph(Surv(os_time, os_status) ~ age + sex + bmi + treatment + stage,
                                         data = df), exponentiate = TRUE),
        broom = tidy(coxph(Surv(os_time, os_status) ~ age + sex + bmi + treatment + stage,
                           data = df), conf.int = TRUE, exponentiate = TRUE),
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Cox Regression Benchmarks ---\n\n")
bench_cox_500 <- benchmark_cox(data_500, df_500, n_times = 15)
bench_cox_1k <- benchmark_cox(data_1k, df_1k, n_times = 10)
bench_cox_5k <- benchmark_cox(data_5k, df_5k, n_times = 5)

bench_cox_all <- rbindlist(list(
    as.data.table(bench_cox_500)[, n := 500], as.data.table(bench_cox_1k)[, n := 1000],
    as.data.table(bench_cox_5k)[, n := 5000]
))
bench_cox_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_cox_all <- order_by_speed(bench_cox_all)

p_cox_box <- ggplot(bench_cox_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Cox Regression: fit() vs Alternatives",
         subtitle = "Execution time in milliseconds (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_cox.png", p_cox_box, width = 9, height = 5)
cat("Saved: benchmark_cox.png\n")

## =============================================================================
## BENCHMARK 6: MIXED-EFFECTS MODELS
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 6: MIXED-EFFECTS MODELS\n")
cat(strrep("=", 70), "\n\n")

benchmark_mixed <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    
    microbenchmark(
        summata = fit(data, outcome = "los_days", predictors = c("age", "sex", "treatment", "(1|site)"),
                      model_type = "lmer"),
        summata_minimal = fit(data, outcome = "los_days", predictors = c("age", "sex", "treatment", "(1|site)"),
                              model_type = "lmer", show_n = FALSE, show_events = FALSE, reference_rows = FALSE),
        finalfit = fit2df(lmmixed(df, "los_days", c("age", "sex", "treatment"), "(1|site)")),
        gtsummary = tbl_regression(lmer(los_days ~ age + sex + treatment + (1|site), data = df)),
        broom = tidy(lmer(los_days ~ age + sex + treatment + (1|site), data = df), conf.int = TRUE),
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Mixed-Effects Benchmarks ---\n\n")
bench_mixed_500 <- benchmark_mixed(data_500, df_500, n_times = 15)
bench_mixed_1k <- benchmark_mixed(data_1k, df_1k, n_times = 10)
bench_mixed_5k <- benchmark_mixed(data_5k, df_5k, n_times = 5)

bench_mixed_all <- rbindlist(list(
    as.data.table(bench_mixed_500)[, n := 500], as.data.table(bench_mixed_1k)[, n := 1000],
    as.data.table(bench_mixed_5k)[, n := 5000]
))
bench_mixed_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_mixed_all <- order_by_speed(bench_mixed_all)

p_mixed_box <- ggplot(bench_mixed_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Mixed-Effects Models: fit() with lmer vs Alternatives",
         subtitle = "Linear mixed model with random intercepts (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_mixed.png", p_mixed_box, width = 9, height = 5)
cat("Saved: benchmark_mixed.png\n")

## =============================================================================
## BENCHMARK 7: UNIVARIABLE SCREENING
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 7: UNIVARIABLE SCREENING (Multiple Predictors)\n")
cat(strrep("=", 70), "\n\n")

benchmark_uniscreen <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    predictors <- c("age", "sex", "bmi", "race", "smoking", "hypertension", 
                    "diabetes", "treatment", "stage", "grade", 
                    "biomarker1", "biomarker2", "creatinine", "hemoglobin")
    
    microbenchmark(
        summata = uniscreen(data, outcome = "response", predictors = predictors,
                            model_type = "glm", family = "binomial"),
        summata_minimal = uniscreen(data, outcome = "response", predictors = predictors,
                                    model_type = "glm", family = "binomial",
                                    show_n = FALSE, show_events = FALSE, reference_rows = FALSE),
        finalfit = fit2df(glmuni(df, "response", predictors), estimate_suffix = " (univariable)"),
        gtsummary = tbl_uvregression(df[, c("response", predictors)], method = glm, y = response,
                                     method.args = list(family = binomial), exponentiate = TRUE),
        arsenal = modelsum(response ~ age + sex + bmi + race + smoking + hypertension + 
                               diabetes + treatment + stage + grade + biomarker1 + biomarker2 + creatinine + hemoglobin,
                           data = df, family = binomial),
        broom_loop = rbindlist(lapply(predictors, function(pred) {
            tidy(glm(as.formula(paste("response ~", pred)), data = df, family = binomial),
                 conf.int = TRUE, exponentiate = TRUE)
        })),
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Univariable Screening Benchmarks ---\n\n")
bench_uni_500 <- benchmark_uniscreen(data_500, df_500, n_times = 15)
bench_uni_1k <- benchmark_uniscreen(data_1k, df_1k, n_times = 10)
bench_uni_5k <- benchmark_uniscreen(data_5k, df_5k, n_times = 5)

bench_uni_all <- rbindlist(list(
    as.data.table(bench_uni_500)[, n := 500], as.data.table(bench_uni_1k)[, n := 1000],
    as.data.table(bench_uni_5k)[, n := 5000]
))
bench_uni_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_uni_all <- order_by_speed(bench_uni_all)

p_uni_box <- ggplot(bench_uni_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Univariable Screening: uniscreen() vs Alternatives",
         subtitle = "Screening 14 predictors (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_uniscreen.png", p_uni_box, width = 10, height = 5)
cat("Saved: benchmark_uniscreen.png\n")

## =============================================================================
## BENCHMARK 8: COMPLETE WORKFLOW (Univariable + Multivariable)
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 8: COMPLETE ANALYSIS WORKFLOW\n")
cat(strrep("=", 70), "\n\n")

benchmark_workflow <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    predictors <- c("age", "sex", "bmi", "smoking", "hypertension", "diabetes", "treatment", "stage")
    
    microbenchmark(
        summata = fullfit(data, outcome = "response", predictors = predictors, method = "all",
                          model_type = "glm", family = "binomial", columns = "both"),
        summata_minimal = fullfit(data, outcome = "response", predictors = predictors, method = "all",
                                  model_type = "glm", family = "binomial", columns = "both",
                                  show_n = FALSE, show_events = FALSE, reference_rows = FALSE),
        finalfit = finalfit(df, dependent = "response", explanatory = predictors),
        gtsummary = {
            uni_table <- tbl_uvregression(df[, c("response", predictors)], method = glm, y = response,
                                          method.args = list(family = binomial), exponentiate = TRUE)
            multi_model <- glm(response ~ age + sex + bmi + smoking + hypertension + 
                                   diabetes + treatment + stage, data = df, family = binomial)
            multi_table <- tbl_regression(multi_model, exponentiate = TRUE)
            tbl_merge(tbls = list(uni_table, multi_table), tab_spanner = c("Univariable", "Multivariable"))
        },
        manual = {
            uni_dt <- rbindlist(lapply(predictors, function(pred) {
                res <- tidy(glm(as.formula(paste("response ~", pred)), data = df, family = binomial),
                            conf.int = TRUE, exponentiate = TRUE)
                res$predictor <- pred
                res
            }))
            multi_formula <- as.formula(paste("response ~", paste(predictors, collapse = " + ")))
            multi_dt <- as.data.table(tidy(glm(multi_formula, data = df, family = binomial),
                                           conf.int = TRUE, exponentiate = TRUE))
            list(univariable = uni_dt, multivariable = multi_dt)
        },
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Complete Workflow Benchmarks ---\n\n")
bench_wf_500 <- benchmark_workflow(data_500, df_500, n_times = 15)
bench_wf_1k <- benchmark_workflow(data_1k, df_1k, n_times = 10)
bench_wf_5k <- benchmark_workflow(data_5k, df_5k, n_times = 5)

bench_wf_all <- rbindlist(list(
    as.data.table(bench_wf_500)[, n := 500], as.data.table(bench_wf_1k)[, n := 1000],
    as.data.table(bench_wf_5k)[, n := 5000]
))
bench_wf_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_wf_all <- order_by_speed(bench_wf_all)

p_wf_box <- ggplot(bench_wf_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Complete Workflow: fullfit() vs Alternatives",
         subtitle = "Combined univariable + multivariable analysis (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_workflow.png", p_wf_box, width = 10, height = 5)
cat("Saved: benchmark_workflow.png\n")

## =============================================================================
## BENCHMARK 9: FOREST PLOTS
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK 9: FOREST PLOT GENERATION\n")
cat(strrep("=", 70), "\n\n")

benchmark_forest <- function(data, df, n_times = 10) {
    cat(sprintf("Testing with %s observations...\n", format(nrow(data), big.mark = ",")))
    cox_model <- coxph(Surv(os_time, os_status) ~ age + sex + bmi + treatment + stage, data = df)
    
    microbenchmark(
        summata = coxforest(x = cox_model, data = df, title = "Cox Regression Results"),
        survminer = ggforest(cox_model, data = df),
        manual = {
            coef_data <- as.data.table(tidy(cox_model, conf.int = TRUE, exponentiate = TRUE))
            coef_data <- coef_data[term != "(Intercept)"]
            ggplot(coef_data, aes(x = estimate, y = term)) +
                geom_point() +
                geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
                geom_vline(xintercept = 1, linetype = "dashed") +
                scale_x_log10() + theme_minimal()
        },
        times = n_times, unit = "ms"
    )
}

cat("\n--- Running Forest Plot Benchmarks ---\n\n")
bench_forest_500 <- benchmark_forest(data_500, df_500, n_times = 15)
bench_forest_1k <- benchmark_forest(data_1k, df_1k, n_times = 10)
bench_forest_5k <- benchmark_forest(data_5k, df_5k, n_times = 5)

bench_forest_all <- rbindlist(list(
    as.data.table(bench_forest_500)[, n := 500], as.data.table(bench_forest_1k)[, n := 1000],
    as.data.table(bench_forest_5k)[, n := 5000]
))
bench_forest_all[, time_ms := time / 1e6]

## Order packages from fastest to slowest
bench_forest_all <- order_by_speed(bench_forest_all)

p_forest_box <- ggplot(bench_forest_all, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~ factor(n, labels = paste0("n = ", scales::comma(c(500, 1000, 5000)))), 
               scales = "free_y", nrow = 1) +
    scale_y_log10() + scale_fill_manual(values = pkg_colors) +
    labs(title = "Forest Plots: coxforest() vs Alternatives",
         subtitle = "Execution time in milliseconds (log scale)",
         x = NULL, y = "Time (ms)", fill = "Package") +
    theme_benchmark() + guides(fill = guide_legend(nrow = 1))
save_plot("benchmark_forest.png", p_forest_box, width = 9, height = 5)
cat("Saved: benchmark_forest.png\n")

## =============================================================================
## SUMMARY VISUALIZATION: RELATIVE PERFORMANCE
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("GENERATING SUMMARY VISUALIZATIONS\n")
cat(strrep("=", 70), "\n\n")

## Combine all benchmarks with labels in vignette order:
## Descriptive tables, Survival tables, Univariable screen, Logistic regression,
## Poisson regression, Linear regression, Cox regression, Mixed effects,
## Complete workflow, Forest plots
all_benchmarks <- rbindlist(list(
    bench_desc_all[, benchmark := "Descriptive Tables"],
    bench_surv_all[, benchmark := "Survival Tables"],
    bench_uni_all[, benchmark := "Univariable Screening"],
    bench_log_all[, benchmark := "Logistic Regression"],
    bench_pois_all[, benchmark := "Poisson Regression"],
    bench_lm_all[, benchmark := "Linear Regression"],
    bench_cox_all[, benchmark := "Cox Regression"],
    bench_mixed_all[, benchmark := "Mixed-Effects"],
    bench_wf_all[, benchmark := "Complete Workflow"],
    bench_forest_all[, benchmark := "Forest Plots"]
), fill = TRUE)

## Set benchmark factor order to match vignette
benchmark_order <- c("Descriptive Tables", "Survival Tables", "Univariable Screening",
                     "Logistic Regression", "Poisson Regression", "Linear Regression",
                     "Cox Regression", "Mixed-Effects", "Complete Workflow", "Forest Plots")
all_benchmarks[, benchmark := factor(benchmark, levels = benchmark_order)]

## Exclude summata_minimal from speedup calculations
speedup_data <- all_benchmarks[!expr %in% c("summata", "summata_minimal")]
summata_times <- all_benchmarks[expr == "summata", .(summata_time = median(time_ms)), by = .(benchmark, n)]
speedup_summary <- speedup_data[, .(median_ms = median(time_ms)), by = .(benchmark, n, expr)]
speedup_summary <- speedup_summary[summata_times, on = .(benchmark, n)]
speedup_summary[, speedup := median_ms / summata_time]

p_speedup <- ggplot(speedup_summary, aes(x = factor(n), y = speedup, fill = expr)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
    facet_wrap(~ benchmark, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = pkg_colors) +
    labs(title = "Relative Performance: Time Relative to summata",
         subtitle = "Values > 1 indicate summata is faster (red line = equal performance)",
         x = "Dataset Size (n)", y = "Relative Time", fill = "Package") +
    theme_benchmark() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))
save_plot("benchmark_speedup.png", p_speedup, width = 10, height = 14)
cat("Saved: benchmark_speedup.png\n")

ff_comparison <- all_benchmarks[expr %in% c("summata", "finalfit"), 
                                .(median_time = median(time_ms)), by = .(benchmark, expr, n)]

p_ff <- ggplot(ff_comparison, aes(x = factor(n), y = median_time, fill = expr)) +
    geom_col(position = "dodge", alpha = 0.8) +
    facet_wrap(~ benchmark, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("summata" = "#FC8D62", "finalfit" = "#5597DF")) +
    labs(title = "summata vs finalfit: Direct Comparison",
         subtitle = "summata builds upon finalfit with performance optimizations",
         x = "Dataset Size (n)", y = "Median Time (ms)", fill = "Package") +
    theme_benchmark() +
    theme(legend.position = "bottom")
save_plot("benchmark_summata_vs_finalfit.png", p_ff, width = 10, height = 14)
cat("Saved: benchmark_summata_vs_finalfit.png\n")

## =============================================================================
## SUMMARY TABLES
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("FINAL SUMMARY TABLES\n")
cat(strrep("=", 70), "\n\n")

final_summary <- dcast(
    all_benchmarks[, .(median_ms = round(median(time_ms), 2)), by = .(benchmark, n, expr)],
    benchmark + n ~ expr, value.var = "median_ms"
)
setorder(final_summary, benchmark, n)
cat("Median execution time in milliseconds:\n\n")
print(final_summary)

## Exclude summata_minimal from speedup table
speedup_for_table <- all_benchmarks[!expr %in% c("summata_minimal")]
speedup_table <- dcast(
    speedup_for_table[, .(median_ms = median(time_ms)), by = .(benchmark, n, expr)][
      , summata_time := median_ms[expr == "summata"], by = .(benchmark, n)][
      , speedup := round(median_ms / summata_time, 2)],
    benchmark + n ~ expr, value.var = "speedup"
)
setorder(speedup_table, benchmark, n)
cat("\n\nSpeedup factors (alternative / summata, >1 means summata is faster):\n\n")
print(speedup_table)

fwrite(final_summary, "benchmark_summary_times.csv")
fwrite(speedup_table, "benchmark_summary_speedup.csv")
cat("\nSaved: benchmark_summary_times.csv, benchmark_summary_speedup.csv\n")

## =============================================================================
## SESSION INFO
## =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("SESSION INFORMATION\n")
cat(strrep("=", 70), "\n\n")
sessionInfo()

cat("\n\nBenchmarking complete!\n\nOutput files:\n")
cat(" - benchmark_desctable.png\n - benchmark_desctable_scaling.png\n")
cat(" - benchmark_survtable.png\n")
cat(" - benchmark_logistic.png\n - benchmark_linear.png\n - benchmark_poisson.png\n")
cat(" - benchmark_cox.png\n - benchmark_mixed.png\n - benchmark_uniscreen.png\n")
cat(" - benchmark_workflow.png\n - benchmark_forest.png\n")
cat(" - benchmark_speedup.png\n - benchmark_summata_vs_finalfit.png\n")
cat(" - benchmark_summary_times.csv\n - benchmark_summary_speedup.csv\n")
