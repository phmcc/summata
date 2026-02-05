#' Simulated Clinical Trial Dataset
#'
#' A simulated dataset representing a multi-site oncology clinical trial with 
#' 850 patients. The dataset includes realistic correlations, interaction effects,
#' site-level clustering, and a comprehensive set of outcomes suitable for 
#' demonstrating regression analysis, survival modeling, and multi-outcome analyses.
#'
#' @format A data frame with 850 rows and 32 variables:
#' \describe{
#'   \item{patient_id}{Character. Unique patient identifier (PT0001-PT0850)}
#'   
#'   \strong{Demographics:}
#'   \item{age}{Numeric. Age in years (18-90)}
#'   \item{sex}{Factor. Sex (Female, Male)}
#'   \item{race}{Factor. Race (White, Black, Asian, Other)}
#'   \item{ethnicity}{Factor. Ethnicity (Non-Hispanic, Hispanic)}
#'   
#'   \strong{Clinical Characteristics:}
#'   \item{bmi}{Numeric. Body mass index in kg/m² (15-50)}
#'   \item{smoking}{Factor. Smoking status (Never, Former, Current)}
#'   \item{hypertension}{Factor. Hypertension diagnosis (No, Yes)}
#'   \item{diabetes}{Factor. Diabetes diagnosis (No, Yes)}
#'   \item{ecog}{Factor. ECOG performance status (0, 1, 2, 3)}
#'   
#'   \strong{Laboratory Values:}
#'   \item{creatinine}{Numeric. Baseline serum creatinine in mg/dL}
#'   \item{hemoglobin}{Numeric. Baseline hemoglobin in g/dL (7-18)}
#'   \item{biomarker_x}{Numeric. Tumor biomarker X in ng/mL}
#'   \item{biomarker_y}{Numeric. Tumor biomarker Y in U/L}
#'   
#'   \strong{Disease Characteristics:}
#'   \item{site}{Factor. Study site (Site Alpha through Site Kappa)}
#'   \item{grade}{Factor. Tumor grade (Well-differentiated, Moderately differentiated, 
#'     Poorly differentiated)}
#'   \item{stage}{Factor. Disease stage (I, II, III, IV)}
#'   
#'   \strong{Treatment:}
#'   \item{treatment}{Factor. Treatment arm (Control, Drug A, Drug B)}
#'   \item{surgery}{Factor. Surgical resection performed (No, Yes)}
#'   
#'   \strong{Hospital Outcomes:}
#'   \item{los_days}{Numeric. Length of hospital stay in days}
#'   
#'   \strong{Complication Outcomes:}
#'   \item{any_complication}{Factor. Any complication occurred (No, Yes)}
#'   \item{wound_infection}{Factor. Wound or site infection (No, Yes)}
#'   \item{readmission_30d}{Factor. 30-day hospital readmission (No, Yes)}
#'   \item{icu_admission}{Factor. ICU admission required (No, Yes)}
#'   
#'   \strong{Recovery Outcomes:}
#'   \item{pain_score}{Numeric. Postoperative pain score on 0-10 visual analog scale}
#'   \item{recovery_days}{Numeric. Days to functional recovery}
#'   
#'   \strong{Count Outcomes:}
#'   \item{ae_count}{Integer. Number of adverse events during study period. 
#'     Overdispersed count suitable for negative binomial or quasipoisson regression.
#'     Associated with age, ECOG, diabetes, treatment, surgery, and stage.}
#'   \item{fu_count}{Integer. Number of follow-up clinic visits. Equidispersed
#'     count suitable for standard Poisson regression. Associated with stage, ECOG,
#'     treatment, and surgery.}
#'   
#'   \strong{Survival Outcomes:}
#'   \item{pfs_months}{Numeric. Progression-free survival time in months}
#'   \item{pfs_status}{Numeric. PFS event indicator (0=censored, 1=progression/death)}
#'   \item{os_months}{Numeric. Overall survival time in months}
#'   \item{os_status}{Numeric. OS event indicator (0=censored, 1=death)}
#' }
#'
#' @details
#' The dataset is designed to demonstrate various statistical modeling scenarios:
#' 
#' \strong{Treatment Effects:}
#' \itemize{
#'   \item Drug A generally shows protective effects (reduced complications, 
#'     faster recovery, improved survival, fewer adverse events)
#'   \item Drug B shows increased toxicity (more complications, longer recovery,
#'     more adverse events)
#'   \item Treatment effects vary by disease stage (interaction effects)
#' }
#' 
#' \strong{Count Outcomes:}
#' \itemize{
#'   \item \code{ae_count}: Generated using negative binomial distribution with
#'     overdispersion (variance > mean). Use for demonstrating quasipoisson or
#'     negative binomial regression.
#'   \item \code{fu_count}: Generated using Poisson distribution (variance ≈ mean).
#'     Use for demonstrating standard Poisson regression.
#' }
#' 
#' \strong{Site Clustering:}
#' \itemize{
#'   \item Outcomes are correlated within sites (random intercepts)
#'   \item Treatment effects vary by site (random slopes)
#'   \item Suitable for mixed-effects modeling with \code{(1|site)} random effects
#' }
#' 
#' \strong{Key Associations:}
#' \itemize{
#'   \item Diabetes is strongly associated with complications, infections, and readmission
#'   \item Age and ECOG performance status affect most outcomes
#'   \item Surgery increases complication risk but may improve survival
#'   \item 30-day readmission is associated with age, sex, diabetes, smoking, 
#'     ECOG status, disease stage, treatment arm, and surgery
#'   \item Multiple interaction effects (age × diabetes, treatment × stage, etc.)
#' }
#' 
#' \strong{Missing Data:}
#' \itemize{
#'   \item Realistic missing data patterns (~2\% across most variables)
#'   \item Some outcomes have outcome-specific missingness (e.g., LOS for 
#'     outpatient procedures)
#' }
#'
#' @source Simulated data generated for package demonstration purposes.
#'
#' @seealso 
#' \code{\link{clintrial_labels}} for variable labels suitable for tables and plots
#'
#' @examples
#' # Load the dataset
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Basic exploration
#' str(clintrial)
#' summary(clintrial$treatment)
#' 
#' # Logistic regression with forest plot
#' model <- glm(os_status ~ age + sex + treatment + stage, 
#'              data = clintrial, family = binomial)
#' glmforest(model, data = clintrial, labels = clintrial_labels)
#' 
#' # Multi-outcome analysis
#' multi <- multifit(
#'     data = clintrial,
#'     outcomes = c("postop_complication", "wound_infection", 
#'                  "readmission_30d", "icu_admission"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "diabetes"),
#'     labels = clintrial_labels
#' )
#' multiforest(multi)
#' 
#' # Mixed-effects model accounting for site clustering
#' library(lme4)
#' me_model <- glmer(os_status ~ age + treatment + (1|site), 
#'                   data = clintrial, family = binomial)
#' glmforest(me_model, data = clintrial, labels = clintrial_labels)
#' 
#' # Survival analysis
#' library(survival)
#' cox_model <- coxph(Surv(os_months, os_status) ~ age + sex + treatment + stage,
#'                    data = clintrial)
#' coxforest(cox_model, data = clintrial, labels = clintrial_labels)
#'
#' @name clintrial
#' @docType data
#' @keywords datasets
NULL

set.seed(71)
n <- 850

create_clintrial <- function() {

    # ============================================================================
    # SITE AND CLUSTERING STRUCTURE
    # ============================================================================
    
    # Site assignment
    greek_letters <- c("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Iota", "Kappa")
    site <- factor(
        paste0("Site ", sample(greek_letters, n, replace = TRUE)),
        levels = paste0("Site ", greek_letters))
    
    # Site-level random effects
    n_sites <- length(unique(site))
    site_names <- levels(site)
    site_intercepts <- rnorm(n_sites, mean = 0, sd = 0.5)
    site_drugA_effect <- rnorm(n_sites, mean = 0, sd = 0.3)
    site_drugB_effect <- rnorm(n_sites, mean = 0, sd = 0.3)
    site_age_slope <- rnorm(n_sites, mean = 1, sd = 0.2)
    names(site_intercepts) <- site_names
    names(site_drugA_effect) <- site_names
    names(site_drugB_effect) <- site_names
    names(site_age_slope) <- site_names
    
    # ============================================================================
    # DEMOGRAPHICS
    # ============================================================================
    
    # Age
    age <- round(rnorm(n, mean = 60, sd = 12))
    age[age < 18] <- 18
    age[age > 90] <- 90
    
    # Sex
    sex <- factor(sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.48, 0.52)),
                  levels = c("Female", "Male"))
    
    # Race
    race <- factor(sample(c("White", "Black", "Asian", "Other"), n, replace = TRUE,
                          prob = c(0.65, 0.20, 0.10, 0.05)),
                   levels = c("White", "Black", "Asian", "Other"))
    
    # Ethnicity
    ethnicity <- factor(sample(c("Non-Hispanic", "Hispanic"), n, replace = TRUE, prob = c(0.85, 0.15)),
                        levels = c("Non-Hispanic", "Hispanic"))
    
    # ============================================================================
    # CLINICAL CHARACTERISTICS
    # ============================================================================
    
    # BMI (correlated with age and diabetes risk)
    bmi <- round(rnorm(n, mean = 27 + (age - 60) * 0.05, sd = 5), 1)
    bmi[bmi < 15] <- 15
    bmi[bmi > 50] <- 50
    
    # Smoking
    smoke_prob <- ifelse(sex == "Male", 0.35, 0.25)
    smoking <- factor(
        sample(c("Never", "Former", "Current"), n, replace = TRUE,
               prob = c(0.45, 0.35, 0.20)),
        levels = c("Never", "Former", "Current"))
    
    # Hypertension (correlated with age and BMI)
    htn_prob <- pmin(0.9, plogis(-2 + 0.04 * age + 0.03 * bmi))
    hypertension <- factor(ifelse(runif(n) < htn_prob, "Yes", "No"), levels = c("No", "Yes"))
    
    # Diabetes (correlated with age, BMI, and race)
    dm_prob <- pmin(0.9, plogis(-3.5 + 0.03 * age + 0.08 * bmi +
                                    (race == "Black") * 0.5 + (race == "Asian") * 0.3))
    diabetes <- factor(ifelse(runif(n) < dm_prob, "Yes", "No"), levels = c("No", "Yes"))
    
    # ECOG performance status (correlated with age)
    ecog_base_prob <- plogis(-2 + 0.03 * age)
    ecog_probs <- cbind(
        1 - ecog_base_prob * 1.5,
        ecog_base_prob * 0.8,
        ecog_base_prob * 0.5,
        ecog_base_prob * 0.2
    )
    ecog_probs <- ecog_probs / rowSums(ecog_probs)
    ecog <- factor(apply(ecog_probs, 1, function(p) sample(0:3, 1, prob = p)), levels = 0:3)
    
    # ============================================================================
    # LABORATORY VALUES
    # ============================================================================
    
    # Creatinine (correlated with age and diabetes)
    creatinine <- round(pmax(0.5, rnorm(n, mean = 0.9 + 0.008 * age + 
                                            (diabetes == "Yes") * 0.3 + 
                                            (sex == "Male") * 0.15, sd = 0.25)), 2)
    
    # Hemoglobin (correlated with sex)
    hemoglobin <- round(rnorm(n, mean = ifelse(sex == "Male", 14.5, 13.0), sd = 1.5), 1)
    hemoglobin[hemoglobin < 7] <- 7
    hemoglobin[hemoglobin > 18] <- 18
    
    # Biomarker X (log-normal distribution, correlated with stage)
    biomarker_x_base <- rlnorm(n, meanlog = 2, sdlog = 0.8)
    biomarker_x <- round(biomarker_x_base, 1)
    
    # Biomarker Y (correlated with age and will be correlated with stage)
    biomarker_y_base <- round(rlnorm(n, meanlog = 3 + 0.01 * age, sdlog = 0.6), 1)
    biomarker_y <- biomarker_y_base
    
    # ============================================================================
    # DISEASE CHARACTERISTICS
    # ============================================================================
    
    # Tumor grade (correlated with biomarker X)
    grade_prob_well <- pmax(0.1, 0.6 - 0.05 * log(biomarker_x))
    grade_prob_poor <- pmin(0.5, 0.1 + 0.04 * log(biomarker_x))
    grade_prob_mod <- 1 - grade_prob_well - grade_prob_poor
    grade_probs <- cbind(grade_prob_well, grade_prob_mod, grade_prob_poor)
    grade <- factor(
        apply(grade_probs, 1, function(p) {
            sample(c("Well-differentiated", "Moderately differentiated", "Poorly differentiated"),
                   1, prob = p)
        }),
        levels = c("Well-differentiated", "Moderately differentiated", "Poorly differentiated")
    )
    
    # Disease stage (correlated with grade, age, biomarkers)
    stage_logit <- -2 +
        (grade == "Moderately differentiated") * 0.5 +
        (grade == "Poorly differentiated") * 1.2 +
        0.02 * age +
        0.15 * log(biomarker_x) +
        0.1 * log(biomarker_y)
    stage_prob <- plogis(stage_logit)
    stage_cats <- cut(stage_prob, breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                      labels = c("I", "II", "III", "IV"), include.lowest = TRUE)
    stage <- factor(stage_cats, levels = c("I", "II", "III", "IV"))
    
    # Update biomarker Y to reflect final stage assignment
    biomarker_y <- round(biomarker_y_base * 
                             ifelse(stage == "I", 0.7,
                                    ifelse(stage == "II", 1.0,
                                           ifelse(stage == "III", 1.4, 2.0))), 1)
    
    # ============================================================================
    # TREATMENT ASSIGNMENT
    # ============================================================================
    
    # Treatment (stratified randomization by stage)
    treatment_list <- lapply(levels(stage), function(s) {
        n_stage <- sum(stage == s, na.rm = TRUE)
        sample(c("Control", "Drug A", "Drug B"), n_stage, replace = TRUE)
    })
    treatment <- factor(unsplit(treatment_list, stage), levels = c("Control", "Drug A", "Drug B"))
    
    # Surgery (correlated with stage and ECOG)
    surgery_prob <- plogis(-1 + (stage == "I") * 1.5 + (stage == "II") * 1.0 +
                               (stage == "III") * 0.3 - as.numeric(ecog) * 0.4)
    surgery <- factor(ifelse(runif(n) < surgery_prob, "Yes", "No"), levels = c("No", "Yes"))
    
    # ============================================================================
    # ICU ADMISSION
    # ============================================================================
    
    icu_logit <- -3.5 +
        0.04 * age +
        (sex == "Male") * 0.2 +
        (diabetes == "Yes") * 0.8 +
        (smoking == "Current") * 0.5 +
        as.numeric(ecog) * 0.6 +
        (stage == "III") * 0.4 +
        (stage == "IV") * 0.8 +
        (treatment == "Drug A") * (-0.4) +
        (treatment == "Drug B") * 0.3 +
        (surgery == "Yes") * 1.2 +
        site_intercepts[as.character(site)] * 0.5
    
    icu_prob <- plogis(icu_logit)
    icu_admission <- factor(ifelse(runif(n) < icu_prob, "Yes", "No"), levels = c("No", "Yes"))
    
    # ============================================================================
    # PAIN SCORE
    # ============================================================================
    
    pain_score_mean <- 3 +
        (surgery == "Yes") * 2.5 +
        (sex == "Female") * 0.5 +
        as.numeric(ecog) * 0.3 +
        (treatment == "Drug A") * (-0.5) +
        (treatment == "Drug B") * 0.4 +
        site_intercepts[as.character(site)] * 0.3
    
    pain_score <- round(pmax(0, pmin(10, rnorm(n, mean = pain_score_mean, sd = 1.5))), 1)
    
    # ============================================================================
    # WOUND INFECTION
    # ============================================================================
    
    wound_logit <- -4.0 +
        0.03 * age +
        (diabetes == "Yes") * 1.5 +
        (smoking == "Current") * 0.7 +
        0.05 * bmi +
        (surgery == "Yes") * 1.8 +
        (treatment == "Drug A") * (-0.3) +
        (treatment == "Drug B") * 0.5 +
        site_intercepts[as.character(site)] * 0.4
    
    wound_prob <- plogis(wound_logit)
    wound_infection <- factor(ifelse(runif(n) < wound_prob, "Yes", "No"), levels = c("No", "Yes"))
    
    # ============================================================================
    # ANY COMPLICATION
    # ============================================================================
    
    comp_logit <- -3.2 +
        0.04 * age +
        (diabetes == "Yes") * 1.2 +
        (hypertension == "Yes") * 0.5 +
        (smoking == "Current") * 0.6 +
        as.numeric(ecog) * 0.7 +
        (stage == "III") * 0.3 +
        (stage == "IV") * 0.6 +
        (surgery == "Yes") * 1.5 +
        (treatment == "Drug A") * (-0.5) +
        (treatment == "Drug B") * 0.6 +
        (wound_infection == "Yes") * 2.0 +
        (icu_admission == "Yes") * 1.5 +
        site_intercepts[as.character(site)] * 0.6
    
    comp_prob <- plogis(comp_logit)
    any_complication <- factor(ifelse(runif(n) < comp_prob, "Yes", "No"), levels = c("No", "Yes"))
    
    # ============================================================================
    # 30-DAY READMISSION
    # ============================================================================
    
    readm_logit <- -3.8 +
        0.03 * age +
        (sex == "Male") * 0.3 +
        (diabetes == "Yes") * 1.0 +
        (smoking == "Current") * 0.5 +
        as.numeric(ecog) * 0.5 +
        (stage == "III") * 0.4 +
        (stage == "IV") * 0.8 +
        (treatment == "Drug A") * (-0.2) +
        (treatment == "Drug B") * 0.4 +
        (surgery == "Yes") * 0.6 +
        (any_complication == "Yes") * 1.8 +
        (wound_infection == "Yes") * 1.2 +
        site_intercepts[as.character(site)] * 0.5
    
    readm_prob <- plogis(readm_logit)
    readmission_30d <- factor(ifelse(runif(n) < readm_prob, "Yes", "No"), levels = c("No", "Yes"))
    
    # ============================================================================
    # LENGTH OF STAY
    # ============================================================================
    
    los_mean <- exp(1.5 +
                        (surgery == "Yes") * 0.6 +
                        0.008 * age +
                        (diabetes == "Yes") * 0.3 +
                        as.numeric(ecog) * 0.2 +
                        (any_complication == "Yes") * 0.5 +
                        (icu_admission == "Yes") * 0.7 +
                        (treatment == "Drug A") * (-0.15) +
                        (treatment == "Drug B") * 0.2 +
                        site_intercepts[as.character(site)] * 0.3)
    
    los_days <- round(pmax(0, rnorm(n, mean = los_mean, sd = los_mean * 0.3)), 1)
    los_days[surgery == "No" & runif(n) < 0.3] <- NA
    
    # ============================================================================
    # RECOVERY DAYS
    # ============================================================================
    
    recovery_mean <- exp(2.5 +
                             (surgery == "Yes") * 0.5 +
                             0.01 * age +
                             (diabetes == "Yes") * 0.3 +
                             as.numeric(ecog) * 0.3 +
                             (any_complication == "Yes") * 0.4 +
                             (treatment == "Drug A") * (-0.2) +
                             (treatment == "Drug B") * 0.25 +
                             site_intercepts[as.character(site)] * 0.3)
    
    recovery_days <- round(pmax(1, rnorm(n, mean = recovery_mean, sd = recovery_mean * 0.25)))
    
    # ============================================================================
    # ADVERSE EVENT COUNT (overdispersed)
    # ============================================================================
    
    ae_lambda <- exp(0.5 +
                         0.015 * age +
                         as.numeric(ecog) * 0.3 +
                         (diabetes == "Yes") * 0.4 +
                         (stage == "II") * 0.2 +
                         (stage == "III") * 0.4 +
                         (stage == "IV") * 0.6 +
                         (treatment == "Drug A") * (-0.3) +
                         (treatment == "Drug B") * 0.5 +
                         (surgery == "Yes") * 0.3 +
                         site_intercepts[as.character(site)] * 0.4)
    
    theta <- 2
    ae_count <- rnbinom(n, size = theta, mu = ae_lambda)
    
    # ============================================================================
    # FOLLOW-UP VISIT COUNT (equidispersed)
    # ============================================================================
    
    fu_lambda <- exp(1.8 +
                         (stage == "II") * 0.15 +
                         (stage == "III") * 0.3 +
                         (stage == "IV") * 0.5 +
                         as.numeric(ecog) * 0.15 +
                         (treatment == "Drug A") * 0.1 +
                         (treatment == "Drug B") * 0.15 +
                         (surgery == "Yes") * 0.2)
    
    fu_count <- rpois(n, lambda = fu_lambda)
    
    # ============================================================================
    # OVERALL SURVIVAL
    # ============================================================================
    
    os_hazard <- exp(-1.5 +
                         site_age_slope[as.character(site)] * 0.04 * age +
                         (sex == "Male") * 0.25 +
                         (diabetes == "Yes") * 0.4 +
                         (smoking == "Current") * 0.5 +
                         (smoking == "Former") * 0.2 +
                         as.numeric(ecog) * 0.6 +
                         (grade == "Moderately differentiated") * 0.3 +
                         (grade == "Poorly differentiated") * 0.7 +
                         (stage == "II") * 0.4 +
                         (stage == "III") * 0.9 +
                         (stage == "IV") * 1.5 +
                         0.3 * log(biomarker_x) +
                         0.15 * log(biomarker_y) +
                         (surgery == "Yes") * (-0.6) +
                         (treatment == "Drug A") * (-0.5 + site_drugA_effect[as.character(site)]) +
                         (treatment == "Drug B") * (0.1 + site_drugB_effect[as.character(site)]) +
                         ((treatment == "Drug A") & (stage %in% c("III", "IV"))) * (-0.4) +
                         ((treatment == "Drug B") & (stage == "I")) * (-0.3) +
                         (age > 70 & diabetes == "Yes") * 0.6 +
                         (as.numeric(ecog) >= 3 & treatment == "Drug A") * 0.4 +
                         site_intercepts[as.character(site)] +
                         rnorm(n, 0, 0.3))
    
    os_months_raw <- round(rexp(n, rate = os_hazard) * 12, 1)
    
    censor_time <- runif(n, min = 12, max = 60)
    os_time <- pmin(os_months_raw, censor_time)
    os_status <- as.numeric(os_months_raw <= censor_time)
    
    # ============================================================================
    # PROGRESSION-FREE SURVIVAL
    # ============================================================================
    
    pfs_hazard <- exp(-1.2 +
                          site_age_slope[as.character(site)] * 0.035 * age +
                          (sex == "Male") * 0.2 +
                          (diabetes == "Yes") * 0.3 +
                          (smoking == "Current") * 0.4 +
                          as.numeric(ecog) * 0.5 +
                          (grade == "Moderately differentiated") * 0.35 +
                          (grade == "Poorly differentiated") * 0.8 +
                          (stage == "II") * 0.5 +
                          (stage == "III") * 1.0 +
                          (stage == "IV") * 1.6 +
                          0.12 * biomarker_x +
                          (surgery == "Yes") * (-0.5) +
                          (treatment == "Drug A") * (-0.4 + site_drugA_effect[as.character(site)]) +
                          (treatment == "Drug B") * (-0.1 + site_drugB_effect[as.character(site)]) +
                          ((treatment == "Drug A") & (stage %in% c("III", "IV"))) * (-0.3) +
                          ((treatment == "Drug B") & (stage == "I")) * (-0.2) +
                          (as.numeric(ecog) >= 3 & treatment == "Drug A") * 0.5 +
                          site_intercepts[as.character(site)] * 0.8 +
                          rnorm(n, 0, 0.3))
    
    pfs_months_raw <- round(rexp(n, rate = pfs_hazard) * 12, 1)
    pfs_censor_time <- censor_time
    pfs_time_initial <- pmin(pfs_months_raw, pfs_censor_time)
    pfs_status_initial <- as.numeric(pfs_months_raw <= pfs_censor_time)
    
    # Ensure PFS <= OS
    pfs_time <- ifelse(pfs_status_initial == 1 & pfs_time_initial > os_time,
                       os_time * runif(n, 0.5, 0.95),
                       pfs_time_initial)
    pfs_time <- pmin(pfs_time, os_time)
    pfs_status <- ifelse(pfs_time == os_time & os_status == 1 & pfs_status_initial == 0,
                         0,
                         pfs_status_initial)
    pfs_time <- round(pfs_time, 1)
    
    # ============================================================================
    # MISSING DATA
    # ============================================================================
    
    missing_idx <- sample(n, size = n * 0.02)
    bmi[sample(missing_idx, 12)] <- NA
    smoking[sample(missing_idx, 17)] <- NA
    hypertension[sample(missing_idx, 15)] <- NA
    diabetes[sample(missing_idx, 16)] <- NA
    ecog[sample(missing_idx, 8)] <- NA
    creatinine[sample(missing_idx, 10)] <- NA
    hemoglobin[sample(missing_idx, 10)] <- NA
    biomarker_x[sample(missing_idx, 8)] <- NA
    biomarker_y[sample(missing_idx, 8)] <- NA
    grade[sample(missing_idx, 10)] <- NA
    stage[sample(missing_idx, 3)] <- NA
    
    # ============================================================================
    # ASSEMBLE DATASET
    # ============================================================================
    
    trial_data <- data.frame(
        patient_id = sprintf("PT%04d", 1:n),
        site = site,
        age = age,
        sex = sex,
        race = race,
        ethnicity = ethnicity,
        bmi = bmi,
        smoking = smoking,
        hypertension = hypertension,
        diabetes = diabetes,
        ecog = ecog,
        creatinine = creatinine,
        hemoglobin = hemoglobin,
        biomarker_x = biomarker_x,
        biomarker_y = biomarker_y,
        grade = grade,
        stage = stage,
        treatment = treatment,
        surgery = surgery,
        icu_admission = icu_admission,
        pain_score = pain_score,
        wound_infection = wound_infection,
        any_complication = any_complication,
        readmission_30d = readmission_30d,
        los_days = los_days,
        recovery_days = recovery_days,
        ae_count = ae_count,
        fu_count = fu_count,
        pfs_months = pfs_time,
        pfs_status = pfs_status,
        os_months = round(os_time, 1),
        os_status = os_status,
        stringsAsFactors = FALSE
    )
    
    return(trial_data)
}

clintrial <- create_clintrial()

#' Variable Labels for Clinical Trial Dataset
#'
#' A named character vector providing human-readable labels for all variables
#' in the \code{\link{clintrial}} dataset. These labels are designed for use
#' in publication-ready tables and forest plots.
#'
#' @format A named character vector with 32 elements. Names correspond to 
#'   variable names in \code{clintrial}, values are display labels.
#'
#' @details
#' Labels follow standard clinical reporting conventions:
#' \itemize{
#'   \item Units are included in parentheses where applicable
#'   \item Survival outcomes include both time and status labels
#'   \item Composite labels for Surv() objects are provided
#' }
#'
#' @seealso \code{\link{clintrial}} for the dataset these labels describe
#'
#' @examples
#' data(clintrial_labels)
#' 
#' # View all labels
#' clintrial_labels
#' 
#' # Use in forest plots
#' model <- glm(os_status ~ age + sex + treatment, 
#'              data = clintrial, family = binomial)
#' glmforest(model, data = clintrial, labels = clintrial_labels)
#' 
#' # Use in fit() tables
#' fit(data = clintrial, 
#'     outcome = "os_status", 
#'     predictors = c("age", "sex", "treatment"),
#'     labels = clintrial_labels)
#'
#' @name clintrial_labels
#' @docType data
#' @keywords datasets
NULL

clintrial_labels <- c(
    "patient_id" = "Patient ID",
    "site" = "Study Site",
    "age" = "Age (years)",
    "sex" = "Sex",
    "race" = "Race",
    "ethnicity" = "Ethnicity",
    "bmi" = "Body Mass Index (kg/m²)",
    "smoking" = "Smoking Status",
    "hypertension" = "Hypertension",
    "diabetes" = "Diabetes",
    "ecog" = "ECOG Performance Status",
    "creatinine" = "Baseline Creatinine (mg/dL)",
    "hemoglobin" = "Baseline Hemoglobin (g/dL)",
    "biomarker_x" = "Biomarker X (ng/mL)",
    "biomarker_y" = "Biomarker Y (U/L)",
    "grade" = "Tumor Grade",
    "stage" = "Disease Stage",
    "treatment" = "Treatment Group",
    "surgery" = "Surgical Resection",
    "icu_admission" = "ICU Admission Required",
    "pain_score" = "Postoperative Pain Score (0-10)",
    "wound_infection" = "Wound/Site Infection",
    "any_complication" = "Any Complication",
    "readmission_30d" = "30-Day Readmission",
    "los_days" = "Length of Hospital Stay (days)",
    "recovery_days" = "Days to Functional Recovery",
    "ae_count" = "Adverse Event Count",
    "fu_count" = "Follow-Up Visit Count",
    "pfs_months" = "Progression-Free Survival Time (months)",
    "pfs_status" = "Progression or Death Event",
    "os_months" = "Overall Survival Time (months)",
    "os_status" = "Death Event",
    "Surv(pfs_months, pfs_status)" = "Progression-Free Survival (months)",
    "Surv(os_months, os_status)" = "Overall Survival (months)"
)

usethis::use_data(clintrial, overwrite = TRUE)
usethis::use_data(clintrial_labels, overwrite = TRUE)
