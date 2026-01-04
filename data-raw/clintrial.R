#' Simulated Clinical Trial Dataset
#'
#' A simulated dataset representing a multi-site oncology clinical trial with 
#' 850 patients. The dataset includes realistic correlations, interaction effects,
#' site-level clustering, and a comprehensive set of outcomes suitable for 
#' demonstrating regression analysis, survival modeling, and multi-outcome analyses.
#'
#' @format A data frame with 850 rows and 30 variables:
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
#'   \item{postop_complication}{Factor. Any complication occurred (No, Yes)}
#'   \item{wound_infection}{Factor. Wound or site infection (No, Yes)}
#'   \item{readmission_30d}{Factor. 30-day hospital readmission (No, Yes)}
#'   \item{icu_admission}{Factor. ICU admission required (No, Yes)}
#'   
#'   \strong{Recovery Outcomes:}
#'   \item{pain_score}{Numeric. Postoperative pain score on 0-10 visual analog scale}
#'   \item{recovery_days}{Numeric. Days to functional recovery}
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
#'     faster recovery, improved survival)
#'   \item Drug B shows increased toxicity (more complications, longer recovery)
#'   \item Treatment effects vary by disease stage (interaction effects)
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
#' \itemize
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

set.seed(71)  # For reproducibility
n <- 850

create_clintrial <- function() {

    ## ------------------------
    ## Baseline characteristics
    ## ------------------------
    
    ## Site (generate early for clustering effects)
    greek_letters <- c("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Iota", "Kappa")
    site <- factor(
        paste0("Site ", sample(greek_letters, n, replace = TRUE)),
        levels = paste0("Site ", greek_letters))
    
    ## Generate site-specific random effects
    n_sites <- length(unique(site))
    site_names <- levels(site)
    
    ## Site random intercepts (some sites have better/worse outcomes)
    site_intercepts <- rnorm(n_sites, mean = 0, sd = 0.5)
    names(site_intercepts) <- site_names
    
    ## Site-specific treatment effects (some sites better with Drug A, others with Drug B)
    site_drugA_effect <- rnorm(n_sites, mean = 0, sd = 0.3)
    site_drugB_effect <- rnorm(n_sites, mean = 0, sd = 0.3)
    names(site_drugA_effect) <- site_names
    names(site_drugB_effect) <- site_names
    
    ## Site-specific age effects (age matters more at some sites)
    site_age_slope <- rnorm(n_sites, mean = 1, sd = 0.2)
    names(site_age_slope) <- site_names
    
    ## Basic demographics
    age <- round(rnorm(n, mean = 60, sd = 12))
    age[age < 18] <- 18
    age[age > 90] <- 90
    
    sex <- factor(sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.55, 0.45)))
    
    race <- factor(sample(c("White", "Black", "Asian", "Other"), n, 
                          replace = TRUE, prob = c(0.70, 0.15, 0.10, 0.05)),
                   levels = c("White", "Black", "Asian", "Other"))
    
    ethnicity <- factor(sample(c("Non-Hispanic", "Hispanic"), n, 
                               replace = TRUE, prob = c(0.88, 0.12)),
                        levels = c("Non-Hispanic", "Hispanic"))
    
    ## Clinical characteristics
    bmi <- round(rnorm(n, mean = 28, sd = 5), 1)
    bmi[bmi < 15] <- 15
    bmi[bmi > 50] <- 50
    
    smoking <- factor(sample(c("Never", "Former", "Current"), n,
                             replace = TRUE, prob = c(0.45, 0.35, 0.20)),
                      levels = c("Never", "Former", "Current"))
    
    hypertension <- factor(sample(c("No", "Yes"), n, replace = TRUE,
                                  prob = c(0.60, 0.40)))
    
    diabetes <- factor(sample(c("No", "Yes"), n, replace = TRUE, 
                              prob = c(0.75, 0.25)))
    
    ## Performance status (generate early for use in other variables)
    ecog <- factor(sample(0:3, n, replace = TRUE, 
                          prob = c(0.30, 0.40, 0.25, 0.05)),
                   levels = 0:3)
    
    ## Disease characteristics
    grade <- factor(sample(c("Well-differentiated", "Moderately differentiated", 
                             "Poorly differentiated"), n,
                           replace = TRUE, prob = c(0.20, 0.50, 0.30)),
                    levels = c("Well-differentiated", "Moderately differentiated", 
                               "Poorly differentiated"))

    stage <- factor(sample(c("I", "II", "III", "IV"), n,
                           replace = TRUE, prob = c(0.25, 0.30, 0.30, 0.15)),
                    levels = c("I", "II", "III", "IV"))
    
    ## Biomarkers (continuous with some correlation to outcome)
    biomarker_x <- round(rgamma(n, shape = 2, rate = 0.5) + 
                         as.numeric(stage) * 0.5 + 
                         (age > 65) * 2 +
                         rnorm(n, 0, 0.5), 2)
    
    biomarker_y <- round(rlnorm(n, meanlog = 3, sdlog = 1) + 
                         (sex == "Male") * 20 + rnorm(n, 0, 10), 1)
    biomarker_y[biomarker_y < 0] <- 0

    ## Laboratory values
    creatinine <- round(rlnorm(n, meanlog = 0, sdlog = 0.3), 2)
    creatinine[creatinine > 5] <- 5

    hemoglobin <- round(rnorm(n, mean = 13, sd = 2), 1)
    hemoglobin[hemoglobin < 7] <- 7
    hemoglobin[hemoglobin > 18] <- 18

    ## -------------
    ## Interventions
    ## -------------
    
    ## Treatment assignment with bias for Drug B
    treatment <- character(n)
    for (i in 1:n) {
        ## Drug B biased toward older, sicker patients (max additional 0.2)
        drug_b_bias <- 0
        drug_b_bias <- drug_b_bias + (age[i] > 65) * 0.10
        drug_b_bias <- drug_b_bias + (as.numeric(ecog[i]) >= 3) * 0.10
        drug_b_bias <- drug_b_bias + (grade[i] == "Poorly differentiated") * 0.05
        drug_b_bias <- drug_b_bias + (stage[i] %in% c("III", "IV")) * 0.10
        
        ## Drug A has slight opposite bias (max additional 0.1)
        drug_a_bias <- 0
        drug_a_bias <- drug_a_bias + (age[i] <= 55) * 0.05
        drug_a_bias <- drug_a_bias + (as.numeric(ecog[i]) <= 2) * 0.05
        
        ## Base probabilities
        base_prob <- 1/3
        
        ## Calculate final probabilities (ensuring they sum to 1)
        probs <- c(
            Control = base_prob - drug_a_bias/2 - drug_b_bias/2,
            DrugA = base_prob + drug_a_bias,
            DrugB = base_prob + drug_b_bias
        )
        
        ## Ensure no negative probabilities
        probs[probs < 0.05] <- 0.05  # Minimum 5% chance for each
        probs <- probs / sum(probs)  # Normalize to sum to 1
        
        treatment[i] <- sample(c("Control", "Drug A", "Drug B"), 1, prob = probs)
    }
    treatment <- factor(treatment, levels = c("Control", "Drug A", "Drug B"))
    
    ## Surgery (biased toward younger, healthier patients)
    surgery_prob <- plogis(-1 + 
                           (65 - age) * 0.05 +  # Younger more likely
                           (ecog == "0") * 1.2 +  # ECOG 0 much more likely
                           (ecog == "1") * 0.6 +  # ECOG 1 somewhat more likely
                           (ecog == "2") * (-0.5) +  # ECOG 2 less likely
                                         (ecog == "3") * (-2) +  # ECOG 3 very unlikely
                                                       (treatment == "Drug A") * 0.8 +  # Drug A more likely
                                                       (treatment == "Drug B") * (-0.8) +  # Drug B less likely
                                                                               (stage == "I") * 1 +  # Stage I more likely
                                                                               (stage == "II") * 0.5 +  #  Stage II slightly more likely
                                                                               (stage == "IV") * (-1.5) +  # Stage IV stage less likely
                                                                                               rnorm(n, 0, 0.3))
    
    surgery <- factor(ifelse(runif(n) < surgery_prob, "Yes", "No"),
                      levels = c("No", "Yes"))

    ## -----------------------
    ## Post-treatment outcomes
    ## -----------------------
    
    ## Hospital length of stay in days (influenced by multiple factors with interactions)
    los_days <- 5 +  # Baseline stay
        ## Demographics
        (age - 60) * 0.1 +  # Older patients stay longer
        (sex == "Male") * 0.5 +  # Males slightly longer stay
        
        ## Clinical factors (STRONG effects)
        as.numeric(ecog) * 2.5 +  # Performance status major factor
        (smoking == "Current") * 1.2 +
        (smoking == "Former") * 0.5 +
        (diabetes == "Yes") * 1.8 +  # Diabetes increases LOS
        (hypertension == "Yes") * 0.8 +
        
        ## Disease characteristics
        (stage == "II") * 1.5 +
        (stage == "III") * 3.0 +
        (stage == "IV") * 4.5 +
        (grade == "Moderately differentiated") * 1.0 +
        (grade == "Poorly differentiated") * 2.0 +
        
        ## Treatment effects
        (surgery == "Yes") * 4.5 +  # Surgery increases LOS significantly
        (treatment == "Drug A") * (-0.8) +  # Drug A reduces LOS slightly
                                    (treatment == "Drug B") * 1.2 +  # Drug B increases LOS (side effects)
                                    
                                    ## INTERACTION EFFECTS
                                        # Treatment × Stage interactions
                                    (treatment == "Drug A" & stage %in% c("III", "IV")) * (-2.0) +  # Drug A particularly effective in late stage
                                                                                            (treatment == "Drug B" & stage %in% c("III", "IV")) * 1.5 +  # Drug B less effective in late stage
                                                                                            
                                        # Age × Diabetes interaction
                                            ((age > 70) & (diabetes == "Yes")) * 3.0 +  # Elderly diabetics stay much longer
                                            
                                        # Sex × Treatment interaction
                                        ((sex == "Male") & (treatment == "Drug B")) * 1.8 +  # Males have more side effects with Drug B
                                        
                                        ## Lab values
                                        (creatinine - 1) * 1.5 +  # Higher creatinine = longer stay
                                        (13 - hemoglobin) * 0.4 +  # Lower hemoglobin = longer stay
                                        
                                        ## Biomarkers with age interaction
                                        biomarker_x * 0.3 +
                                        ((age > 65) * biomarker_x * 0.2) +  # Biomarker X matters more in elderly
                                        
                                        ## Site clustering effect
                                        site_intercepts[as.character(site)] * 2 +  # Site random effect
                                        
                                        ## Random variation
                                        rnorm(n, 0, 2)
    
    ## Clean up LOS values
    los_days <- round(los_days, 1)
    los_days[los_days < 1] <- 1  # Minimum 1 day
    los_days[los_days > 60] <- 60  # Cap at 60 days
    
    ## Add some missing values for LOS (e.g., outpatient procedures)
    los_days[sample(which(surgery == "No"), size = 20)] <- NA

    ## Generate site-specific complication rates (some sites have better/worse outcomes)
    site_complication_effect <- rnorm(n_sites, mean = 0, sd = 0.4)
    names(site_complication_effect) <- site_names
    
    ## Any Postoperative Complication (composite)
    ## Strong associations with: age, ECOG, diabetes, surgery, treatment
    complication_prob <- plogis(-2.0 +
                                ## Demographics
                                (age - 60) * 0.025 +  # Older = higher risk
                                (sex == "Male") * 0.2 +
                                
                                ## Comorbidities (STRONG effects)
                                (diabetes == "Yes") * 0.6 +  # Diabetes major risk factor
                                (hypertension == "Yes") * 0.3 +
                                (smoking == "Current") * 0.5 +
                                (smoking == "Former") * 0.2 +
                                as.numeric(ecog) * 0.4 +  # Poor performance status = higher risk
                                
                                ## Disease characteristics
                                (stage == "III") * 0.3 +
                                (stage == "IV") * 0.6 +
                                
                                ## Treatment effects
                                (surgery == "Yes") * 1.2 +  # Surgery is major risk factor
                                (treatment == "Drug A") * (-0.3) +  # Drug A protective
                                                        (treatment == "Drug B") * 0.4 +  # Drug B increases complications
                                                        
                                                        ## Lab values
                                                        (creatinine - 1) * 0.5 +  # Renal function
                                                        (13 - hemoglobin) * 0.15 +  # Anemia
                                                        
                                                        ## INTERACTIONS
                                                        ((age > 70) & (diabetes == "Yes")) * 0.5 +  # Elderly diabetics highest risk
                                                        ((surgery == "Yes") & (treatment == "Drug B")) * 0.4 +  # Drug B + surgery = more complications
                                                        ((smoking == "Current") & (surgery == "Yes")) * 0.3 +  # Smokers + surgery
                                                        
                                                        ## Site clustering
                                                        site_complication_effect[as.character(site)] +
                                                        
                                                        rnorm(n, 0, 0.2))
    
    any_complication <- factor(ifelse(runif(n) < complication_prob, "Any complication", "No complication"),
                               levels = c("No complication", "Any complication"))
    ## Note: All patients can have complications (surgical or treatment-related)
    ## Surgery is already accounted for in the probability model
    
    ## Wound Infection (SSI)
    ## Subset of complications, more specific risk factors
    wound_prob <- plogis(-3.5 +
                         ## Demographics  
                         (age - 60) * 0.015 +
                         
                         ## Strong risk factors for SSI
                         (diabetes == "Yes") * 0.8 +  # Diabetes = major SSI risk
                         (bmi - 28) * 0.05 +  # Obesity increases risk
                         (smoking == "Current") * 0.7 +  # Smoking impairs healing
                         (smoking == "Former") * 0.2 +
                         
                         ## Surgery-related
                         (surgery == "Yes") * 2.0 +  # Only surgical patients get SSI
                         
                         ## Treatment
                         (treatment == "Drug B") * 0.4 +  # Drug B immunosuppressive
                         
                         ## Lab values
                         (hemoglobin < 10) * 0.5 +  # Severe anemia
                         
                         ## INTERACTIONS
                         ((diabetes == "Yes") & (bmi > 35)) * 0.6 +  # Obese diabetics
                         ((smoking == "Current") & (diabetes == "Yes")) * 0.4 +
                         
                         ## Site clustering (some sites have better sterile technique)
                         site_complication_effect[as.character(site)] * 0.8 +
                         
                         rnorm(n, 0, 0.15))
    
    wound_infection <- factor(ifelse(runif(n) < wound_prob, "Wound infection", "No wound infection"),
                              levels = c("No wound infection", "Wound infection"))
    ## Note: Non-surgical patients have very low SSI risk (from procedures, lines, etc.)
    ## Surgery effect is captured in the probability model (+2.0 log-odds)
    
    ## 30-Day Readmission
    ## Related to complications but also social/system factors
    ## NOTE: Strengthened effects to ensure detectable associations in regression
    readmit_prob <- plogis(-2.8 +
                           ## Demographics (stronger effects)
                           (age - 60) * 0.035 +  # Age effect strengthened
                           (sex == "Male") * 0.35 +  # Male sex increases risk
                           
                           ## Clinical factors (stronger effects)
                           (diabetes == "Yes") * 0.7 +  # Diabetes strongly predicts readmission
                           (hypertension == "Yes") * 0.25 +
                           as.numeric(ecog) * 0.45 +  # ECOG performance status
                           (smoking == "Current") * 0.5 +  # Current smoking
                           (smoking == "Former") * 0.15 +
                           
                           ## Prior complications increase readmission
                           (any_complication == "Yes") * 1.2 +  # Strong predictor
                           (wound_infection == "Yes") * 0.9 +
                           
                           ## Disease severity (stronger stage effects)
                           (stage == "II") * 0.2 +
                           (stage == "III") * 0.5 +
                           (stage == "IV") * 0.85 +
                           
                           ## Treatment effects (more differentiated)
                           (treatment == "Drug A") * (-0.4) +  # Drug A patients do better
                                                   (treatment == "Drug B") * 0.5 +  # Drug B increases readmission
                                                   
                                                   ## Surgery effect
                                                   (surgery == "Yes") * 0.6 +  # Surgical patients more likely to be readmitted
                                                   
                                                   ## Lab values
                                                   (creatinine - 1) * 0.4 +  # Renal function matters
                                                   (hemoglobin < 10) * 0.35 +  # Anemia
                                                   
                                                   ## LOS effect (very short or very long LOS = readmission risk)
                                                   ifelse(!is.na(los_days), (abs(los_days - 7) * 0.04), 0) +
                                                   
                                                   ## INTERACTION EFFECTS
                                                   ((diabetes == "Yes") & (age > 70)) * 0.4 +  # Elderly diabetics
                                                   ((surgery == "Yes") & (stage %in% c("III", "IV"))) * 0.35 +
                                                   
                                                   ## Site clustering (some sites have better discharge planning)
                                                   site_complication_effect[as.character(site)] * 0.7 +
                                                   
                                                   rnorm(n, 0, 0.15))
    
    readmission_30d <- factor(ifelse(runif(n) < readmit_prob, "30-day readmission", "No 30-day readmission"),
                              levels = c("No 30-day readmission", "30-day readmission"))
    
    ## ICU Admission
    ## More severe outcome, stronger associations
    icu_prob <- plogis(-3.0 +
                       ## Demographics
                       (age - 60) * 0.03 +  # Age is important
                       
                       ## Clinical factors (STRONG)
                       as.numeric(ecog) * 0.5 +
                       (diabetes == "Yes") * 0.4 +
                       (hypertension == "Yes") * 0.3 +
                       
                       ## Disease severity
                       (stage == "III") * 0.5 +
                       (stage == "IV") * 1.0 +
                       (grade == "Poorly differentiated") * 0.3 +
                       
                       ## Surgery major factor
                       (surgery == "Yes") * 1.5 +
                       
                       ## Treatment
                       (treatment == "Drug B") * 0.5 +  # Drug B more toxic
                       
                       ## Lab values
                       (creatinine - 1) * 0.8 +  # Renal function critical
                       (hemoglobin < 9) * 0.6 +  # Severe anemia
                       
                       ## INTERACTIONS
                       ((age > 75) & (surgery == "Yes")) * 0.6 +  # Elderly surgical patients
                       ((stage == "IV") & (surgery == "Yes")) * 0.5 +
                       
                       ## Site clustering
                       site_complication_effect[as.character(site)] * 0.5 +
                       
                       rnorm(n, 0, 0.15))
    
    icu_admission <- factor(ifelse(runif(n) < icu_prob, "ICU admission", "No ICU admission"),
                            levels = c("No ICU admission", "ICU admission"))
    
    ## Postoperative Pain Score (0-10)
    ## Visual analog scale, varies by patient and treatment
    pain_score <- 4.0 +  # Baseline moderate pain
        ## Demographics
        (age - 60) * (-0.02) +  # Older patients report less pain (counterintuitive but observed)
                       (sex == "Female") * 0.5 +  # Women report more pain
                       
                       ## Clinical factors
                       as.numeric(ecog) * 0.3 +
                       (smoking == "Current") * 0.4 +
                       
                       ## Disease/procedure factors
                       (surgery == "Yes") * 2.0 +  # Surgery = more pain
                       (stage == "III") * 0.3 +
                       (stage == "IV") * 0.5 +
                       
                       ## Treatment effects
                       (treatment == "Drug A") * (-0.8) +  # Drug A has analgesic properties
                                                   (treatment == "Drug B") * 0.5 +  # Drug B causes more pain
                                                   
                                                   ## Complications increase pain
                                                   ifelse(!is.na(any_complication) & any_complication == "Yes", 1.5, 0) +
                                                   ifelse(!is.na(wound_infection) & wound_infection == "Yes", 1.0, 0) +
                                                   
                                                   ## INTERACTIONS
                                                   ((sex == "Female") & (treatment == "Drug B")) * 0.4 +  # Women + Drug B
                                                   ((surgery == "Yes") & (stage == "IV")) * 0.5 +
                                                   
                                                   ## Site clustering (some sites have better pain management)
                                                   site_intercepts[as.character(site)] * (-0.5) +
                                                                                           
                                                                                           rnorm(n, 0, 1.5)
    
    ## Clean up pain scores
    pain_score <- round(pain_score, 1)
    pain_score[pain_score < 0] <- 0
    pain_score[pain_score > 10] <- 10
    
    ## Days to Functional Recovery
    ## Time to return to baseline function, log-normal distribution
    recovery_days <- exp(2.0 +  # Baseline ~7 days
                         ## Demographics
                         (age - 60) * 0.015 +  # Older = slower recovery
                         
                         ## Clinical factors
                         as.numeric(ecog) * 0.15 +  # Baseline function matters
                         (diabetes == "Yes") * 0.2 +
                         (smoking == "Current") * 0.15 +
                         
                         ## Disease/procedure factors
                         (surgery == "Yes") * 0.5 +  # Surgery extends recovery
                         (stage == "III") * 0.1 +
                         (stage == "IV") * 0.2 +
                         
                         ## Treatment effects
                         (treatment == "Drug A") * (-0.15) +  # Drug A faster recovery
                                                 (treatment == "Drug B") * 0.2 +  # Drug B slower recovery
                                                 
                                                 ## Complications extend recovery significantly
                                                 ifelse(!is.na(any_complication) & any_complication == "Yes", 0.4, 0) +
                                                 ifelse(!is.na(wound_infection) & wound_infection == "Yes", 0.5, 0) +
                                                 (icu_admission == "Yes") * 0.6 +
                                                 
                                                 ## INTERACTIONS
                                                 ((age > 70) & (surgery == "Yes")) * 0.2 +  # Elderly surgical patients
                                                 ((diabetes == "Yes") & (wound_infection == "Yes")) * 0.3 +
                                                 
                                                 ## Site clustering
                                                 site_intercepts[as.character(site)] * 0.2 +
                                                 
                                                 rnorm(n, 0, 0.3))
    
    ## Clean up recovery days
    recovery_days <- round(recovery_days, 0)
    recovery_days[recovery_days < 1] <- 1
    recovery_days[recovery_days > 90] <- 90
    
    ## Add some missing values for recovery (lost to follow-up)
    recovery_days[sample(n, size = 15)] <- NA
    pain_score[sample(n, size = 10)] <- NA

    ## -----------------
    ## Survival outcomes
    ## -----------------
    
    ## Survival outcomes (correlated with predictors including interactions)
    hazard <- exp(-3 + 
                  ## Disease characteristics (STRONGER effects)
                  (stage == "II") * 0.4 +
                  (stage == "III") * 0.8 +  
                  (stage == "IV") * 1.5 +   
                  (grade == "Moderately differentiated") * 0.3 +
                  (grade == "Poorly differentiated") * 0.6 +  
                  
                  ## Demographics and clinical
                  (age - 60) * 0.02 * site_age_slope[as.character(site)] +  # Age effect varies by site
                  as.numeric(ecog) * 0.3 +
                  (sex == "Male") * 0.15 +
                  (smoking == "Current") * 0.25 +
                  
                  ## Biomarkers
                  biomarker_x * 0.1 +
                  
                  ## Treatments (ENSURE Drug A benefit)
                  (surgery == "Yes") * (-0.6) +
                                     (treatment == "Drug A") * (-0.5 + site_drugA_effect[as.character(site)]) +  # Drug A effect varies by site
                                                             (treatment == "Drug B") * (-0.15 + site_drugB_effect[as.character(site)]) +  # Drug B effect varies by site
                                                                                     
                                                                                     ## INTERACTION EFFECTS for survival
                                        # Treatment × Stage interaction (Drug A works better in advanced stages)
                                        ((treatment == "Drug A") & (stage %in% c("III", "IV"))) * (-0.4) +
                                                                                                ((treatment == "Drug B") & (stage %in% c("III", "IV"))) * 0.3 +
                                                                                                
                                        # Sex × Treatment interaction
                                        ((sex == "Male") & (treatment == "Drug A")) * 0.2 +  # Drug A less effective in males
                                        ((sex == "Female") & (treatment == "Drug B")) * (-0.3) +  # Drug B more effective in females
                                                                                      
                                        # Age × Biomarker interaction
                                        ((age > 65) * (biomarker_x > median(biomarker_x, na.rm = TRUE))) * 0.5 +  # High biomarker worse in elderly
                                        
                                        # Smoking × Diabetes interaction
                                        ((smoking == "Current") & (diabetes == "Yes")) * 0.4 +  # Compound risk
                                        
                                        ## Site clustering effect
                                        site_intercepts[as.character(site)] +  # Site random intercept
                                        
                                        rnorm(n, 0, 0.3))
    
    ## Generate survival times
    surv_months <- round(rexp(n, rate = hazard) * 12, 1)
    surv_months[surv_months > 120] <- 120  # Cap at 10 years
    
    ## Generate censoring (administrative at 60 months + random)
    censor_time <- pmin(60, rexp(n, rate = 0.02) * 12)
    
    ## Observed time and status
    os_time <- pmin(surv_months, censor_time)
    os_status <- as.numeric(surv_months <= censor_time)
    
    ## Progression-free survival (PFS must be <= OS)
    ## PFS hazard is higher than OS hazard (progression happens before death)
    pfs_hazard <- exp(-2.5 +  # Higher baseline hazard than OS
                      ## Disease characteristics (similar effects as OS)
                      (stage == "II") * 0.45 +
                      (stage == "III") * 0.9 +
                      (stage == "IV") * 1.7 +
                      (grade == "Moderately differentiated") * 0.35 +
                      (grade == "Poorly differentiated") * 0.7 +
                      
                      ## Demographics and clinical
                      (age - 60) * 0.015 * site_age_slope[as.character(site)] +  # Site-specific age effect
                      as.numeric(ecog) * 0.35 +
                      (sex == "Male") * 0.1 +
                      (smoking == "Current") * 0.2 +
                      
                      ## Biomarkers
                      biomarker_x * 0.12 +
                      
                      ## Treatments (Drug A shows benefit on PFS)
                      (surgery == "Yes") * (-0.5) +
                                         (treatment == "Drug A") * (-0.4 + site_drugA_effect[as.character(site)]) +  # Site-specific effect
                                                                 (treatment == "Drug B") * (-0.1 + site_drugB_effect[as.character(site)]) +  # Site-specific effect
                                                                                         
                                                                                         ## INTERACTION EFFECTS for PFS
                                        # Treatment × Stage interaction
                                        ((treatment == "Drug A") & (stage %in% c("III", "IV"))) * (-0.3) +
                                                                                                ((treatment == "Drug B") & (stage == "I")) * (-0.2) +  # Drug B works better in early stage
                                                                                                                                           
                                        # ECOG × Treatment interaction
                                        (as.numeric(ecog) >= 3 & treatment == "Drug A") * 0.5 +  # Drug A less effective in poor PS
                                        
                                        ## Site clustering effect
                                        site_intercepts[as.character(site)] * 0.8 +  # Site effect (slightly dampened vs OS)
                                        
                                        rnorm(n, 0, 0.3))
    
    ## Generate progression times
    pfs_months_raw <- round(rexp(n, rate = pfs_hazard) * 12, 1)
    
    ## Generate PFS censoring (same as OS censoring)
    pfs_censor_time <- censor_time
    
    ## Observed PFS time and status
    pfs_time_initial <- pmin(pfs_months_raw, pfs_censor_time)
    pfs_status_initial <- as.numeric(pfs_months_raw <= pfs_censor_time)
    
    ## Ensure PFS time is always <= OS time
    ## For patients with progression events, PFS must be before OS
    pfs_time <- ifelse(pfs_status_initial == 1 & pfs_time_initial > os_time,
                       os_time * runif(n, 0.5, 0.95),  # Set PFS to 50-95% of OS time
                       pfs_time_initial)
    
    ## Also ensure no PFS time exceeds OS time even for censored observations
    pfs_time <- pmin(pfs_time, os_time)
    
    ## Update PFS status: if OS event occurred and PFS time equals OS time, patient died without documented progression
    pfs_status <- ifelse(pfs_time == os_time & os_status == 1 & pfs_status_initial == 0,
                         0,  # Died without progression
                         pfs_status_initial)
    
    ## Round final PFS time
    pfs_time <- round(pfs_time, 1)
    
    ## * Add some missing values realistically
    missing_idx <- sample(n, size = n * 0.02)  # 2% missing
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
    
    ## * Create dataset
    trial_data <- data.frame(
        ## Identifiers
        patient_id = sprintf("PT%04d", 1:n),
        site = site,
        ## Baseline characteristics
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
        ## Interventions
        treatment = treatment,
        surgery = surgery,
        ## Post-treatment outcomes
        icu_admission = icu_admission,
        pain_score = pain_score,
        wound_infection = wound_infection,
        any_complication = any_complication,
        readmission_30d = readmission_30d,
        los_days = los_days,
        recovery_days = recovery_days,
        ## Survival outcomes
        pfs_months = pfs_time,
        pfs_status = pfs_status,
        os_months = round(os_time, 1),
        os_status = os_status,
        stringsAsFactors = FALSE
    )
    
    return(trial_data)
}

## Create the dataset
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

## Create variable labels
clintrial_labels <- c(
    ## Identifiers
    "patient_id" = "Patient ID",
    "site" = "Study Site",
    ## Baseline characteristics
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
    ## Interventions
    "treatment" = "Treatment Group",
    "surgery" = "Surgical Resection",
    ## Post-treatment outcomes
    "icu_admission" = "ICU Admission Required",
    "pain_score" = "Postoperative Pain Score (0-10)",
    "wound_infection" = "Wound/Site Infection",
    "any_complication" = "Any Complication",
    "readmission_30d" = "30-Day Readmission",
    "los_days" = "Length of Hospital Stay (days)",
    "recovery_days" = "Days to Functional Recovery",
    ## Survival outcomes
    "pfs_months" = "Progression-Free Survival Time (months)",
    "pfs_status" = "Progression or Death Event",
    "os_months" = "Overall Survival Time (months)",
    "os_status" = "Death Event",
    "Surv(pfs_months, pfs_status)" = "Progression-Free Survival (months)",
    "Surv(os_months, os_status)" = "Overall Survival (months)"
)

## Save both objects to "data/" directory
usethis::use_data(clintrial, overwrite = TRUE)
usethis::use_data(clintrial_labels, overwrite = TRUE)
