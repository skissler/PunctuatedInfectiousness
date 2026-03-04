# ==============================================================================
# Master analysis script — orchestrates all analyses
#
# Run from the project root:
#   Rscript code/run_analysis.R
#
# Scripts are ordered to match the writeup (writeup/punctuatedinfectiousness.md):
#
#   0. Setup              — ensure output directories exist
#   1. Validation         — test_equivalence: verify sim_stochastic_fast and
#                           factory profiles match reference implementations
#   2. simulate_epidemics — uncontrolled epidemic dynamics under different
#                           individual infectiousness profiles (Sections 3–4)
#   3. early_growth_analysis — first-transmission times, burstiness, and
#                           early growth rates (Sections 3–4 continued)
#   4. visualize_profiles_gamma — Gamma convolutional model visualization
#                           (Section 5)
#   5. test_based_screening — impact of punctuated infectiousness on
#                           test-based screening (Section 6)
#   6. analyze_overdispersion — biological/contact decomposition and
#                           overdispersion from periodic contacts (Sections 7–8)
#   7. simulate_trials    — identifiability of kappa from household/contact-
#                           tracing data (Section 9)
# ==============================================================================

cat("==============================================================\n")
cat("  Punctuated Infectiousness — Full Analysis Pipeline\n")
cat("==============================================================\n\n")

# ==============================================================================
# 0. Setup
# ==============================================================================

cat("--- 0. Setup ---\n")
if (!dir.exists("figures")) dir.create("figures")

# ==============================================================================
# 1. Validation
# ==============================================================================

cat("\n--- 1. Validating simulation infrastructure ---\n\n")
source("code/test_equivalence.R")

# ==============================================================================
# 2. Uncontrolled epidemic simulations (Sections 3–4)
# ==============================================================================

cat("\n--- 2. Simulating uncontrolled epidemics ---\n\n")
source("code/simulate_epidemics.R")

# ==============================================================================
# 3. Early growth analysis (Sections 3–4 continued)
# ==============================================================================

cat("\n--- 3. Early growth analysis ---\n\n")
source("code/early_growth_analysis.R")

# ==============================================================================
# 4. Gamma profile visualization (Section 5)
# ==============================================================================

cat("\n--- 4. Visualizing Gamma infectiousness profiles ---\n\n")
source("code/visualize_profiles_gamma.R")

# ==============================================================================
# 5. Test-based screening (Section 6)
# ==============================================================================

cat("\n--- 5. Test-based screening analysis ---\n\n")
source("code/test_based_screening.R")

# ==============================================================================
# 6. Overdispersion from periodic contacts (Sections 7–8)
# ==============================================================================

cat("\n--- 6. Overdispersion analysis ---\n\n")
source("code/analyze_overdispersion.R")

# ==============================================================================
# 7. Identifiability / trial simulations (Section 9)
# ==============================================================================

cat("\n--- 7. Trial simulations for identifiability ---\n\n")
source("code/simulate_trials.R")

# ==============================================================================
# 8. Gathering size restrictions (Section 10)
# ==============================================================================

cat("\n--- 8. Gathering size restrictions analysis ---\n\n")
source("code/gathering_size_restrictions.R")

# ==============================================================================
# 9. Analytical W moments and generalized gamma matching (Section 19)
# ==============================================================================

cat("\n--- 9. Analytical W moments and generalized gamma matching ---\n\n")
source("code/analytical_W_moments.R")

# ==============================================================================
# Done
# ==============================================================================

cat("\n==============================================================\n")
cat("  Pipeline complete.\n")
cat("  Figures saved to figures/\n")
cat("==============================================================\n")
