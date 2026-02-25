# ==============================================================================
# Master analysis script — orchestrates all analyses
#
# Run from the project root:
#   Rscript code/run_analysis.R
#
# Scripts are ordered to match the writeup (writeup/punctuatedinfectiousness.md):
#
#   1. Setup              — load libraries, ensure output directories exist
#   2. simulate_epidemics — uncontrolled epidemic dynamics under different
#                           individual infectiousness profiles (Sections 3–4)
#   3. visualize_profiles_gamma — Gamma convolutional model visualization
#                           (Section 5)
#   4. test_based_screening — impact of punctuated infectiousness on
#                           test-based screening (Section 6)
#   5. analyze_overdispersion — biological/contact decomposition and
#                           overdispersion from periodic contacts (Sections 7–8)
#   6. simulate_trials    — identifiability of kappa from household/contact-
#                           tracing data (Section 9)
# ==============================================================================

cat("==============================================================\n")
cat("  Punctuated Infectiousness — Full Analysis Pipeline\n")
cat("==============================================================\n\n")

# ==============================================================================
# 0. Setup
# ==============================================================================

cat("--- Setup ---\n")
if (!dir.exists("figures")) dir.create("figures")

# ==============================================================================
# 1. Uncontrolled epidemic simulations (Sections 3–4)
# ==============================================================================

cat("\n--- 1. Simulating uncontrolled epidemics ---\n\n")
source("code/simulate_epidemics.R")

# ==============================================================================
# 2. Gamma profile visualization (Section 5)
# ==============================================================================

cat("\n--- 2. Visualizing Gamma infectiousness profiles ---\n\n")
source("code/visualize_profiles_gamma.R")

# ==============================================================================
# 3. Test-based screening (Section 6)
# ==============================================================================

cat("\n--- 3. Test-based screening analysis ---\n\n")
source("code/test_based_screening.R")

# ==============================================================================
# 4. Overdispersion from periodic contacts (Sections 7–8)
# ==============================================================================

cat("\n--- 4. Overdispersion analysis ---\n\n")
source("code/analyze_overdispersion.R")

# ==============================================================================
# 5. Identifiability / trial simulations (Section 9)
# ==============================================================================

cat("\n--- 5. Trial simulations for identifiability ---\n\n")
source("code/simulate_trials.R")

# ==============================================================================
# Done
# ==============================================================================

cat("\n==============================================================\n")
cat("  Pipeline complete.\n")
cat("  Figures saved to figures/\n")
cat("==============================================================\n")
