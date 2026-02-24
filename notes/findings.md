# Findings log

Running log of results for the manuscript. All simulations use SEIR parameters ($e\_dur = 2$, $i\_dur = 3$, $R_0 = 2$, $N = 1000$) unless otherwise noted. The three individual infectiousness profiles — smooth (all-equal), stepwise, and spike (delta function) — all share the same population-level infectiousness kernel $A(\tau)$ and the same $R_0$.

---

## 1. Offspring distribution differs across profiles despite identical $R_0$

The smooth and spike cases both produce Poisson($R_0$) secondary infections per individual. The stepwise case produces a compound Poisson — Poisson($\beta \cdot D$) where $D \sim \text{Exp}(1/i\_dur)$ — which is overdispersed relative to Poisson($R_0$). Concretely, $\text{Var}[\text{offspring}] = R_0 + R_0^2$ for the stepwise case (a negative binomial with $k = 1$), versus $\text{Var}[\text{offspring}] = R_0$ for the smooth and spike cases.

This means the stepwise case — which is the implicit individual-level interpretation of the standard SEIR model — carries superspreading-like overdispersion even without any explicit heterogeneity in contact rates or susceptibility.

## 2. Extinction probability: identical for smooth vs. spike, elevated for stepwise

Because the smooth and spike cases share the same offspring distribution (Poisson($R_0$)), they have *exactly* the same probability of stochastic extinction in the branching process approximation. This is not merely an empirical observation — it follows from the fact that the extinction probability depends on the offspring distribution, not on the timing of transmission events.

Confirmed empirically (5,000 simulations each, $N = 1000$):
- Smooth: ~19.5% extinction (final size $< 10$)
- Spike: ~20.1% extinction

The stepwise case should have a higher extinction probability due to its overdispersed offspring distribution (a well-known result: for a fixed mean, higher variance in the offspring distribution increases extinction probability).

## 3. Spike case takes off slower: a "minimum of order statistics" mechanism

Conditional on establishing, the spike case takes ~0.8–1 day longer to reach early case-count thresholds than the smooth case. The mechanism is clean:

**Single-infector level.** In the smooth case, an infector with $k$ contacts generates $k$ independent transmission times drawn from the generation interval distribution. The *first* onward transmission is $\min(\tau_1, \ldots, \tau_k)$, which is stochastically earlier than a single draw. In the spike case, all $k$ contacts occur at the same time $\tau_i^*$, so the first transmission time is just a single generation-interval draw regardless of $k$.

First-transmission time from a single infector, stratified by number of contacts ($10^5$ replicates; $E[\text{GenInterval}] = e\_dur + i\_dur = 5.0$ days):

| Contacts ($k$) | Smooth | Spike | Difference |
|-----------------|--------|-------|------------|
| 1 | 4.98 | 4.97 | ~0 (both are single draws) |
| 2 | 3.12 | 5.03 | 1.9 days |
| 3 | 2.39 | 4.99 | 2.6 days |
| 4 | 1.98 | 4.99 | 3.0 days |
| 5 | 1.75 | 4.99 | 3.2 days |

The spike first-transmission time equals $E[\text{GenInterval}]$ regardless of $k$, confirming the mechanism. Averaging over $k \sim \text{Poisson}(R_0)$, the smooth case produces a first transmission 0.88 days earlier overall (3.34 vs. 4.23 days; $p < 10^{-15}$, KS test).

**Early epidemic level.** This per-infector advantage compounds: the delay concentrates in the very first inter-event time (infection 1 $\to$ 2: 3.3 days for smooth vs. 5.1 days for spike), then rapidly washes out as multiple infectors become active simultaneously. By infection 5 $\to$ 6, inter-event times are indistinguishable (~1.14 days for both).

## 4. Establishment timing variance is higher for spike

The same burstiness that slows early takeoff also makes the timing more variable. The spike case shows 22–48% higher variance in the time to reach early thresholds (5,000 simulations):

| Threshold | Smooth | Spike | Variance ratio |
|-----------|--------|-------|----------------|
| $N = 5$   | $8.5 \pm 4.7$ | $9.3 \pm 5.7$ | **1.48** ($p < 10^{-15}$) |
| $N = 10$  | $12.7 \pm 5.9$ | $13.6 \pm 6.9$ | **1.41** ($p < 10^{-15}$) |
| $N = 20$  | $17.1 \pm 6.8$ | $18.0 \pm 7.7$ | **1.27** ($p < 10^{-14}$) |
| $N = 50$  | $22.9 \pm 7.6$ | $23.7 \pm 8.4$ | **1.22** ($p < 10^{-9}$) |

The effect is strongest at the smallest thresholds and diminishes as the law of large numbers smooths out individual-level variation.

## 5. Spike epidemics are burstier: higher temporal clustering of early cases

Even after accounting for the slower start, the spike case produces qualitatively different temporal structure in early case arrivals (5,000 simulations, first 20 infections):

- **Coefficient of variation of inter-event times**: 1.89 (spike) vs. 1.23 (smooth) — a 54% increase. The spike case alternates between quiet gaps and bursts; the smooth case has more regular spacing.
- **Same-day burst fraction**: 64.6% of early infections in the spike case arrive on the same calendar day as the previous infection, versus 47.5% for smooth.
- **Exponential growth rate** (log-linear fit, cases 5–50): nearly identical means (0.166 vs. 0.170), but the spike case has 13% higher SD (0.053 vs. 0.047), consistent with the higher-variance theme.

Analysis code: `code/early_growth_analysis.R`.

## 6. Large epidemics converge to the same deterministic trajectory

By construction, all three profiles share the same population-level $A(\tau)$ and therefore the same mean-field ODE dynamics. Empirically, the mean cumulative infection curves for established epidemics converge to the SEIR ODE solution for all three profiles. Differences between profiles are purely stochastic and vanish in the large-population limit.

## 7. Within-individual correlation structure differs between smooth and spike

In the smooth case, a single infector generates contacts at multiple different generation intervals (spanning the full generation interval distribution). An individual who transmits early also transmits late. In the spike case, all of a person's transmissions occur at the same generation interval $\tau_i^*$. This means:

- In the smooth case, generation intervals from a common infector are positively correlated (they are draws from the same profile, conditional on the infector existing).
- In the spike case, generation intervals from a common infector are identical (all equal $\tau_i^*$), but generation intervals from *different* infectors are independent.

This distinction should produce different phylogenetic tree structures (more synchronized branching in the spike case) and different temporal autocorrelation patterns in incidence data, even when aggregate case counts are similar.

## 8. Intervention impact variance differs across profiles

Consider a one-time intervention that removes a fraction of infectious individuals at time $t$ (e.g., a quarantine sweep). In the smooth case, each removed individual has a predictable amount of remaining infectiousness. In the spike case, a removed individual has either already transmitted all $R_0$ of their infections (if $t > t_i + \tau_i^*$) or has transmitted none of them (if $t < t_i + \tau_i^*$). The expected impact of removal is similar across profiles, but the *variance* in impact is higher for the spike case — interventions are a higher-variance gamble when infectiousness is punctuated.

*Status: theoretical prediction, not yet tested in simulation.*

---

*Last updated: 2026-02-23*
