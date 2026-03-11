## Offspring distribution differs across profiles despite identical $R_0$

The smooth and spike cases both produce Poisson($R_0$) secondary infections per individual. The stepwise case produces a compound Poisson — Poisson($\beta \cdot D$) where $D \sim \text{Exp}(1/i_{dur})$ — which is overdispersed relative to Poisson($R_0$). Concretely, $\text{Var}[\text{offspring}] = R_0 + R_0^2$ for the stepwise case (a negative binomial with $k = 1$), versus $\text{Var}[\text{offspring}] = R_0$ for the smooth and spike cases.

This means the stepwise case — which is the implicit individual-level interpretation of the standard SEIR model — carries superspreading-like overdispersion even without any explicit heterogeneity in contact rates or susceptibility.

## Extinction probability: identical for smooth vs. spike, elevated for stepwise

Because the smooth and spike cases share the same offspring distribution (Poisson($R_0$)), they have *exactly* the same probability of stochastic extinction in the branching process approximation. This is not merely an empirical observation — it follows from the fact that the extinction probability depends on the offspring distribution, not on the timing of transmission events.

## Spike case takes off slower: a "minimum of order statistics" mechanism

Conditional on establishing, the spike case takes ~0.8–1 day longer to reach early case-count thresholds than the smooth case. The mechanism is clean:

**Single-infector level.** In the smooth case, an infector with $k$ contacts generates $k$ independent transmission times drawn from the generation interval distribution. The *first* onward transmission is $\min(\tau_1, \ldots, \tau_k)$, which is stochastically earlier than a single draw. In the spike case, all $k$ contacts occur at the same time $\tau_i^*$, so the first transmission time is just a single generation-interval draw regardless of $k$.

First-transmission time from a single infector, stratified by number of contacts ($10^5$ replicates; $E[\text{GenInterval}] = e_{dur} + i_{dur} = 5.0$ days):

| Contacts ($k$) | Smooth | Spike | Difference |
|-----------------|--------|-------|------------|
| 1 | 4.98 | 4.97 | ~0 (both are single draws) |
| 2 | 3.12 | 5.03 | 1.9 days |
| 3 | 2.39 | 4.99 | 2.6 days |
| 4 | 1.98 | 4.99 | 3.0 days |
| 5 | 1.75 | 4.99 | 3.2 days |

The spike first-transmission time equals $E[\text{GenInterval}]$ regardless of $k$, confirming the mechanism. Averaging over $k \sim \text{Poisson}(R_0)$, the smooth case produces a first transmission 0.88 days earlier overall (3.34 vs. 4.23 days; $p < 10^{-15}$, KS test).

## Early epidemic growth: smooth vs. spike across $R_0$

Punctuated (spike) profiles produce slower early epidemic growth than smooth profiles, via a minimum-of-order-statistics mechanism. In the smooth case, an infector's $k$ infection attempts race independently across the generation interval, so the first successful transmission happens at time $\min(\tau_1, \ldots, \tau_k)$, which shrinks with $k$. In the spike case, all $k$ infection attempts occur at the same time, so the first-transmission time is a single draw from the generation interval regardless of $k$.

### Time to 100 infections (established epidemics)

Simulations with $N = 1000$, $e_\text{dur} = 2$, $i_\text{dur} = 3$ (mean generation interval = 5 days), 2000--3000 runs per ($R_0$, profile) pair, conditioning on epidemics reaching $\geq 100$ infections:

| $R_0$ | $\bar{t}_{10}$ smooth | $\bar{t}_{10}$ spike | delay (t10) | $\bar{t}_{100}$ smooth | $\bar{t}_{100}$ spike | delay (t100) |
|-------|----------------------|---------------------|-------------|----------------------|---------------------|-------------|
| 1.2   | 17.7                 | 17.2                | ~0%         | 65.5                 | 63.7                | ~0%         |
| 1.5   | 15.7                 | 16.9                | +7.6%       | 43.1                 | 44.1                | +2.1%       |
| 2.0   | 12.7                 | 13.5                | +6.4%       | 27.6                 | 28.4                | +2.8%       |
| 3.0   | 8.8                  | 10.2                | +16.7%      | 17.1                 | 18.5                | +8.4%       |
| 5.0   | 5.6                  | 7.4                 | +32.4%      | 10.5                 | 12.3                | +17.6%      |

### Variability in time to 100 infections

The spike case also produces substantially more **variable** epidemic timelines, with the effect again scaling with $R_0$:

| $R_0$ | SD($t_{100}$) smooth | SD($t_{100}$) spike | spike/smooth ratio |
|-------|---------------------|--------------------|--------------------|
| 1.2   | 22.9                | 22.8               | 1.00               |
| 1.5   | 14.6                | 14.9               | 1.03               |
| 2.0   | 8.4                 | 9.0                | 1.07               |
| 3.0   | 4.4                 | 6.0                | 1.34               |
| 5.0   | 2.4                 | 4.3                | 1.82               |

At $R_0 = 5$, the spike case has nearly twice the standard deviation in time-to-100 as the smooth case. The mechanism is the same order-statistics effect operating on variance rather than the mean: in the smooth case, each infector's first-transmission time is $\min(\tau_1, \ldots, \tau_k)$, which concentrates around its mean as $k$ grows (the minimum of many independent draws has lower variance than a single draw). In the spike case, the first-transmission time is always a single draw from the generation interval, preserving the full variance regardless of $k$. This per-generation variance compounds across the chain of transmission, producing more variable epidemic trajectories.

### R0 dependence

The growth delay from punctuation **increases** with $R_0$, contrary to a naive expectation. The mechanism: the order-statistics advantage of smooth profiles scales with the expected number of infection attempts per infector, $E[k] = R_0$. At $R_0 = 1.2$, most infectors draw only $k = 1$ attempt (P(k=1) $\approx 0.36$, P(k=0) $\approx 0.30$), so there is no order-statistics effect to exploit — the two profiles produce indistinguishable growth. At $R_0 = 5$, infectors typically draw $k = 4 \text{-} 6$ attempts, and $\min(\tau_1, \ldots, \tau_5)$ is substantially earlier than a single $\tau$, giving smooth profiles a compounding head start at each generation.

The effect is strongest at the earliest milestones: the t10 delay at $R_0 = 5$ is 32%, attenuating to 18% by t100 as later generations compound and susceptible depletion begins to equalize growth rates.

### Intuition

At high $R_0$, smooth profiles benefit from a "fastest horse wins" dynamic — many independent infection attempts race across the generation interval, and only the fastest one matters for determining the next generation's timing. Spike profiles forfeit this advantage by concentrating all attempts at one time point. The result is that punctuation slows early growth most precisely when $R_0$ is large, i.e., when the epidemic would otherwise grow fastest.

## Establishment timing variance is higher for spike

The same burstiness that slows early takeoff also makes the timing more variable. The spike case shows 22–48% higher variance in the time to reach early thresholds (5,000 simulations):

| Threshold | Smooth | Spike | Variance ratio |
|-----------|--------|-------|----------------|
| $N = 5$   | $8.5 \pm 4.7$ | $9.3 \pm 5.7$ | **1.48** ($p < 10^{-15}$) |
| $N = 10$  | $12.7 \pm 5.9$ | $13.6 \pm 6.9$ | **1.41** ($p < 10^{-15}$) |
| $N = 20$  | $17.1 \pm 6.8$ | $18.0 \pm 7.7$ | **1.27** ($p < 10^{-14}$) |
| $N = 50$  | $22.9 \pm 7.6$ | $23.7 \pm 8.4$ | **1.22** ($p < 10^{-9}$) |

## Spike epidemics have more variable empirical exponential growth rates

- **Exponential growth rate** : Nearly identical means (0.166 vs. 0.170), but the spike case has 13% higher SD (0.053 vs. 0.047), consistent with the higher-variance theme. Based on log-linear fit, cases 5–50. 

## Overdispersion from punctuated infectiousness x periodic contacts

We decompose individual infectiousness as $a_i(\tau) = R_0 \cdot b_i(\tau) \cdot z(t_i + \tau)$, where:

- $b_i(\tau) = f_\kappa(\tau - s_i)$ is the biological timing density (a Gamma($\kappa$, $r$) density shifted by individual onset $s_i \sim \text{Gamma}(\alpha - \kappa, r)$)
- $z(t) = 1 + \epsilon \sin(2\pi t / T)$ is a periodic contact rate multiplier with period $T = 7$ days

The parameter $\kappa$ controls punctuation: small $\kappa$ gives narrow, spike-like individual profiles; large $\kappa$ gives broad profiles that resemble the population average. The parameter $\epsilon$ controls the amplitude of contact variation ($\epsilon = 0$ is constant contacts, $\epsilon = 0.9$ is strong weekly oscillation).

The population-level kernel $A(\tau) = R_0 \cdot \text{Gamma}(\tau; \alpha, r)$ is invariant across $\kappa$ by the Gamma additivity property. Parameters: $R_0 = 2$, $\alpha = 10$, $\mu = 5$ days, so $r = 2$.

### Part A: Branching process offspring distributions

For each ($\kappa$, $\epsilon$) pair, we drew 10,000 individuals with random infection phase $t_i \sim \text{Uniform}(0, T)$ and computed the realized reproduction number $\nu_i = R_0 \int b_i(\tau) \, z(t_i + \tau) \, d\tau$ via numerical integration, then drew offspring $n_i \sim \text{Poisson}(\nu_i)$.

#### Verification: mean offspring

The mean offspring count is $\approx 2.0$ across all ($\kappa$, $\epsilon$) combinations (range: 1.91--2.02). This confirms that periodic contacts redistribute infectiousness across individuals but do not change the mean. The modulation is purely a second-order effect.

#### Verification: Poisson baseline

At $\epsilon = 0$ (constant contacts), the variance/mean ratio is $\approx 1.0$ for all $\kappa$ values (range: 0.98--1.01). This confirms that without contact variation, the offspring distribution is Poisson regardless of how punctuated the biological profile is. Punctuation alone does not create overdispersion in a branching process with constant contacts.

#### Variance/mean ratio

The core finding: overdispersion arises from the *interaction* between punctuated infectiousness and periodic contacts. Neither ingredient alone is sufficient.

| $\kappa$ | $\epsilon = 0$ | $\epsilon = 0.3$ | $\epsilon = 0.5$ | $\epsilon = 0.7$ | $\epsilon = 0.9$ |
|----------|----------------|-------------------|-------------------|-------------------|-------------------|
| 0.5      | 1.00           | 1.10              | 1.20              | 1.41              | 1.69              |
| 1        | 0.98           | 1.08              | 1.19              | 1.41              | 1.71              |
| 2        | 0.99           | 1.05              | 1.21              | 1.37              | 1.60              |
| 4        | 0.98           | 1.00              | 1.12              | 1.27              | 1.42              |
| 6        | 1.00           | 1.03              | 1.06              | 1.17              | 1.24              |
| 8        | 1.00           | 1.01              | 1.06              | 1.12              | 1.16              |
| 9.5      | 1.01           | 1.05              | 1.04              | 1.07              | 1.12              |

The overdispersion is concentrated in the low-$\kappa$, high-$\epsilon$ corner: punctuated profiles ($\kappa \leq 2$) with strong contact oscillation ($\epsilon \geq 0.7$) yield variance/mean ratios of 1.4--1.7. Smooth profiles ($\kappa = 9.5$) show at most a ratio of 1.12 even at maximum contact amplitude.

#### Negative binomial dispersion parameter $k$

Fitting $k$ via method of moments ($k = \mu^2 / (\text{Var} - \mu)$):

| $\kappa$ | $\epsilon = 0$ | $\epsilon = 0.5$ | $\epsilon = 0.7$ | $\epsilon = 0.9$ |
|----------|----------------|-------------------|-------------------|-------------------|
| 0.5      | $\infty$       | 9.7               | 4.7               | 2.8               |
| 1        | $\infty$       | 10.3              | 4.9               | 2.8               |
| 2        | $\infty$       | 9.6               | 5.4               | 3.4               |
| 4        | $\infty$       | 17.2              | 7.3               | 4.7               |
| 6        | $\infty$       | 30.9              | 11.8              | 8.4               |
| 8        | $\infty$       | 32.5              | 16.5              | 12.4              |
| 9.5      | $\infty$       | 51.0              | 27.7              | 16.8              |

At $\epsilon = 0$, all $\kappa$ values give $k = \infty$ (Poisson). At $\epsilon = 0.9$ with $\kappa = 0.5$, we get $k \approx 2.8$ -- substantial overdispersion from a purely mechanistic origin, comparable to values estimated for some respiratory infections. The gradient from $k \approx 17$ (smooth + periodic) to $k \approx 3$ (punctuated + periodic) demonstrates that punctuation amplifies contact-driven overdispersion by roughly a factor of 6 in the dispersion parameter.

#### Mechanism

The mechanism is transparent: when the biological profile $b_i$ is a narrow spike (small $\kappa$), the integral $\int b_i(\tau) \, z(t_i + \tau) \, d\tau$ essentially samples $z$ at one point, so $\nu_i$ varies between $R_0(1 - \epsilon)$ and $R_0(1 + \epsilon)$ depending on when the spike falls relative to the contact cycle. When $b_i$ is broad (large $\kappa$), it averages over the contact cycle, and $\nu_i \approx R_0$ for everyone.

This is a sampling-vs-averaging effect: punctuated profiles *sample* the contact function; smooth profiles *average* over it.

### Part B: Figures

Generated in `figures/`:

- **`fig_overdispersion_heatmap.pdf`**: Var/mean ratio across the ($\kappa$, $\epsilon$) grid. Overdispersion is concentrated in the low-$\kappa$, high-$\epsilon$ corner, confirming the interaction.
- **`fig_overdispersion_lines.pdf`**: Var/mean ratio vs $\kappa$ for selected $\epsilon$ values. All curves converge to $\approx 1$ at high $\kappa$, and fan out at low $\kappa$ proportional to $\epsilon$.
- **`fig_overdispersion_histograms.pdf`**: Offspring distributions for four corner cases compared to Poisson($R_0$). The punctuated + periodic case has visibly heavier tails than the Poisson reference; the other three are nearly indistinguishable from Poisson.
- **`fig_overdispersion_nbk.pdf`**: Negative binomial $k$ parameter across the grid (lower = more overdispersed).

### Part C: Full epidemic simulations

200 stochastic epidemics per scenario ($N = 1000$, $R_0 = 2$), using `sim_stochastic_fast` with the `make_profile_gamma_contacts` factory.

#### Establishment probability

| Scenario              | P(established) |
|-----------------------|----------------|
| Smooth + constant     | 0.75           |
| Punctuated + constant | 0.77           |
| Smooth + periodic     | 0.63           |
| Punctuated + periodic | 0.55           |

The punctuated + periodic scenario has the lowest establishment probability (0.55 vs 0.75 for the baseline). This is consistent with the overdispersed offspring distribution: higher variance in $\nu_i$ means more individuals with $\nu_i < 1$ who fail to sustain transmission, increasing the probability of early stochastic extinction.

#### Final size (established epidemics only)

| Scenario              | N established | Mean final size | SD   |
|-----------------------|---------------|-----------------|------|
| Smooth + constant     | 150           | 796             | 19.6 |
| Punctuated + constant | 153           | 796             | 20.9 |
| Smooth + periodic     | 126           | 727             | 26.7 |
| Punctuated + periodic | 110           | 730             | 30.0 |

Among established epidemics, the final size is comparable between constant-contact scenarios ($\approx 796$) but lower for periodic scenarios ($\approx 728$). This is because the time-varying contact rate creates periods of suppressed transmission that slow propagation and allow more individuals to remain susceptible when the epidemic wanes. The punctuated + periodic scenario also has the highest final-size variability (SD = 30 vs 20 for baseline), consistent with its overdispersed dynamics.

#### Epidemic trajectories

Example curves are in `figures/fig_epidemic_curves.pdf`. The constant-contact scenarios (both smooth and punctuated) produce classic sigmoid epidemic curves. The periodic scenarios show a subtle modulation in growth rate. The key difference is not in the shape of individual trajectories but in the *distribution* across simulation runs: the punctuated + periodic case has more variable outcomes, reflected in the wider spread of final sizes and the lower establishment probability.

## Summary of overdispersion and early-growth findings

1. **Punctuation alone does not create overdispersion** in a branching process with constant contacts. The offspring distribution is Poisson($R_0$) for all $\kappa$ when $\epsilon = 0$.

2. **Contact variation alone creates mild overdispersion** even with smooth profiles. At $\epsilon = 0.9$ and $\kappa = 9.5$, the variance/mean ratio is 1.12 ($k \approx 17$).

3. **The interaction amplifies overdispersion**: punctuated profiles + periodic contacts yield variance/mean $\approx 1.7$ and $k \approx 3$ -- a 6-fold reduction in the dispersion parameter compared to smooth profiles with the same contact variation.

4. **The mechanism is sampling vs averaging**: narrow individual profiles sample the contact function at a random phase; broad profiles average over it. This is a direct consequence of the decomposition $a_i(\tau) = R_0 \cdot b_i(\tau) \cdot z(t_i + \tau)$.

5. **Epidemic consequences**: the overdispersion reduces establishment probability (0.55 vs 0.75) and increases variability in epidemic outcomes, without changing the mean offspring count.

6. **Punctuation slows early growth, amplified by $R_0$**: spike profiles reach 100 infections 3--18% slower than smooth profiles (at $R_0 = 2 \text{-} 5$), via a minimum-of-order-statistics mechanism. The effect grows with $R_0$ because more infection attempts per infector means a larger order-statistics advantage for smooth profiles. At $R_0 \leq 1.2$, the effect vanishes because most infectors only generate one attempt.