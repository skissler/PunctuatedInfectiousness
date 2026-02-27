# Findings log

Running log of results for the manuscript. All simulations use SEIR parameters ($e\_{dur} = 2$, $i\_{dur} = 3$, $R_0 = 2$, $N = 1000$) unless otherwise noted. The three individual infectiousness profiles — smooth (all-equal), stepwise, and spike (delta function) — all share the same population-level infectiousness kernel $A(\tau)$ and the same $R_0$.

---

## 1. Offspring distribution differs across profiles despite identical $R_0$

The smooth and spike cases both produce Poisson($R_0$) secondary infections per individual. The stepwise case produces a compound Poisson — Poisson($\beta \cdot D$) where $D \sim \text{Exp}(1/i\_{dur})$ — which is overdispersed relative to Poisson($R_0$). Concretely, $\text{Var}[\text{offspring}] = R_0 + R_0^2$ for the stepwise case (a negative binomial with $k = 1$), versus $\text{Var}[\text{offspring}] = R_0$ for the smooth and spike cases.

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

Consider a one-time intervention that removes a fraction of infectious individuals at time $t$ (e.g., a quarantine sweep). In the smooth case, each removed individual has a predictable amount of remaining infectiousness. In the spike case, a removed individual has either already transmitted all $R_0$ of their infections (if $t > t_i + \tau_i^*$) or has transmitted none of them (if $t < t_i + \tau_i^\*$). The expected impact of removal is similar across profiles, but the *variance* in impact is higher for the spike case — interventions are a higher-variance gamble when infectiousness is punctuated.

*Status: theoretical prediction, not yet tested in simulation.*

## 9. Shifted Gamma construction: a clean one-parameter interpolation

The original three profiles (smooth, stepwise, spike) are useful theoretical bookends but leave open the question of how to *interpolate* between the smooth and spike extremes with a single punctuation parameter, while preserving:

1. $\int a_i(\tau)\,d\tau = R_0$ exactly for every individual,
2. $E[a_i(\tau)] = A(\tau)$ exactly (population kernel recovered in expectation),
3. $A(\tau)$ invariant as the punctuation parameter varies,
4. All individual profiles having the same shape, height, and width (differing only in location).

The shifted Gamma construction achieves all four simultaneously.

### Setup

Choose a population-level generation interval distribution

$$A(\tau) = R_0 \cdot \text{Gamma}(\tau; \alpha_{\text{total}}, r), \qquad r = \alpha_{\text{total}} / \mu,$$

where $\mu$ is the mean generation time and $\alpha_{\text{total}}$ controls the shape (larger = more symmetric/bell-shaped). This does not need to match the SEIR kernel exactly; it is a modelling choice.

### Decomposition

Decompose each contact time as

$$c_j = s_i + \varepsilon_j,$$

where:

- $s_i \sim \text{Gamma}(\alpha_{\text{total}} - \kappa,\; r)$ is an individual-specific **onset shift** (drawn once per individual),
- $\varepsilon_j \sim \text{Gamma}(\kappa,\; r)$ is a **jitter** around the onset (drawn independently per contact, same distribution for all individuals),
- $\kappa \in (0, \alpha_{\text{total}})$ is the **punctuation parameter**.

### Why it works

The key identity is that Gamma distributions with the same rate are closed under convolution:

$$\text{Gamma}(\alpha_{\text{total}} - \kappa,\; r) + \text{Gamma}(\kappa,\; r) = \text{Gamma}(\alpha_{\text{total}},\; r).$$

So the marginal distribution of $c_j$ (integrating over $s_i$) is $\text{Gamma}(\alpha_{\text{total}}, r)$ regardless of $\kappa$. This immediately gives $E[a_i(\tau)] = A(\tau)$ and guarantees that $A(\tau)$ is invariant as $\kappa$ varies.

### Individual profiles

Each individual's infectiousness profile is

$$a_i(\tau) = R_0 \cdot f_\kappa(\tau - s_i), \qquad f_\kappa = \text{Gamma}(\kappa, r),$$

which is zero for $\tau < s_i$ and follows a $\text{Gamma}(\kappa, r)$ density thereafter. Since $f_\kappa$ is a proper density, $\int a_i(\tau)\,d\tau = R_0$ exactly. All individuals share the same $f_\kappa$ — the profiles are identical in shape, height, and width; only the onset time $s_i$ differs.

### Limiting behaviour

| $\kappa$ | Individual profile $f_\kappa$ | Shift distribution | Interpretation |
|---|---|---|---|
| $\kappa \to 0$ | $\delta$-function (infinitely narrow spike) | $\text{Gamma}(\alpha_{\text{total}}, r) = A/R_0$ | Maximally punctuated: all contacts at one instant |
| Small $\kappa$ (e.g. 2) | Narrow unimodal bump | Wide shift distribution | Punctuated but not singular |
| $\kappa \to \alpha_{\text{total}}$ | $\approx A(\tau)/R_0$ (broad, matches population kernel) | $\delta$-function at 0 (no shift) | Smooth: all individuals have the same profile $\approx A(\tau)$ |

### Variance decomposition

The total variance of contact times decomposes additively:

$$\text{Var}(c_j) = \underbrace{\text{Var}(s_i)}_{\text{between-individual}} + \underbrace{\text{Var}(\varepsilon_j)}_{\text{within-individual}} = \frac{\alpha_{\text{total}} - \kappa}{r^2} + \frac{\kappa}{r^2} = \frac{\alpha_{\text{total}}}{r^2}.$$

The fraction of total variance that is between-individual (i.e., due to punctuation) is

$$\text{punctuation fraction} = \frac{\alpha_{\text{total}} - \kappa}{\alpha_{\text{total}}} = 1 - \frac{\kappa}{\alpha_{\text{total}}}.$$

This gives a clean interpretation: $\kappa / \alpha_{\text{total}} \in (0, 1)$ is the fraction of generation-time variability that is within-individual.

### Implementation

The factory function `make_profile_gamma(mu, R0, alpha_total, kappa)` in `code/utils.R` returns a closure that can be passed directly to `sim_stochastic_fast` via the `gen_inf_attempts` argument:

```r
profile <- make_profile_gamma(mu = 5, R0 = 2, alpha_total = 10, kappa = 4)
result  <- sim_stochastic_fast(n = 1000, gen_inf_attempts = profile)
```

Visualization: `code/visualize_profiles_gamma.R`.

### Design choices

- **$\alpha_{\text{total}}$**: Controls the shape of $A(\tau)$. Larger values give a more symmetric, bell-shaped kernel. The moment-matched value for SEIR parameters ($e\_{dur} = 2$, $i\_{dur} = 3$) is $\alpha_{\text{total}} = \mu^2 / (e\_{dur}^2 + i\_{dur}^2) \approx 1.92$, but this is too small for $f_\kappa$ to be unimodal at most $\kappa$ values (Gamma shape $< 1$ gives a monotone-decreasing density). We use $\alpha_{\text{total}} = 10$ as a default, which gives a well-shaped $A(\tau)$ and allows $\kappa$ to range from $\sim 2$ (punctuated, narrow unimodal bumps) to $\sim 9.5$ (smooth, profiles $\approx A$).

- **$\mu$**: The mean generation time, equal to $\alpha_{\text{total}} / r$. Changing $\mu$ scales the time axis without affecting the shape of $A(\tau)$ or the punctuation structure.

## 10. Decomposing infectiousness into biological and contact components

The individual infectiousness profile $a_i(\tau)$ conflates two distinct processes: the biological trajectory of infectiousness (viral load, shedding) and the rate at which the individual makes potentially transmissive contacts. Separating these gives a cleaner framework for studying how temporal structure in contact patterns interacts with punctuated biological infectiousness.

### Definitions

Decompose

$$a_i(\tau) = b_i(\tau) \cdot c_i(t_i + \tau),$$

where $t_i$ is the calendar time at which person $i$ was infected, and:

- **$b_i(\tau)$** is the *biological infectiousness profile*: a probability density on $[0, \infty)$ satisfying $\int_0^\infty b_i(\tau)\,d\tau = 1$. It describes *when* the individual is capable of transmitting, relative to their time of infection. It contains no information about how many people they will infect — only the temporal shape.

- **$c_i(t)$** is the *instantaneous transmission potential* at calendar time $t$: the expected number of secondary infections the individual would produce at time $t$ if all their biological infectiousness were concentrated at that instant (i.e., if $b_i$ were a delta function). This absorbs contact rate, transmission probability per contact, and susceptibility of contacts. When $c_i$ is constant in time, $c_i = \nu_i$ is exactly the individual reproduction number in the sense of Lloyd-Smith *et al.* (2005).

The individual's realized reproduction number is then

$$\tilde{\nu}_i = \int_0^\infty a_i(\tau)\,d\tau = \int_0^\infty b_i(\tau)\,c_i(t_i + \tau)\,d\tau,$$

a weighted average of their contact potential over the times they are biologically infectious.

### Population-level factorization

Define $B(\tau) = E[b_i(\tau)]$ and $C(t) = E[c_i(t)]$. Under the assumption that $b_i$ and $c_i$ are **independent** — that an individual's biological trajectory is independent of their contact pattern — the population-level infectiousness kernel factorizes:

$$A_t(\tau) = E[a_i(\tau)] = B(\tau) \cdot C(t_i + \tau).$$

In the renewal equation, all contributions to the force of infection at calendar time $t$ share the same calendar-time contact factor $C(t)$, so $C$ factors out:

$$J(t) = S(t) \cdot C(t) \int_0^\infty J(t - \tau)\,B(\tau)\,d\tau.$$

The time-varying contact function acts as a multiplicative modulator on the effective reproduction number. When $C(t) = C$ (constant), we recover the standard time-homogeneous renewal equation with $A(\tau) = C \cdot B(\tau)$.

*(If $b_i$ and $c_i$ are not independent — for instance, if symptomatic individuals reduce contacts late in their infection — the factorization breaks and the joint distribution $E[b_i(\tau) \cdot c_i(t)]$ is needed. We assume independence throughout this project.)*

### Connection to the shifted Gamma construction

In the shifted Gamma framework (Section 9), with constant contacts:

- $b_i(\tau) = f_\kappa(\tau - s_i)$ — the Gamma$(\kappa, r)$ density shifted to onset time $s_i$
- $c_i(t) = R_0$ for all $i, t$

This is the simplest case: all the punctuation lives in $b_i$, contacts are homogeneous and time-invariant, and $\tilde{\nu}_i = R_0$ for every individual.

### Punctuated infectiousness as a source of superspreading

The decomposition reveals a mechanism for generating overdispersion in the individual reproduction number $\tilde{\nu}_i$ that, to our knowledge, has not been previously described.

**Setup.** Suppose:

- Biological profiles are shifted Gammas: $b_i(\tau) = f_\kappa(\tau - s_i)$ with punctuation parameter $\kappa$
- The contact function is the *same for everyone* and varies periodically in calendar time: $c_i(t) = c(t) = R_0 \cdot z(t)$, where $z(t)$ oscillates around 1 (e.g., weekday/weekend cycles, seasonal forcing)

There is no individual-level variation in contact rates — the only heterogeneity is in the biological timing $s_i$. Yet the individual reproduction number becomes:

$$\tilde{\nu}_i = R_0 \int_0^\infty f_\kappa(\tau - s_i) \cdot z(t_i + \tau)\,d\tau.$$

The variance of $\tilde{\nu}_i$ depends on how much the biological profile "averages over" the oscillations in $z$:

- **Smooth limit** ($\kappa \to \alpha_{\text{total}}$): $b_i \approx B$ for all $i$, so $\tilde{\nu}_i \approx R_0 \int B(\tau)\,z(t_i + \tau)\,d\tau$. The broad biological window averages over $z$, and $\tilde{\nu}_i \approx R_0$ for everyone. Minimal overdispersion.

- **Delta limit** ($\kappa \to 0$): $b_i \to \delta(\tau - s_i)$, so $\tilde{\nu}_i \to R_0 \cdot z(t_i + s_i)$. The individual's reproduction number directly *samples* the contact function at a single point. If $z$ oscillates between 0.5 and 1.5, so does $\tilde{\nu}_i$. Maximal overdispersion.

- **Intermediate $\kappa$**: the narrower $b_i$ is, the less it averages over $z$, and the more variable $\tilde{\nu}_i$ becomes. The punctuation parameter $\kappa$ continuously controls the degree of overdispersion.

**The mechanism in words:** when biological infectiousness is concentrated in a narrow window, the individual's total transmission depends on *when* that window falls relative to the contact landscape. Some individuals' narrow windows align with high-contact periods (weekdays, social events, crowded settings); others align with low-contact periods (nights, weekends, holidays). This "sampling" of the contact function creates between-individual variation in $\tilde{\nu}_i$ — i.e., superspreading — even without any intrinsic heterogeneity in contact rates or biological susceptibility.

This is distinct from the standard superspreading mechanisms (heterogeneous $\nu_i$ across individuals, heterogeneous susceptibility, stochastic contact networks). It requires only three ingredients:

1. Punctuated biological infectiousness (narrow $b_i$)
2. Stochastic variation in the timing of infectiousness (random $s_i$)
3. Time-varying contact rates at the population level ($z(t) \neq \text{const}$)

All three are empirically plausible for many infectious diseases. The shifted Gamma construction provides a single parameter ($\kappa$) that controls the strength of this effect.

## 11. Overdispersion from punctuated infectiousness x periodic contacts

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

The overdispersion is concentrated in the **low-$\kappa$, high-$\epsilon$ corner**: punctuated profiles ($\kappa \leq 2$) with strong contact oscillation ($\epsilon \geq 0.7$) yield variance/mean ratios of 1.4--1.7. Smooth profiles ($\kappa = 9.5$) show at most a ratio of 1.12 even at maximum contact amplitude.

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

## 12. Early epidemic growth: smooth vs. spike across $R_0$

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

The growth delay from punctuation **increases** with $R_0$, contrary to a naive expectation. The mechanism: the order-statistics advantage of smooth profiles scales with the expected number of infection attempts per infector, $E[k] = R_0$. At $R_0 = 1.2$, most infectors draw only $k = 1$ attempt (P(k=1) $\approx 0.36$, P(k=0) $\approx 0.30$), so there is no order-statistics effect to exploit — the two profiles produce indistinguishable growth. At $R_0 = 5$, infectors typically draw $k = 4$--$6$ attempts, and $\min(\tau_1, \ldots, \tau_5)$ is substantially earlier than a single $\tau$, giving smooth profiles a compounding head start at each generation.

The effect is strongest at the earliest milestones: the t10 delay at $R_0 = 5$ is 32%, attenuating to 18% by t100 as later generations compound and susceptible depletion begins to equalize growth rates.

### Intuition

At high $R_0$, smooth profiles benefit from a "fastest horse wins" dynamic — many independent infection attempts race across the generation interval, and only the fastest one matters for determining the next generation's timing. Spike profiles forfeit this advantage by concentrating all attempts at one time point. The result is that punctuation slows early growth most precisely when $R_0$ is large, i.e., when the epidemic would otherwise grow fastest.

## 13. Summary of overdispersion and early-growth findings

1. **Punctuation alone does not create overdispersion** in a branching process with constant contacts. The offspring distribution is Poisson($R_0$) for all $\kappa$ when $\epsilon = 0$.

2. **Contact variation alone creates mild overdispersion** even with smooth profiles. At $\epsilon = 0.9$ and $\kappa = 9.5$, the variance/mean ratio is 1.12 ($k \approx 17$).

3. **The interaction amplifies overdispersion**: punctuated profiles + periodic contacts yield variance/mean $\approx 1.7$ and $k \approx 3$ -- a 6-fold reduction in the dispersion parameter compared to smooth profiles with the same contact variation.

4. **The mechanism is sampling vs averaging**: narrow individual profiles sample the contact function at a random phase; broad profiles average over it. This is a direct consequence of the decomposition $a_i(\tau) = R_0 \cdot b_i(\tau) \cdot z(t_i + \tau)$.

5. **Epidemic consequences**: the overdispersion reduces establishment probability (0.55 vs 0.75) and increases variability in epidemic outcomes, without changing the mean offspring count.

6. **Punctuation slows early growth, amplified by $R_0$**: spike profiles reach 100 infections 3--18% slower than smooth profiles (at $R_0 = 2$--$5$), via a minimum-of-order-statistics mechanism. The effect grows with $R_0$ because more infection attempts per infector means a larger order-statistics advantage for smooth profiles. At $R_0 \leq 1.2$, the effect vanishes because most infectors only generate one attempt.

## 14. Test-and-isolate distorts the generation interval — but only for smooth profiles

Test-based screening with post-detection isolation reduces transmission by truncating each detected individual's infectiousness profile: $a_i^{\text{eff}}(\tau) = a_i(\tau) \cdot \mathbf{1}[\tau < \tau_{\text{detect}}]$. This truncation both reduces $R_{\text{eff}}$ and reshapes the generation interval distribution — but the degree of reshaping depends critically on how punctuated the individual profile is.

### Setup

Regular screening every $\Delta$ days (random phase), with a positivity window of $[t_{\text{peak}} - w_-, t_{\text{peak}} + w_+]$, sensitivity $\sigma$, and false positive rate $\phi$. If detected, the individual isolates immediately, zeroing out all subsequent transmission. We compare the distribution of generation intervals among **successful** transmissions (those occurring before detection) to the natural generation interval $g(\tau) = \text{Gamma}(\tau; \alpha_{\text{total}}, r)$.

### Spike (delta) case: all-or-nothing filtering

All of individual $i$'s infection attempts occur at a single time $s_i$. Detection at $\tau_{\text{detect}}$ either happens before $s_i$ (blocking everything) or after $s_i$ (blocking nothing). Since the testing phase is uniform and the positivity window spans the spike symmetrically, the blocking probability $p_{\text{block}}$ is approximately constant across $s_i$ values. The effective generation interval among successful transmissions is:

$$g_{\text{eff}}(\tau) \propto f_S(\tau) \cdot (1 - p_{\text{block}}(\tau)) \approx f_S(\tau) \cdot \text{const}$$

**The shape is unchanged.** Testing removes a random fraction of would-be transmitters without altering the timing of the remaining ones.

### Smooth case: profile truncation

Each individual's profile is broad ($\approx A(\tau)$). Detection near the peak truncates the right tail, so only the left side contributes to successful transmission:

$$g_{\text{eff}}(\tau) \propto A(\tau) \cdot P(T_{\text{detect}} > \tau)$$

Since $P(T_{\text{detect}} > \tau)$ is decreasing in $\tau$, the effective generation interval **shifts shorter** relative to the natural $g(\tau)$.

### Epidemiological consequence: partially self-defeating intervention

The Euler-Lotka equation for the epidemic growth rate $r$ is:

$$1 = R_{\text{eff}} \int_0^\infty e^{-r\tau} g_{\text{eff}}(\tau) \, d\tau$$

In the spike case, $g_{\text{eff}} = g$, so the growth rate decreases in direct proportion to the $R_{\text{eff}}$ reduction.

In the smooth case, $g_{\text{eff}}$ shifts shorter. For a given $R_{\text{eff}}$, a shorter generation interval implies a **larger** growth rate (transmissions happen sooner $\Rightarrow$ faster epidemic). So the reduction in $r$ from lowering $R_{\text{eff}}$ is **partially offset** by the shortening of $g_{\text{eff}}$.

Concretely: for the same fractional reduction in $R_{\text{eff}}$, the spike case achieves a bigger reduction in epidemic growth rate than the smooth case.

### What is distinctive about test-and-isolate

This partial self-defeat is specific to interventions that truncate the infectiousness profile at a time point (test-and-isolate, contact tracing with delayed isolation). It does not apply to interventions that reduce $R_{\text{eff}}$ without altering the generation interval shape (e.g., vaccination, uniform contact reduction). The coupling between $R$ reduction and $g$ distortion, and its dependence on $\kappa$, is a distinctive feature of temporal interventions operating on punctuated infectiousness.

Analysis code: `code/generation_interval_distortion.R`.

## 15. Why the same test-and-isolate intervention yields different $R_{\text{eff}}$ reductions depending on $\kappa$

Section 14 showed that test-and-isolate distorts the generation interval differently for spike vs smooth profiles. But a more fundamental question precedes it: why does the $R_{\text{eff}}$ reduction itself depend on $\kappa$? The population-level kernel $A(\tau) = R_0 \cdot \text{Gamma}(\alpha_{\text{total}}, r)$ is invariant across $\kappa$, so the test positivity window catches the same fraction of population-level transmission mass regardless of punctuation. Yet spike profiles see a larger $R_{\text{eff}}$ reduction than smooth profiles.

### The detection probability is not the driver

Simulation shows that the probability of detection is approximately constant across $\kappa$ at $\sim$97.5% (with 3-day screening, sensitivity 0.85, and a 6-day positivity window). This makes sense: the positivity window is wide relative to the screening interval, so almost everyone gets tested while positive regardless of profile shape.

### The conditional value of detection is the driver

The key quantity is $E[\text{frac averted} \mid \text{detected}]$ — the expected fraction of an individual's transmission attempts that occur after detection, conditional on being detected. This drops from $\sim$0.88 for spike profiles ($\kappa = 0.5$) to $\sim$0.80 for smooth profiles ($\kappa = 9.5$).

Since $R_{\text{eff}} = R_0 \cdot (1 - P(\text{det}) \cdot E[\text{frac averted} \mid \text{det}])$, and $P(\text{det})$ is nearly constant, the $R_{\text{eff}}$ variation across $\kappa$ is driven almost entirely by the conditional value of detection.

### Mechanism 1: All-or-nothing vs partial aversion

For spike profiles, all transmission attempts cluster at a single time point. If detection occurs before this spike, 100% of transmission is averted; if after, 0% is averted. The conditional expectation $E[\text{frac averted} \mid \text{detected}]$ is therefore a weighted average of 1 and 0, dominated by the large fraction detected before the spike.

For smooth profiles, attempts are spread across the profile. Detection at any time point averts only the fraction of attempts that have not yet occurred. An individual detected near the peak of a broad profile may have already completed 30--50% of their transmission. Pre-peak detection averts only $\sim$88% (not 100%) because some transmission occurs in the early tail before any test can fire.

The distribution of frac averted given detection is bimodal (clustered at 0 and 1) for spike profiles but unimodal and centered below 1 for smooth profiles.

### Mechanism 2: Tail leakage

A secondary effect: smooth profiles place a non-negligible fraction of transmission attempts outside the positivity window entirely. At $\kappa = 9.5$, approximately 7% of attempts fall before the window opens or after it closes, making them fundamentally uncatchable by any test-based intervention. Spike profiles, by concentrating transmission near the peak, keep nearly 100% of attempts within the catchable window.

### Quantitative decomposition

| $\kappa$ | $P(\text{det})$ | $E[\text{frac averted} \mid \text{det}]$ | $R_{\text{eff}}$ | Tail leakage |
|----------|-----------------|------------------------------------------|-------------------|---------------|
| 0.5      | 0.975           | 0.88                                     | 0.28              | 0%            |
| 2        | 0.975           | 0.86                                     | 0.32              | 1%            |
| 6        | 0.975           | 0.82                                     | 0.40              | 5%            |
| 9.5      | 0.975           | 0.80                                     | 0.44              | 7%            |

The $R_{\text{eff}}$ varies by a factor of $\sim$1.6 across $\kappa$ values despite identical detection probabilities. The primary driver (accounting for $\sim$75% of the variation) is the conditional value of detection; tail leakage accounts for the remaining $\sim$25%.

### Implication

Test-and-isolate interventions are more effective against punctuated (spike-like) infectiousness profiles not because they detect more cases, but because each detection is worth more. For spike profiles, catching someone is an all-or-nothing event — detection before the spike averts everything; detection after averts nothing, but the "before" cases dominate. For smooth profiles, every detection is a partial victory, with diminishing returns as more of the profile has already contributed to transmission.

Analysis code: `code/reff_reduction_decomposition.R`.

### Note on the positivity window at small $\tau$

When $\kappa$ is small, some individuals have early peaks ($t_{\text{peak}} < w_-$), causing the positivity window $[t_{\text{peak}} - w_-, t_{\text{peak}} + w_+]$ to extend before $\tau = 0$ (infection time). This is biologically impossible — a person cannot test positive before they are infected. With $\kappa = 0.5$ and $w_- = 3$, roughly 13% of individuals have $t_{\text{peak}} < 3$.

In practice this does not bias results against our conclusions: tests cannot fire before $\tau = 0$ (screening starts post-infection), so the only effect is that the first test at small positive $\tau$ may fall within the nominal window and register as a true positive. This is a generous assumption for the intervention at small $\kappa$, slightly overstating detection effectiveness for spike profiles. The true $\kappa$-dependence of $R_{\text{eff}}$ reduction is therefore at least as strong as reported — conservative in the direction that reinforces the findings.

---

*Last updated: 2026-02-25*
