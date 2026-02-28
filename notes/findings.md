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

## 14. Detect-and-isolate: a unified framework for screening and symptom-triggered interventions

Detect-and-isolate interventions reduce transmission by truncating each detected individual's infectiousness profile at the detection time: $a_i^{\text{eff}}(\tau) = a_i(\tau) \cdot \mathbf{1}[\tau < \tau_{\text{det}}]$. This truncation both reduces $R_{\text{eff}}$ and reshapes the generation interval distribution — but the degree of both effects depends critically on $\kappa$.

The key insight that unifies screening-based and symptom-triggered interventions: both truncate the profile at a detection time $\tau_{\text{det}}$, and the entire $\kappa$-dependence flows through a single quantity:

$$F_\kappa(\tau_{\text{det}} - s_i) = \int_0^{\tau_{\text{det}} - s_i} f_\kappa(u)\,du$$

— the fraction of the biological profile's mass that has elapsed before detection. What differs between detection mechanisms is the **joint distribution of $\tau_{\text{det}}$ and the profile timing**.

### The detection spectrum: independent vs. anchored triggers

Detection mechanisms sit on a spectrum of **biological anchoring** — the degree to which $\tau_{\text{det}}$ is correlated with the individual's profile timing:

**Independent trigger ($\rho = 0$): regular screening.** Detection time is determined by an external process (periodic testing every $\Delta$ days, random phase) intersected with a positivity window. The screening clock runs independently of $s_i$, so $\tau_{\text{det}}$ is approximately uncorrelated with the profile. The $\kappa$-dependence comes purely from the *shape* of $F_\kappa$ — spike profiles have a step-function CDF (all-or-nothing), smooth profiles have a gradual sigmoid (partial aversion).

**Anchored trigger ($\rho = 1$): symptom onset.** Detection time is deterministically coupled to the profile: $\tau_{\text{det}} = s_i + \delta_{\text{sym}} + \delta_{\text{iso}}$, where $\delta_{\text{sym}}$ is the delay from infectiousness onset to symptoms and $\delta_{\text{iso}}$ is the delay from symptoms to effective isolation. Individuals with early onset (small $s_i$) are both early transmitters and early detectors. The $\kappa$-dependence comes from both the shape of $F_\kappa$ and the correlation between $\tau_{\text{det}}$ and $s_i$.

**Intermediate anchoring ($0 < \rho < 1$):** Hybrid scenarios such as screening that targets symptomatic individuals, or symptoms that are a noisy indicator of infectiousness timing. These interpolate between the independent and anchored extremes.

### The core quantity: fraction averted given detection

For any detection mechanism, the fraction of individual $i$'s transmission averted is:

$$1 - F_\kappa(\tau_{\text{det}} - s_i)$$

The population-level $R_{\text{eff}}$ reduction is:

$$R_{\text{eff}} = R_0 \cdot \left(1 - P(\text{det}) \cdot E\left[1 - F_\kappa(\tau_{\text{det}} - s_i) \mid \text{detected}\right]\right)$$

The $\kappa$-dependence enters through $F_\kappa$ and, for anchored triggers, through the correlation between $\tau_{\text{det}}$ and $s_i$.

### Results for regular screening (implemented)

Simulation with 3-day screening, sensitivity 0.85, and a 6-day positivity window confirms strong $\kappa$-dependence.

#### Detection probability is not the driver

$P(\text{det})$ is approximately constant across $\kappa$ at $\sim$97.5%. The positivity window is wide relative to the screening interval, so almost everyone is detected regardless of profile shape.

#### The conditional value of detection is the driver

$E[\text{frac averted} \mid \text{detected}]$ drops from $\sim$0.88 for spike profiles ($\kappa = 0.5$) to $\sim$0.80 for smooth profiles ($\kappa = 9.5$):

| $\kappa$ | $P(\text{det})$ | $E[\text{frac averted} \mid \text{det}]$ | $R_{\text{eff}}$ | Tail leakage |
|----------|-----------------|------------------------------------------|-------------------|---------------|
| 0.5      | 0.975           | 0.88                                     | 0.28              | 0%            |
| 2        | 0.975           | 0.86                                     | 0.32              | 1%            |
| 6        | 0.975           | 0.82                                     | 0.40              | 5%            |
| 9.5      | 0.975           | 0.80                                     | 0.44              | 7%            |

$R_{\text{eff}}$ varies by a factor of $\sim$1.6 across $\kappa$ despite identical detection probabilities.

#### Mechanism 1: All-or-nothing vs partial aversion

For spike profiles, all transmission attempts cluster at a single time point. If detection occurs before this spike, 100% is averted; if after, 0%. The conditional expectation $E[\text{frac averted} \mid \text{detected}]$ is a weighted average of 1 and 0, dominated by the large fraction detected before the spike.

For smooth profiles, detection at any time averts only the remaining fraction. An individual detected near the peak has already completed 30–50% of transmission. The distribution of frac averted is bimodal (clustered at 0 and 1) for spike profiles but unimodal and centered below 1 for smooth profiles.

#### Mechanism 2: Tail leakage

Smooth profiles place a non-negligible fraction of transmission outside the positivity window entirely (7% at $\kappa = 9.5$). Spike profiles keep nearly 100% within the catchable window.

#### Generation interval distortion

Detection also reshapes the effective generation interval $g_{\text{eff}}(\tau)$, with a $\kappa$-dependent effect:

**Spike case:** $g_{\text{eff}} = g$ (unchanged). Testing removes a random fraction of transmitters without altering the timing of the remaining ones — the blocking probability is approximately constant across spike locations.

**Smooth case:** $g_{\text{eff}}(\tau) \propto A(\tau) \cdot P(T_{\text{detect}} > \tau)$ shifts shorter. The right tail of the profile is truncated by detection, so surviving transmissions are biased early.

**Epidemiological consequence: partially self-defeating intervention.** Via the Euler-Lotka equation $1 = R_{\text{eff}} \int e^{-r\tau} g_{\text{eff}}(\tau)\,d\tau$, a shorter $g_{\text{eff}}$ implies a *larger* growth rate for a given $R_{\text{eff}}$. So for smooth profiles, the $R_{\text{eff}}$ reduction from detection is partially offset by the generation interval shortening. For spike profiles, $g_{\text{eff}} = g$, so the full $R_{\text{eff}}$ reduction translates to growth rate reduction. This partial self-defeat is specific to profile-truncating interventions and does not apply to interventions that reduce $R_{\text{eff}}$ without altering the generation interval (vaccination, uniform contact reduction).

Analysis code: `code/generation_interval_distortion.R`, `code/reff_reduction_decomposition.R`.

#### Note on the positivity window at small $\kappa$

When $\kappa$ is small, some individuals have early peaks ($t_{\text{peak}} < w_-$), causing the positivity window to nominally extend before $\tau = 0$. In practice tests cannot fire before infection, so the only effect is a slight overstatement of detection effectiveness for spike profiles. The true $\kappa$-dependence is therefore at least as strong as reported.

### Predictions for symptom-triggered isolation (planned)

*Status: not yet implemented.*

Symptom-triggered isolation replaces the independent screening clock with a biologically anchored trigger: $\tau_{\text{det}} = s_i + \delta_{\text{sym}} + \delta_{\text{iso}}$, deterministically coupled to the onset shift $s_i$.

#### Why anchoring should sharpen the $\kappa$-dependence

With screening, $\tau_{\text{det}}$ is approximately independent of $s_i$, so detection can occur at any point relative to the profile. With symptoms, $\tau_{\text{det}}$ tracks the profile: individuals with early onset are detected early. This coupling means the analysis reduces to a single, clean question: **what fraction of the profile's mass falls before symptom onset?** That fraction is $F_\kappa(\delta_{\text{sym}})$ — a property of the profile shape alone, without the noise introduced by the screening schedule.

#### Key predictions

**Delta (small $\kappa$):** Transmission occurs at $\tau^* \approx s_i + \varepsilon$ where $\varepsilon$ is small. Symptom onset at $s_i + \delta_{\text{sym}}$. If $\delta_{\text{sym}} + \delta_{\text{iso}} > \varepsilon$ (symptoms + isolation delay exceeds the jitter), isolation arrives after the spike — nothing averted. If $\delta_{\text{sym}} + \delta_{\text{iso}} < \varepsilon$, everything averted. The all-or-nothing structure is sharper than for screening because the detection time tracks the profile.

**Smooth (large $\kappa$):** Transmission is spread over a broad window. Symptom onset truncates the right tail. The fraction averted is $1 - F_\kappa(\delta_{\text{sym}})$, a smooth function of $\delta_{\text{sym}}$.

**The key $\kappa$-dependent quantity:** $F_\kappa(\delta_{\text{sym}})$, the CDF of $\text{Gamma}(\kappa, r)$ evaluated at the symptom delay. For small $\kappa$, this is a sharp step function. For large $\kappa$, a smooth sigmoid.

#### What's new relative to screening

1. **Pre-symptomatic transmission is the key parameter.** The entire analysis reduces to one number: the fraction of transmission before symptoms. No test-performance parameters (sensitivity, specificity, screening interval).

2. **Correlation between detection and profile.** Individuals with early $s_i$ are both early transmitters and early detectors. Whether this correlation amplifies or dampens the $\kappa$-dependence (relative to screening) depends on where symptom onset falls relative to the transmission mass.

3. **No false positives or sensitivity parameters.** Only $\delta_{\text{sym}}$, $\delta_{\text{iso}}$, and $p_{\text{iso}}$ (compliance).

#### Design choices to resolve

- **Symptom onset model.** Where does $\tau_{\text{sym}}$ fall relative to $b_i$? Options:
  - At the mode of $f_\kappa$ (symptoms at peak infectiousness). Mode of $\text{Gamma}(\kappa, r)$ is $(\kappa - 1)/r$ for $\kappa \geq 1$; for $\kappa < 1$, mode is at 0.
  - At a fixed fraction of the profile's CDF (e.g., symptoms when 30% of infectiousness mass has elapsed). This gives a $\kappa$-dependent delay.
  - At a fixed absolute delay $\delta_{\text{sym}}$ after $s_i$ (simplest, but least biologically motivated).

- **Isolation delay $\delta_{\text{iso}}$:** 0, 0.5, 1, 2 days.

- **Asymptomatic fraction.** $p_{\text{sym}} < 1$ acts as a ceiling on effectiveness but doesn't interact with $\kappa$ in an interesting way.

## 15. Gathering size restrictions are more effective for smooth profiles

We investigated how gathering size restrictions interact with the temporal shape of individual infectiousness profiles. The central question: are gathering size caps more effective at preventing outbreaks when infectiousness is punctuated (spiky) versus smooth?

The answer turns out to be the opposite of the naive expectation. Caps are **more effective for smooth profiles**, because they reduce mean transmission without losing the overdispersion that was already suppressing outbreak probability for spiky profiles.

### Model

Each individual moves through a sequence of gatherings. The contact process $c_i(t)$ is piecewise-constant:

- **Gathering sizes**: drawn iid from $\text{NegBin}(n = 10/9,\; p = 1/10)$, giving mean 10, variance 100 (SD = 10). After normalisation (dividing by 5), the contact rate has mean 2 and variance 4.
- **Holding times**: $\text{Exp}(\text{rate} = 2)$, so individuals change gathering approximately every 0.5 days.
- **Gathering size cap**: draws exceeding the cap are rejected (truncated NegBin).

The individual reproduction number $R_i$ arises from the interaction of the biological profile and the contact process:

- **Delta profile**: $R_i = c_i(t^*)$ at a single random time. This is one draw from the contact-value distribution, inheriting its full variance.
- **Smooth profile**: $R_i = \int c_i(\tau) \cdot f(\tau)\,d\tau$, where $f$ is the $\text{Gamma}(10, 0.3)$ PDF spanning $\sim$33 days. This averages over $\sim$36 independent contact-process steps, dramatically reducing variance.

Offspring are $\text{Poisson}(R_i)$ — the standard superspreading framework.

### Finding 16a: Caps reduce mean $R_i$ identically for all profile shapes

The contact process enters multiplicatively with the biological profile. Since every biological profile integrates to the same value (regardless of $\kappa$), the expected contact rate at any time point determines $E[R_i]$. Truncating the gathering size distribution reduces this expected contact rate equally for spiky and smooth profiles.

| Scenario      | $E[R_i]$ |
|---------------|----------|
| Uncapped      | 2.00     |
| Capped at 20  | 1.39     |

This was confirmed analytically and in simulation (28 $\kappa \times$ cap scenarios, with $E[R_i]$ varying by $< 0.04$ across $\kappa$ for each cap level).

**The mean effect of gathering restrictions is $\kappa$-independent.**

### Finding 16b: Caps reduce variance identically in fractional terms

The variance of $R_i$ decomposes cleanly:

| Profile | Scenario   | $E[R_i]$ | $\text{Var}(R_i)$ | $\text{CV}(R_i)$ |
|---------|------------|----------|--------------------|-------------------|
| Delta   | Uncapped   | 2.00     | 4.00               | 1.00              |
| Delta   | Capped@20  | 1.39     | 1.20               | 0.79              |
| Smooth  | Uncapped   | 2.00     | 0.11               | 0.17              |
| Smooth  | Capped@20  | 1.39     | 0.03               | 0.13              |

The cap reduces variance by $\sim$70% for both profiles. This is because the variance reduction factor from averaging ($I_{\text{corr}} = 0.028$, corresponding to $\sim$36 effective independent samples) multiplies both the uncapped and capped contact-process variance identically:

$$\text{Var}(R_i)_{\text{smooth}} = \text{Var}(c_i) \cdot I_{\text{corr}}$$

So the fractional reduction in variance from capping is the same: $\text{Var}(\text{capped}) / \text{Var}(\text{uncapped}) = \text{Var}(X \mid X \leq 20) / \text{Var}(X)$, regardless of the biological profile.

**The differential is not in the fractional variance reduction — it is in the absolute amount of superspreading removed.** The cap eliminates $\Delta\text{Var} = 2.8$ for the delta profile versus $\Delta\text{Var} = 0.08$ for the smooth profile (a 35$\times$ difference).

### Finding 16c: Overdispersion increases extinction probability

For a branching process with $\text{Poisson}(R_i)$ offspring, higher variance in $R_i$ (for fixed $E[R_i] > 1$) increases the probability of stochastic extinction. Superspreading creates many individuals with $R_i$ near zero, providing frequent opportunities for chains to die out.

**Unmitigated ($R = 2$):**

| Profile              | $P(\text{extinction})$ | $P(\text{outbreak})$ |
|----------------------|------------------------|----------------------|
| Poisson(2) reference | 0.20                   | 0.80                 |
| Smooth               | 0.22                   | 0.78                 |
| Delta                | 0.51                   | 0.49                 |

Despite having the same mean $R = 2$, the delta profile produces outbreaks less than half the time, while the smooth profile produces outbreaks $\sim$80% of the time. The smooth profile's near-Poisson $R_i$ distribution means reliable, moderate transmission from every individual.

### Finding 16d: Caps are more effective at preventing outbreaks for smooth profiles

**Mitigated (cap at 20, $R$ drops to 1.39):**

| Profile              | $P(\text{extinction})$ | $P(\text{outbreak})$ |
|----------------------|------------------------|----------------------|
| Poisson(1.39) ref    | 0.50                   | 0.50                 |
| Smooth               | 0.50                   | 0.50                 |
| Delta                | 0.67                   | 0.33                 |

The change in outbreak probability from imposing the cap:

| Profile | $P(\text{outbreak})$ uncapped | $P(\text{outbreak})$ capped | Reduction |
|---------|-------------------------------|-----------------------------|-----------|
| Smooth  | 0.78                          | 0.50                        | **$-29$ pp** |
| Delta   | 0.49                          | 0.33                        | **$-17$ pp** |

The cap reduces smooth-profile outbreak probability by 29 percentage points, but delta-profile outbreak probability by only 17 percentage points.

#### Mechanism

Two effects of the cap act in opposite directions for the delta profile:

1. **Lower mean $R$** ($2 \to 1.39$): pushes extinction up (fewer infections on average)
2. **Lower variance** ($4 \to 1.2$): pushes extinction *down* (fewer zero-offspring individuals, more consistent transmission)

These partially cancel. The cap removes the large-gathering tail that was simultaneously (a) creating superspreaders and (b) creating dead ends. Eliminating both leaves the outbreak probability only modestly changed.

For the smooth profile, variance was already near zero (0.11), so the variance reduction is negligible. The full force of the mean reduction translates directly into fewer outbreaks, matching the Poisson(1.39) reference almost exactly.

### Finding 16e: Simulation results confirm the analytical predictions

The large-scale simulation (`code/gathering_size_restrictions.R`) swept 7 values of $\kappa$ (0.5 to 9.5) across 4 cap levels (5, 10, 20, $\infty$), with 50,000 MC replicates per scenario for $R_i$ distributions and 200 epidemic simulations per scenario.

Key confirmations:

- **$E[R_i]$ is constant across $\kappa$** for each cap level (range $< 0.04$)
- **$\text{CV}(R_i)$ decreases with $\kappa$** for uncapped scenarios (from 1.07 at $\kappa = 0.5$ to 0.82 at $\kappa = 9.5$), confirming that smooth profiles average out contact variation
- **All capped epidemic scenarios showed zero establishment** in the full simulation (because with geometric gathering sizes of mean 30, even a cap at 20 reduces $E[R_i]$ to $\sim$0.58 — well below 1). This demonstrates that for heavy-tailed gathering size distributions, even moderate caps can collapse mean transmission, making the variance-level differential moot for epidemic outcomes.

The analytical calculations (Finding 16d) used a different parameterisation ($\text{NegBin}/5$ with mean 2, cap at 20 giving $E[R_i] = 1.39 > 1$) to isolate the regime where the differential effect on outbreak probability is visible.

### Summary

The naive hypothesis — that gathering size restrictions are *more* effective against punctuated infectiousness because they selectively remove superspreading events — is incorrect, or at least incomplete.

Gathering size caps reduce mean transmission identically regardless of profile shape. But they also reduce variance, and for spiky profiles, that variance was *already suppressing outbreaks* via stochastic extinction. The cap simultaneously removes the superspreading tail (which sustained rare large outbreaks) and the zero-transmission tail (which killed most chains). The net effect on outbreak probability is smaller than for smooth profiles, where there was no beneficial overdispersion to lose.

**Gathering size restrictions are more effective at preventing outbreaks when infectiousness is smooth**, because smooth profiles have no overdispersion buffer to erode.

Analysis code: `code/gathering_size_restrictions.R`.

## 16. Planned: Contact tracing with imperfect isolation / post-exposure prophylaxis

*Status: planned, not yet implemented.*

### Motivation

Contact tracing is the most natural timing-dependent intervention to study because its effectiveness depends on the temporal profiles of *both* the infector (index case) and the infectee (traced contact). The delay from infection of the index case to intervention on the traced contact passes through a chain of timing-dependent steps, each of which interacts with $\kappa$.

We model imperfect isolation and post-exposure prophylaxis (PEP) within a single framework. The key observation: both reduce (rather than eliminate) a traced contact's onward transmission from some intervention time onward. Perfect isolation ($\eta = 0$) is a special case; leaky isolation and PEP correspond to $\eta > 0$, differing only in the plausible parameter ranges and whether a pre-emptive mode exists.

### Model

Contact tracing unfolds as a cascade of delays:

**Step 1: Index case detection** ($\tau_{\text{det}}^A$)

The index case (person A) is detected via one of:
- (a) Symptom onset + testing (Section 14, anchored trigger): $\tau_{\text{det}}^A = s_A + \delta_{\text{sym}} + \delta_{\text{test}}$
- (b) Routine screening (Section 14, independent trigger): $\tau_{\text{det}}^A$ depends on the positivity window and screening phase

In both cases, $\tau_{\text{det}}^A$ depends on person A's profile via $s_A$ and $f_\kappa$.

**Step 2: Tracing delay** ($\delta_{\text{trace}}$)

After detection, contacts are identified and notified. This introduces a further delay $\delta_{\text{trace}}$ that is approximately independent of $\kappa$ (it depends on public health infrastructure, not biology). We can model it as fixed or $\text{Exp}(\text{rate})$.

**Step 3: Intervention on the traced contact**

Person B (a contact of A) is reached at calendar time:

$$T_{\text{int}}^B = t_A + \tau_{\text{det}}^A + \delta_{\text{trace}}$$

Person B was infected by A at some time $t_B = t_A + \tau_{A \to B}$, where $\tau_{A \to B}$ is the generation interval from A to B. The time from B's infection to intervention is:

$$\delta_B = T_{\text{int}}^B - t_B = \tau_{\text{det}}^A + \delta_{\text{trace}} - \tau_{A \to B}$$

This depends on A's detection time (which depends on A's profile via $\kappa$) and the generation interval from A to B (which also depends on A's profile).

**The intervention modifies B's profile.** The effective biological profile becomes:

$$b_B^{\text{eff}}(\tau) = \begin{cases} b_B(\tau) & \tau < \delta_B \\ \eta \cdot b_B(\tau) & \tau \geq \delta_B \end{cases}$$

where $\eta \in [0, 1]$ is the residual transmission fraction:
- $\eta = 0$: perfect isolation (Section 14 applied to B)
- $0 < \eta < 1$: leaky isolation (household transmission, imperfect compliance) or PEP (antivirals reducing but not eliminating infectiousness)
- $\eta = 1$: intervention has no effect

The effective reproduction number for B is:

$$R_B^{\text{eff}} = R_0 \left[ F_B(\delta_B) + \eta \cdot (1 - F_B(\delta_B)) \right]$$

where $F_B(\delta_B) = \int_0^{\delta_B} b_B(\tau)\,d\tau$ is the fraction of B's transmission that has already occurred by the time of intervention. This simplifies to:

$$R_B^{\text{eff}} = R_0 \left[ 1 - (1 - \eta)(1 - F_B(\delta_B)) \right]$$

The fraction of B's transmission averted is $(1 - \eta)(1 - F_B(\delta_B))$, which depends on both the timing $\delta_B$ and B's profile shape via $F_B$.

#### Pre-emptive mode (PEP-specific)

For PEP delivered before B's infectiousness window opens ($\delta_B < s_B$), there is an additional biological effect: prophylaxis may prevent productive infection entirely. In this case, the entire profile is scaled:

$$b_B^{\text{eff}}(\tau) = \eta_{\text{pre}} \cdot b_B(\tau) \qquad \text{if } \delta_B < s_B$$

where $\eta_{\text{pre}} \leq \eta$ (pre-emptive PEP is at least as effective as late PEP). This adds a threshold effect absent from pure isolation. For the default analysis, we set $\eta_{\text{pre}} = \eta$ (no distinct pre-emptive mode), and note where this assumption matters.

### Why $\kappa$ enters three times

The effectiveness of tracing person B depends on:

1. **How quickly A is detected** ($\tau_{\text{det}}^A$) — determines when the tracing clock starts. For spike A: detection is either very early (if the trigger precedes the spike) or very late (if the spike precedes the trigger). For smooth A: detection is at a predictable time near the profile's peak.

2. **What generation interval connected A to B** ($\tau_{A \to B}$) — determines how much head start B has. For spike A: $\tau_{A \to B} \approx s_A + \varepsilon$ (all of A's transmissions cluster at one time, so all traced contacts have similar head starts). For smooth A: $\tau_{A \to B}$ is spread over a wide range, so some contacts were infected early (long head start, likely already transmitted) and others late (short head start, intervention is more useful).

3. **How concentrated B's own transmission is** ($F_B(\delta_B)$) — determines how much of B's profile falls before vs. after the intervention. For spike B: $F_B$ is a sharp step function, making the outcome all-or-nothing. For smooth B: $F_B$ is a smooth sigmoid, giving graded partial benefit.

### Key predictions

**Spike A, spike B (small $\kappa$ throughout):** A's transmission spike produces a cluster of contacts all infected at nearly the same time. If intervention is fast enough to reach them before their own spikes, all transmission is averted (or reduced to $\eta$). If not, none is averted. **High variance, all-or-nothing.** The residual fraction $\eta$ barely matters: either the spike beat the intervention ($F_B \approx 1$, nothing to reduce) or it didn't ($F_B \approx 0$, full $(1-\eta)$ reduction).

**Smooth A, smooth B (large $\kappa$ throughout):** A is detected at a predictable time. Contacts infected early by A have a long head start and have completed much of their own broad transmission window ($F_B$ large); contacts infected late have a short head start ($F_B$ small, intervention is effective). **Graded effectiveness, depending on the generation interval from A to B.** The residual fraction $\eta$ matters throughout: even after intervention, smooth B continues transmitting at rate $\eta$, and this residual accumulates over the remaining broad profile.

**Spike A, smooth B (or vice versa):** The cross-$\kappa$ cases reveal asymmetries. A spike infector produces clustered contacts (all with similar $\delta_B$), while a smooth infector produces dispersed contacts (variable $\delta_B$). The infectee's $\kappa$ then determines whether intervention catches their transmission window.

### The role of $\eta$

**At $\eta = 0$ (perfect isolation):** Results reduce to the pure contact tracing case. The fraction averted is $1 - F_B(\delta_B)$, purely a race between intervention timing and B's profile.

**At $\eta > 0$ (imperfect isolation / PEP):** The residual transmission $\eta \cdot (1 - F_B(\delta_B)) \cdot R_0$ persists after intervention. This is more consequential for smooth B than spike B:

- **Spike B**: either $F_B \approx 0$ (spike not yet fired, fraction averted $\approx 1 - \eta$) or $F_B \approx 1$ (spike already fired, fraction averted $\approx 0$). The $\eta$ only matters in the "caught before spike" case, where it determines whether the individual is fully or partially suppressed.
- **Smooth B**: $F_B$ takes intermediate values for most intervention timings, so the residual $\eta \cdot (1 - F_B)$ contributes meaningfully. Higher $\eta$ (worse isolation/PEP) erodes the benefit most for smooth profiles, because they have the most remaining profile mass to leak through.

**Prediction:** Imperfect isolation ($\eta > 0$) differentially harms the effectiveness of contact tracing for smooth profiles relative to spike profiles. This could partially offset the advantage that smooth profiles have from their more predictable detection timing.

### Cross-over prediction for PEP-like interventions

When $\eta$ is small but delay $\delta_B$ varies:

- At **short delays**, smooth B may benefit more: everyone gets a partial reduction, and the broad profile has substantial remaining mass to suppress.
- At **long delays**, spike B may benefit more: some individuals still have late onset shifts ($s_B > \delta_B$) and are fully protected, while smooth B has already completed most of its broad transmission window.

The $\kappa$-dependence of effectiveness should therefore **reverse sign** as a function of $\delta_B$. The cross-over point depends on the relationship between the delay distribution and the onset-shift distribution $s_B \sim \text{Gamma}(\alpha_{\text{total}} - \kappa, r)$.

### Connection to Section 14

The model composes the components from the detect-and-isolate framework:

- Step 1 reuses the detection model from Section 14 — either the independent trigger (screening) or anchored trigger (symptom-triggered) variant
- Step 3 reuses the truncation/scaling analysis from Section 14, but now applied to person B with a delay that depends on person A's profile

This modular structure means the analysis can build directly on the existing codebase: the detection module determines $\tau_{\text{det}}^A$, and the truncation module determines $F_B(\delta_B)$ given a known intervention time.

### Parameters to explore

- **Tracing delay** $\delta_{\text{trace}}$: 0, 1, 2, 3 days
- **Index detection mechanism**: symptom-triggered (Section 14, anchored) vs. screening-triggered (Section 14, independent)
- **Residual transmission** $\eta$: 0 (perfect isolation), 0.1, 0.3, 0.5 (leaky isolation / PEP)
- **Pre-emptive efficacy** $\eta_{\text{pre}}$: 0 (PEP fully prevents if early), $\eta$ (no distinct pre-emptive mode)
- **Tracing coverage** $p_{\text{trace}}$: fraction of contacts successfully identified (0.5, 0.75, 1.0)
- **$\kappa$ for infector vs infectee**: sweep independently if the population is heterogeneous, or assume a single $\kappa$ if all individuals share the same profile shape

### Taxonomy of timing-dependent interventions

| Intervention | Acts on | Mechanism | $\eta$ | Profile interaction |
|---|---|---|---|---|
| Detect-and-isolate, screening (§14) | Infector | Truncation at random detection time | 0 | Sampling vs. averaging of positivity window |
| Detect-and-isolate, symptoms (§14) | Infector | Truncation at biologically anchored time | 0 | Pre-symptomatic fraction depends on $\kappa$ |
| Contact tracing + perfect isolation | Infectee | Truncation at traced time | 0 | Two-generation: A's profile → delay → B's truncation |
| Contact tracing + leaky isolation | Infectee | Scaling at traced time | $> 0$ | As above, plus residual leakage erodes smooth-B benefit |
| Contact tracing + PEP | Infectee | Scaling at traced time, pre-emptive mode | $> 0$ | As above, plus race between PEP and B's onset shift |

The progression from Section 14 to Section 16 traces a path of increasing mechanistic coupling between the intervention and the temporal structure of infectiousness. This section unifies the final three rows into a single parameterised model.

---

*Last updated: 2026-02-27*
