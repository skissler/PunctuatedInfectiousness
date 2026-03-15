## Spike case takes off slower: a "minimum of order statistics" mechanism

Conditional on establishing, the spike case takes longer to reach early case-count thresholds than the smooth case. 

In the smooth case, an infector's $\chi_i$ infection attempts have times drawn from the generation interval distribution, $g(\tau)$. The *first* onward transmission is $\min(\tau_1, \ldots, \tau_{\chi_i})$, which is stochastically earlier than a single draw from $g(\tau)$. In the spike case, all $\chi_i$ infection attempts occur at the same time $\tau_i^*$, so the first transmission time is just a single generation-interval draw regardless of $\chi_i$.

### R0 dependence

The growth delay from punctuation **increases** with $R_0$. This is because the order-statistics advantage of smooth profiles scales with the expected number of infection attempts per infector, $E[k] = R_0$. At $R_0 = 1.2$, most infectors draw only $\chi_i = 1$ attempt, so there is no order-statistics effect to exploit — the two profiles produce indistinguishable growth. At $R_0 = 5$, infectors typically draw $k = 4 \text{-} 6$ attempts, and $\min(\tau_1, \ldots, \tau_5)$ is substantially earlier than a single $\tau_j$, giving smooth profiles a compounding head start at each generation. The result is that punctuation slows early growth most precisely when $R_0$ is large, i.e., when the epidemic would otherwise grow fastest.

## Spike epidemics have more variable empirical exponential growth rates

Estimates of the early-epidemic exponential growth rate $r$ have early identical means for smooth and spike infectiousness profiles, but the estimates for the spike case have higher standard deviation, based on log-linear fit to cases 10–100. 

## Overdispersion from punctuated infectiousness with periodic contacts

We decompose individual infectiousness as $a_i(\tau) = b_i(\tau) \cdot R_0 z(t_i + \tau)$, where:

- $b_i(\tau) = f_\psi(\tau - l_i)$ is the biological timing density (a Gamma($\psi \alpha$, $\beta$) density shifted by individual onset $`l_i \sim \text{Gamma}((1-\psi) \alpha, \beta)`$)
- $z(t) = 1 + \epsilon \sin(2\pi t / T)$ is a periodic contact rate multiplier with period $T = 7$ days

The parameter $\psi$ controls punctuation: small $\psi$ gives narrow, spike-like individual profiles; large $\psi$ gives broad profiles that resemble the population average. The parameter $\epsilon$ controls the amplitude of contact variation ($\epsilon = 0$ is constant contacts, $\epsilon = 0.9$ is strong weekly oscillation).

The population-level kernel $A(\tau) = R_0 \cdot \text{Gamma}(\tau; \alpha, \beta)$ is invariant across $\psi$ by the Gamma additivity property. Parameters: $R_0 = 2$, $\alpha = 10$, $\mu = 5$ days, so $r = 2$.

### Part A: Branching process offspring distributions

For each ($\psi$, $\epsilon$) pair, we drew 10,000 individuals with random infection phase $t_i \sim \text{Uniform}(0, T)$ and computed the realized reproduction number $\nu_i = R_0 \int b_i(\tau) \, z(t_i + \tau) \, d\tau$ via numerical integration, then drew offspring $n_i \sim \text{Poisson}(\nu_i)$.

#### Verification: mean offspring

The mean offspring count is $\approx 2.0$ across all ($\psi$, $\epsilon$) combinations (range: 1.91--2.02). This confirms that periodic contacts redistribute infectiousness across individuals but do not change the mean. The modulation is purely a second-order effect.

#### Verification: Poisson baseline

At $\epsilon = 0$ (constant contacts), the variance/mean ratio is $\approx 1.0$ for all $\psi$ values (range: 0.98--1.01). This confirms that without contact variation, the offspring distribution is Poisson regardless of how punctuated the biological profile is. Punctuation alone does not create overdispersion in a branching process with constant contacts.

#### Variance/mean ratio

The core finding: overdispersion arises from the *interaction* between punctuated infectiousness and periodic contacts. Neither ingredient alone is sufficient.

| $\psi$ | $\epsilon = 0$ | $\epsilon = 0.3$ | $\epsilon = 0.5$ | $\epsilon = 0.7$ | $\epsilon = 0.9$ |
|----------|----------------|-------------------|-------------------|-------------------|-------------------|
| 0.5      | 1.00           | 1.10              | 1.20              | 1.41              | 1.69              |
| 1        | 0.98           | 1.08              | 1.19              | 1.41              | 1.71              |
| 2        | 0.99           | 1.05              | 1.21              | 1.37              | 1.60              |
| 4        | 0.98           | 1.00              | 1.12              | 1.27              | 1.42              |
| 6        | 1.00           | 1.03              | 1.06              | 1.17              | 1.24              |
| 8        | 1.00           | 1.01              | 1.06              | 1.12              | 1.16              |
| 9.5      | 1.01           | 1.05              | 1.04              | 1.07              | 1.12              |

The overdispersion is concentrated in the low-$\psi$, high-$\epsilon$ corner: punctuated profiles ($\psi \leq 2$) with strong contact oscillation ($\epsilon \geq 0.7$) yield variance/mean ratios of 1.4--1.7. Smooth profiles ($\psi = 9.5$) show at most a ratio of 1.12 even at maximum contact amplitude.

#### Negative binomial dispersion parameter $k$

Fitting $k$ via method of moments ($k = \mu^2 / (\text{Var} - \mu)$):

| $\psi$ | $\epsilon = 0$ | $\epsilon = 0.5$ | $\epsilon = 0.7$ | $\epsilon = 0.9$ |
|----------|----------------|-------------------|-------------------|-------------------|
| 0.5      | $\infty$       | 9.7               | 4.7               | 2.8               |
| 1        | $\infty$       | 10.3              | 4.9               | 2.8               |
| 2        | $\infty$       | 9.6               | 5.4               | 3.4               |
| 4        | $\infty$       | 17.2              | 7.3               | 4.7               |
| 6        | $\infty$       | 30.9              | 11.8              | 8.4               |
| 8        | $\infty$       | 32.5              | 16.5              | 12.4              |
| 9.5      | $\infty$       | 51.0              | 27.7              | 16.8              |

At $\epsilon = 0$, all $\psi$ values give $k = \infty$ (Poisson). At $\epsilon = 0.9$ with $\psi = 0.5$, we get $k \approx 2.8$ -- substantial overdispersion from a purely mechanistic origin, comparable to values estimated for some respiratory infections. The gradient from $k \approx 17$ (smooth + periodic) to $k \approx 3$ (punctuated + periodic) demonstrates that punctuation amplifies contact-driven overdispersion by roughly a factor of 6 in the dispersion parameter.

#### Mechanism

The mechanism is transparent: when the biological profile $b_i$ is a narrow spike (small $\psi$), the integral $\int b_i(\tau) \, z(t_i + \tau) \, d\tau$ essentially samples $z$ at one point, so $\nu_i$ varies between $R_0(1 - \epsilon)$ and $R_0(1 + \epsilon)$ depending on when the spike falls relative to the contact cycle. When $b_i$ is broad (large $\psi$), it averages over the contact cycle, and $\nu_i \approx R_0$ for everyone.

This is a sampling-vs-averaging effect: punctuated profiles *sample* the contact function; smooth profiles *average* over it.

### Part B: Figures

Generated in `figures/`:

- **`fig_overdispersion_heatmap.pdf`**: Var/mean ratio across the ($\psi$, $\epsilon$) grid. Overdispersion is concentrated in the low-$\psi$, high-$\epsilon$ corner, confirming the interaction.
- **`fig_overdispersion_lines.pdf`**: Var/mean ratio vs $\psi$ for selected $\epsilon$ values. All curves converge to $\approx 1$ at high $\psi$, and fan out at low $\psi$ proportional to $\epsilon$.
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

1. **Punctuation alone does not create overdispersion** in a branching process with constant contacts. The offspring distribution is Poisson($R_0$) for all $\psi$ when $\epsilon = 0$.

2. **Contact variation alone creates mild overdispersion** even with smooth profiles. At $\epsilon = 0.9$ and $\psi = 9.5$, the variance/mean ratio is 1.12 ($k \approx 17$).

3. **The interaction amplifies overdispersion**: punctuated profiles + periodic contacts yield variance/mean $\approx 1.7$ and $k \approx 3$ -- a 6-fold reduction in the dispersion parameter compared to smooth profiles with the same contact variation.

4. **The mechanism is sampling vs averaging**: narrow individual profiles sample the contact function at a random phase; broad profiles average over it. This is a direct consequence of the decomposition $a_i(\tau) = R_0 \cdot b_i(\tau) \cdot z(t_i + \tau)$.

5. **Epidemic consequences**: the overdispersion reduces establishment probability (0.55 vs 0.75) and increases variability in epidemic outcomes, without changing the mean offspring count.

6. **Punctuation slows early growth, amplified by $R_0$**: spike profiles reach 100 infections 3--18% slower than smooth profiles (at $R_0 = 2 \text{-} 5$), via a minimum-of-order-statistics mechanism. The effect grows with $R_0$ because more infection attempts per infector means a larger order-statistics advantage for smooth profiles. At $R_0 \leq 1.2$, the effect vanishes because most infectors only generate one attempt.


## 17. Why punctuated profiles yield slower early growth despite identical mean-field dynamics

### The apparent paradox

Findings 3, 4, 5, and 12 all show that spike (punctuated) profiles produce slower and more variable early epidemic growth than smooth profiles. Yet by construction, both profiles share the same population-level infectiousness kernel $A(\tau)$, and therefore the same deterministic renewal equation, the same Euler–Lotka growth rate $\alpha$, and the same mean-field epidemic trajectory. How can two models with identical $E[Z(t)]$ consistently differ in the time to reach epidemic milestones?

The resolution involves two distinct insights:

1. "Same expected trajectory" does not imply "same expected time to reach a milestone." The two quantities are related by a nonlinear transformation, and Jensen's inequality drives a wedge between them.

2. More fundamentally: the **mean** trajectory is not the **typical** trajectory. The spike case's $Z(t)$ distribution is more right-skewed than the smooth case's — rare explosive realizations pull up $E[Z(t)]$ to match, even though the median and most individual realizations are slower.

### Setup: the Crump–Mode–Jagers branching process

The early epidemic (before susceptible depletion matters) is a Crump–Mode–Jagers (CMJ) branching process. Let $Z(t)$ denote the total number of individuals infected by time $t$ (cumulative incidence), starting from a single index case at $t = 0$. Each individual $i$, infected at time $t_i$, independently generates offspring according to a point process on $[0, \infty)$ with intensity measure $\mu(d\tau) = A(\tau)\,d\tau$.

For the **smooth** model, individual $i$ draws $k_i \sim \text{Poisson}(R_0)$ offspring at independent times $\tau_1, \ldots, \tau_{k_i} \sim g(\tau)$.

For the **spike** model, individual $i$ draws $k_i \sim \text{Poisson}(R_0)$ offspring, all at the same time $\tau_i^* \sim g(\tau)$.

Both models have the same intensity measure:

$$\mu(B) = E[\text{number of offspring in time set } B] = R_0 \int_B g(\tau)\,d\tau = \int_B A(\tau)\,d\tau$$

for any measurable set $B$. This is because $E[k_i \cdot \mathbf{1}(\tau_i^* \in B)] = R_0 \cdot P(\tau^* \in B) = R_0 \int_B g(\tau)\,d\tau$ for the spike case, which matches the smooth case. Thus:

$$E[Z(t)] \text{ is identical for both models at every } t.$$

The Malthusian growth rate $\alpha$, defined by the Euler–Lotka equation

$$1 = \int_0^\infty e^{-\alpha \tau} A(\tau)\,d\tau,$$

is the same. The mean-field dynamics are truly, provably, identical.

**Importantly, this is not a large-numbers approximation.** The proof uses the Campbell theorem for point processes: $E[\sum_j f(\tau_j)] = \int f(s)\,\mu(ds)$ for any point process with intensity measure $\mu$, regardless of the dependence structure between the points $\tau_j$. It does not matter whether offspring times are independent (smooth) or perfectly correlated (spike). The renewal equation for $E[Z(t)]$ inherits this, giving an identical expected cumulative incidence curve at every $t$, even starting from a single case.

### The mean trajectory is not the typical trajectory

If this is true, why do simulations consistently show the spike case reaching milestones later? Because $E[Z(t)]$ — the average across many simulations at each fixed $t$ — is not what we see when we look at individual simulation curves. What we see is closer to the **median** of $Z(t)$, which differs between models even when the mean does not.

The spike case's $Z(t)$ has a more right-skewed distribution than the smooth case's. A few explosive realizations — where early spikes all happen to fall at short generation intervals — produce very large $Z(t)$ values that pull up $E[Z(t)]$ to match the smooth case. But these are rare. The typical (median) realization is slower.

**A concrete toy example.** Suppose at $t = 5$ days, five simulations of each model give:

- Smooth: $Z(5) = 8, 9, 10, 11, 12$. Mean $= 10$. Median $= 10$. All near threshold.
- Spike: $Z(5) = 3, 5, 6, 8, 28$. Mean $= 10$. Median $= 6$. Only 1 of 5 has crossed 10.

Same $E[Z(5)]$. Very different $T_{10}$. The one explosive spike realization ($Z = 28$) does all the work of preserving the mean, while the majority of spike trajectories lag behind the smooth trajectories. This is the "skewness tax."

**What you would see in simulation:** if you ran 10,000 simulations of each model and plotted:

- The **mean curve** (average $Z(t)$ across simulations at each $t$): **identical** for both models, exactly.
- The **median curve**: **shifted right** for the spike case.
- Individual curves: more spread out for the spike case.

The user's intuition that "the spike curves are shifted right" is correct for the typical/median trajectory, and this is what drives the difference in $E[T_n]$. The mean trajectory is unshifted, but $T_n$ is determined by where individual trajectories cross the threshold $n$, not by where the mean trajectory crosses.

### The CMJ limit theorem and the random variable $W$

CMJ theory provides the following almost-sure convergence result: conditional on non-extinction,

$$\frac{Z(t)}{e^{\alpha t}} \xrightarrow{a.s.} W \qquad \text{as } t \to \infty,$$

where $W > 0$ is a random variable that depends on the stochastic early phase of the epidemic. The random variable $W$ captures the "effective head start" of the epidemic: a large $W$ means the epidemic got off to a fast start; a small $W$ means a slow start. Once $Z(t)$ is large, the epidemic grows deterministically at rate $\alpha$, but the trajectory is permanently shifted by the early stochastic phase.

The key properties of $W$:

- **$E[W \mid \text{survival}]$** is determined by the intensity measure $\mu$ and is the same for both models (it depends on $A(\tau)$ and $\alpha$, which are shared).
- **$\text{Var}(W)$** depends on the full offspring process, not just the intensity measure. It differs between the smooth and spike models.

### Why the spike case has higher $\text{Var}(W)$

The variance of $W$ is determined by the second-moment structure of the offspring process. Define the "variance contribution" from each individual:

- **Smooth model:** Individual $i$'s offspring arrive at $k_i$ distinct times, spread across the generation interval. The contribution to $Z(t)$ from $i$'s offspring is a sum of $k_i$ terms at different times. By the independence of offspring times, the variance of this contribution is roughly $k_i$ times the variance of a single-offspring contribution (reduced by the spreading across time).

- **Spike model:** Individual $i$'s offspring all arrive at the same time $\tau_i^*$. The contribution to $Z(t)$ is $k_i$ individuals all launched from the same time point. The variance is higher because offspring are perfectly correlated in timing — they all either arrive "early" (contributing to fast growth) or "late" (contributing to slow growth) together.

Formally, the within-individual covariance of offspring contributions is zero for the smooth model (independent arrival times) but positive for the spike model (identical arrival times). This inflates $\text{Var}(W)$ for the spike case.

### Jensen's inequality: from $E[Z(t)]$ to $E[T_n]$

Define $T_n = \inf\{t : Z(t) \geq n\}$, the time to reach $n$ infections. Using the CMJ limit:

$$Z(T_n) \approx n \implies W \cdot e^{\alpha T_n} \approx n \implies T_n \approx \frac{\log n - \log W}{\alpha}$$

Taking expectations (conditional on survival):

$$E[T_n] \approx \frac{\log n - E[\log W]}{\alpha}$$

Now apply Jensen's inequality. Since $\log$ is concave:

$$E[\log W] \leq \log E[W]$$

with equality only when $W$ is deterministic. The gap $\log E[W] - E[\log W]$ is exactly the **expected log-deviation**, which increases with $\text{Var}(W)$. More precisely, for a random variable $W > 0$:

$$\log E[W] - E[\log W] = \frac{\text{Var}(W)}{2\,(E[W])^2} + O(\text{higher cumulants})$$

(this is the second-order expansion via the delta method, exact for log-normal $W$).

Since $`\text{Var}(W)_{\text{spike}}`$ > $`\text{Var}(W)_{\text{smooth}}`$ while $E[W]$ is the same:

$$E[\log W]_{\text{spike}} < E[\log W]_{\text{smooth}}$$

$$\implies E[T_n]_{\text{spike}} > E[T_n]_{\text{smooth}}$$

**The spike case reaches milestones later, on average, despite having the same expected population size at every time point.** The epidemic doesn't grow slower — it effectively starts later (in a stochastic sense), then grows at the same deterministic rate.

### Quantitative connection to the order-statistics mechanism

The order-statistics mechanism described in Findings 3 and 12 is the *micro-level* manifestation of the variance difference in $W$.

In the smooth case, an infector with $k$ offspring has first-transmission time $\min(\tau_1, \ldots, \tau_k)$, which is stochastically earlier than a single draw from $g(\tau)$. This "fastest horse wins" dynamic reduces the *effective* generation interval for the first link in each chain, compressing the early stochastic phase and reducing $\text{Var}(W)$.

In the spike case, the first-transmission time is $\tau^*$ regardless of $k$, so the effective generation interval equals the population generation interval. There is no order-statistics acceleration, and $\text{Var}(W)$ remains high.

The connection: the per-infector first-transmission advantage compounds multiplicatively across the first few generations of the branching process, translating into a multiplicative factor in $W$. The smooth case's $W$ is a product of terms that are individually less variable (each has a narrower first-transmission distribution), yielding a tighter distribution overall.

### When the effect is large vs. small

The magnitude of the delay $`\Delta T_n`$ = $`E[T_n]_{\text{spike}}`$ - $`E[T_n]_{\text{smooth}}`$ depends on:

1. **$R_0$**: The order-statistics advantage scales with the expected number of offspring $E[k] = R_0$. At $R_0 \leq 1.2$, most infectors produce 0 or 1 offspring, so $\min(\tau_1) = \tau_1$ — no advantage. At $R_0 = 5$, infectors typically produce 4–6 offspring, and $\min(\tau_1, \ldots, \tau_5)$ is much earlier than $\tau^*$. This explains the $R_0$-dependence in Finding 12 (delay grows from ~0% at $R_0 = 1.2$ to ~18% at $R_0 = 5$).

2. **$n$ (the threshold)**: The delay $\Delta T_n$ is approximately *constant* in $n$ (it's a time offset, not a rate difference). As a fraction of $T_n$, it decreases as $\log n$ grows. This explains the finding that percentage delays attenuate from $t_{10}$ to $t_{100}$ (Finding 12).

3. **Initial number of cases**: If the epidemic is seeded with $m$ initial cases (rather than 1), the variance of $W$ decreases by a factor of $\sim m$ (law of large numbers over independent branching processes), and $\Delta T_n \to 0$. The effect is purely a *small-number-of-cases* phenomenon.

4. **$\psi$ (for the Gamma convolutional model)**: At $\psi \to 0$ (delta spike), $\text{Var}(W)$ is maximised. At $\psi \to \alpha_{\text{total}}$ (smooth = $A$), $\text{Var}(W)$ is minimised. The delay interpolates continuously between these extremes.

### Key clarifications

**This is not a finite-population effect.** The branching process operates in an effectively infinite population (no susceptible depletion). The delay arises from starting with a small number of infected individuals, not from the total population being finite. An epidemic seeded into a population of $10^9$ with a single index case would show the same delay as one seeded into a population of $10^3$.

**The expected growth rate is truly identical.** The Malthusian parameter $\alpha$ is the same. The delay is a *time offset*, not a rate reduction. Once $Z(t)$ is large enough for the law of large numbers to apply (roughly $`Z(t) \gtrsim 50 - 100`$), both models grow at rate $\alpha$ and the offset is frozen.

**$E[Z(t)]$ is identical — but the median of $Z(t)$ is not.** This is the core resolution of the apparent paradox. $E[Z(t)]$ is the same at every $t$, starting from a single case, with no approximations. This is an exact consequence of the Campbell theorem (the intensity measure determines $E[Z(t)]$, and the intensity measure is the same for both models). However, the **median** (and other quantiles) of $Z(t)$ differs: the spike case has a lower median at any given $t$ during the early epidemic, compensated by a heavier right tail. When we look at simulation curves, we see the median behavior, not the mean. The mean curve would overlay exactly; the median curve is shifted right for the spike case. The delay in $E[T_n]$ arises because $T_n$ is determined by where individual trajectories cross the threshold — a quantity sensitive to the median/distribution of $Z(t)$, not just its mean. Jensen's inequality on the concave function $\log$ formalises this: higher variance in $W$ (the stochastic head-start factor) produces a later expected milestone time.

**Conditioning on establishment does not eliminate the effect.** Among epidemics that establish (avoid stochastic extinction), the spike case still reaches milestones later. Both models have the same extinction probability (both have Poisson($R_0$) offspring counts — the timing of offspring does not affect the offspring *distribution*). Conditioning on $W > 0$ (survival) changes the distribution of $W$ but preserves the ordering $`\text{Var}(W)_{\text{spike}} > \text{Var}(W)_{\text{smooth}}`$, so Jensen's inequality still applies.

---

## 18. Connection to Morris, Maclean & Black (2024): computing the time-shift distribution

**Reference:** Morris, Maclean & Black (2024), "Computation of random time-shift distributions for stochastic population models," *Journal of Mathematical Biology*.

### The shared framework

Morris et al. formalise exactly the object that drives the growth-delay phenomenon described in Section 17: the random variable $W$ from the branching-process limit theorem, and its associated **time-shift** $\tau$.

Their starting point is the same CMJ/CT-MBP convergence:

$$Z(t) \approx W \cdot v(t) \qquad \text{as } t \to \infty,$$

where $v(t) = c \, e^{\alpha t}$ is the deterministic leading-order solution (with $\alpha$ the Malthusian parameter and $c$ a normalisation constant). Each stochastic trajectory is approximately a *time-shifted* copy of the deterministic trajectory:

$$Z(t) \approx v(t + \tau), \qquad \tau = \frac{\log W - \log E[W]}{\alpha}.$$

The random time-shift $\tau$ is what separates a fast-starting epidemic from a slow one. In our framework:

- **Smooth model:** $\tau$ is tightly distributed (low $\text{Var}(W)$), so trajectories cluster near the deterministic solution.
- **Spike model:** $\tau$ is broadly distributed (high $\text{Var}(W)$), producing a wide spread of trajectory timings.

The growth delay discussed in Section 17 is precisely the statement that $E[\tau]_{\text{spike}}$ < $E[\tau]_{\text{smooth}}$ (i.e., the spike model's typical trajectory is shifted further *behind* the deterministic solution). This follows from Jensen's inequality on $\log W$: since $E[W]$ is the same for both models but $\text{Var}(W)$ is larger for the spike model, $E[\log W]_{\text{spike}}$ < $E[\log W]_{\text{smooth}}$, and hence $E[\tau]_{\text{spike}}$ < $E[\tau]_{\text{smooth}}$.

### What they compute (and how)

Morris et al. develop two numerical methods to compute the *full distribution* of $W$ (and hence $\tau$):

1. **PE (Probability Estimation) method:** Numerically inverts the Laplace-Stieltjes transform (LST) of $W$. The LST satisfies a functional equation derived from the branching structure:

   $$E[e^{-sW}] = g\left(\frac{1}{R_0}\int_0^\infty E[e^{-s \, e^{-\alpha u} W}] \, \mu(du)\right)$$

   (where $g$ is the offspring probability generating function and $\mu$ is the intensity measure). They solve this iteratively and invert to obtain the density of $W$. This gives the exact distribution up to numerical precision.

2. **MM (Moment Matching) method:** Computes the first few moments of $W$ analytically (from the branching process structure) and fits a **generalised gamma distribution** to match them. This is faster and gives a closed-form approximation that is typically very accurate.

### What we do that they don't

Morris et al. treat the offspring process as *given* and ask "what is the distribution of $W$?" Their paper does not consider the possibility of varying the offspring process while keeping the intensity measure (and hence the Malthusian parameter $\alpha$ and the deterministic trajectory $v(t)$) invariant.

This is precisely our contribution. We construct a family of offspring processes — parametrised by $\psi$ — that all share the same $A(\tau)$ (and hence the same $\alpha$, $R_0$, and generation-interval distribution) but differ in the within-individual correlation structure of offspring timing. This family traces a path through the space of branching processes that Morris et al.'s methods could be applied to, generating a *curve* of $W$-distributions indexed by $\psi$.

In their language: we hold the intensity measure $\mu(d\tau) = A(\tau) \, d\tau$ fixed and vary the offspring point process (from independent arrival times at $\psi = \alpha_{\text{total}}$ to perfectly synchronised arrivals at $\psi \to 0$). The deterministic trajectory $v(t)$ is invariant across this family; only the stochastic fluctuation $W$ changes.

### Their framework applied to our problem

**Peak timing and planning.** Morris et al. motivate their work with the problem of predicting when an epidemic will reach its peak — a quantity that depends on $\tau$ and hence on $W$. Our results add a new source of uncertainty to this prediction: even if $R_0$, the generation-interval distribution, and $\alpha$ are perfectly known, uncertainty about $\psi$ (how punctuated individual infectiousness is) translates into uncertainty about $\text{Var}(W)$ and hence about the spread of possible peak times.

**A computational route to $\text{Var}(W)$ for our models.** Section 17 argues qualitatively that $`\text{Var}(W)_{\text{spike}} > \text{Var}(W)_{\text{smooth}}`$, based on the within-individual covariance structure. Morris et al.'s PE or MM methods could be used to compute $\text{Var}(W)$ *exactly* for each $\psi$ in our Gamma convolutional family. This would:

- Quantify the growth delay $\Delta E[T_n]$ analytically (rather than relying on simulation).
- Provide the full distribution of $\tau$, enabling probabilistic statements like "there is a 90% chance the spike epidemic reaches 100 cases between days $X$ and $Y$."
- Validate the approximation $\Delta E[T_n] \approx \text{Var}(W) / (2 \alpha \, (E[W])^2)$ from the delta-method expansion in Section 17.

**The moment-matching route.** Morris et al. show that $W$ is often well-approximated by a generalised gamma distribution. If this holds for our family, the entire effect of $\psi$ on the time-shift distribution could be summarised by three parameters (shape, scale, power) as functions of $\psi$, giving a parsimonious description of how punctuated infectiousness affects epidemic timing uncertainty.

### Technical subtleties

**CT-MBP vs. CMJ.** Morris et al.'s implementation focuses on continuous-time Markovian branching processes (CT-MBP), where each individual has an exponentially-distributed lifetime and gives birth at a constant rate while alive. Our models are age-dependent (CMJ): the infectiousness profile $a_i(\tau)$ is a non-trivial function of infection age $\tau$. The *theory* (convergence to $W$, functional equations for the LST) applies to both, but their specific numerical algorithms would need adaptation for the non-Markovian case. The Gamma convolutional model could potentially be embedded in a Markovian framework via phase-type distributions (the Gamma profile with integer shape parameter is a sum of exponentials), making it amenable to their methods.

**The functional equation for $W$.** For a CMJ process with offspring point process $\xi$ (a random counting measure on $[0, \infty)$), the LST of $W$ satisfies:

$$E[e^{-sW}] = E\left[\prod_{j=1}^{N} E[e^{-s \, e^{-\alpha \tau_j} W}]\right]$$

where $N$ is the number of offspring and $\tau_1, \ldots, \tau_N$ are their birth times. For the smooth model ($\tau_j$ independent given $N$), the product factorises into a simpler form. For the spike model ($\tau_j = \tau^*$ for all $j$), the product collapses to $E[e^{-sW}]^N$ evaluated at a single random time. The difference in these functional equations is what drives the difference in $W$-distributions — and hence the growth delay.

**Independence structure matters.** A key insight connecting our work to theirs: the functional equation for $W$ depends on the *joint* distribution of offspring times, not just their marginal distribution. Two models can have the same number-of-offspring distribution and the same marginal birth-time distribution (i.e., the same $A(\tau)$) but different $W$ distributions if the offspring times are correlated differently. This is exactly the mechanism our $\psi$ parameter controls.

### Summary of the relationship

| Aspect | Morris et al. (2024) | Our framework |
|--------|----------------------|---------------|
| Core object | $W$ and time-shift $\tau$ | Same |
| What varies | The epidemic/population model | The offspring correlation structure ($\psi$), holding $A(\tau)$ fixed |
| What's computed | Full distribution of $W$ | Currently: simulated $E[T_n]$; could use their methods for $W$ |
| Key result | Computational methods (PE, MM) for $W$ | Qualitative ordering: $\text{Var}(W)$ increases as $\psi \to 0$ |
| Practical focus | Peak timing uncertainty | Growth delay, surveillance implications |
| Process type | CT-MBP (Markovian) | CMJ (age-dependent) |

Their computational machinery is complementary to our conceptual framework. We identify *why* $W$ varies across models with the same mean-field dynamics (offspring synchronisation), while they provide the tools to compute *how much* it varies. A natural collaboration or extension would apply their PE/MM methods to our Gamma convolutional family, producing exact $\text{Var}(W)$-vs-$`\psi`$ curves and validated time-shift distributions.

## Analytical computation of Var(W) and generalized gamma moment matching

Section 17 argued qualitatively that the spike model has higher $\text{Var}(W)$ than the smooth model, and Section 18 noted that Morris et al. (2024) provide computational methods for the $W$ distribution. Here we derive **exact closed-form expressions** for $\text{Var}(W)$ as a function of $\psi$ in our Gamma convolutional family and implement **generalized gamma moment matching** to approximate the full distribution of $W$ (and hence the time-shift $\tau$).

### Closed-form Malthusian parameter

For a Gamma($\alpha\_{\text{total}}$, $r$) generation-interval distribution, the Euler-Lotka equation $R_0 \cdot (r/(r + \alpha))^{\alpha\_{\text{total}}} = 1$ yields:

$$\alpha = r \cdot (R_0^{1/\alpha\_{\text{total}}} - 1)$$

With standard parameters ($R_0 = 2$, $\alpha\_{\text{total}} = 10$, $r = 2$): $\alpha \approx 0.1435$.

### Var(W) via the distributional fixed-point equation

$W$ satisfies the distributional fixed-point equation $W \stackrel{d}{=} \sum_j e^{-\alpha \tau_j} W_j$, where the $W_j$ are i.i.d. copies of $W$ independent of the offspring times $\tau_j$. Define $S = \sum_j e^{-\alpha \tau_j}$ (discounted offspring sum). The key variance formula is:

$$\text{CV}^2(W) = \frac{\text{Var}(W)}{(E[W])^2} = \frac{\text{Var}(S)}{1 - m_2}$$

where $m_2 = E[\sum_j e^{-2\alpha \tau_j}] = R_0 \cdot \rho_2^{\alpha\_{\text{total}}}$ (independent of $\psi$), with $\rho = r/(r+\alpha)$ and $\rho_2 = r/(r+2\alpha)$.

### The crucial $\psi$-dependent quantity: Var(S)

Using the decomposition $\tau_j = s + \epsilon_j$ (shared shift + independent jitter), where $s \sim \text{Gamma}(\alpha\_{\text{total}} - \psi, \beta)$ and $\epsilon_j \sim \text{Gamma}(\psi, \beta)$:

$$S = e^{-\alpha s} \cdot T, \qquad T = \sum_j e^{-\alpha \epsilon_j}$$

Since $s \perp T$:

$$\text{Var}(S) = \underbrace{m_2}_{\text{Poisson + jitter}} + \underbrace{R_0^2 \cdot \rho^{2\psi} \cdot [\rho_2^{\alpha\_{\text{total}} - \psi} - \rho^{2(\alpha\_{\text{total}} - \psi)}]}_{\text{synchronisation penalty}}$$

The first term ($m_2$) is constant across $\psi$ — it captures variance from the random offspring count and independent jitter. The second term is the **synchronisation penalty**: it vanishes when $\psi = \alpha\_{\text{total}}$ (smooth, $s = 0$) and is maximised when $\psi \to 0$ (spike, $s$ carries all variance). This confirms the qualitative argument from Section 17 with an exact formula.

### Time-delay formula

The expected time delay of model $\psi$ relative to the smooth limit is:

$$\Delta E[T_n] \approx \frac{\text{CV}^2(W)\_\psi - \text{CV}^2(W)\_{\text{smooth}}}{2\alpha}$$

This is independent of $n$ (the case threshold), consistent with the branching-process limit theorem: the time shift is established during the initial stochastic phase and persists unchanged as the epidemic grows.

### Variance of the time shift

The time shift relative to the smooth limit is $\tau = (\log W - \log E_c)/\alpha$, so $\text{Var}(\tau) = \text{Var}(\log W)/\alpha^2$. Since $\text{Var}(\log W)$ has no exact closed form from moments alone, we use two approaches.

**Delta method (first-order approximation).** Taylor-expanding $\log W$ around $E[W]$ gives $\text{Var}(\log W) \approx \text{CV}^2(W)$, so:

$$\text{Var}(\tau) \approx \frac{\text{CV}^2_{\text{cond}}}{\alpha^2}$$

This separates into a $\psi$-independent baseline plus a $\psi$-dependent penalty:

$$\text{Var}(\tau) \approx A + B\bigl(\sigma^{(1-\psi)\alpha\_{\text{total}}} - 1\bigr)$$

where the constants are:

- $A = V_0 = \dfrac{m_2 - q}{\alpha^2(1 - m_2)}$ — the variance in the smooth limit ($\psi = 1$), arising from random offspring count and independent jitter alone.
- $B = \dfrac{1 - q}{\alpha^2(1 - m_2)}$ — the coefficient of the synchronisation penalty.
- $\sigma = \dfrac{\eta^2}{2\eta - 1}$ where $\eta = R_0^{1/\alpha\_{\text{total}}}$ — the per-stage discount variability ratio $\rho_2/\rho^2$.
- $m_2 = R_0 \cdot \rho_2^{\alpha\_{\text{total}}}$ — the second-moment discount (defined in the Var(S) derivation above).
- $q$ = extinction probability of the Poisson($R_0$) branching process.
- $\alpha = r(\eta - 1)$ = Malthusian growth rate.

Since $\alpha = \alpha\_{\text{total}}(\eta - 1)/T$, we have $\text{Var}(\tau) \propto T^2$: the variance of the time shift scales with the square of the mean generation interval.

**Exact under the Gamma approximation.** If $W \mid \text{surv} \sim \text{Gamma}(k, \lambda)$ (the 2-moment approximation from below), then $\text{Var}(\log W) = \psi_1(k)$, where $\psi_1$ is the trigamma function. This gives:

$$\text{Var}(\tau) = \frac{\psi_1(k)}{\alpha^2}, \qquad k = \frac{1}{\text{CV}^2\_{\text{cond}}}$$

For large $k$, $\psi_1(k) \approx 1/k = \text{CV}^2_{\text{cond}}$, recovering the delta method. For small $k$ (small $\psi$, highly skewed $W$), the trigamma gives a more accurate answer — the relative error of the delta method is roughly $1/(2k)$.

### Third moment and generalized gamma matching

$E[W^3]$ requires three auxiliary quantities (all closed-form for the Gamma convolutional model):

- $m_3 = R_0 \cdot \rho_3^{\alpha\_{\text{total}}}$ where $\rho_3 = r/(r + 3\alpha)$
- $m\_{21} = R_0^2 \cdot \rho_2^\psi \cdot \rho^\psi \cdot \rho_3^{\alpha\_{\text{total}} - \psi}$
- $m\_{111} = R_0^3 \cdot \rho^{3\psi} \cdot \rho_3^{\alpha\_{\text{total}} - \psi}$

Then: $E[W^3] = (3 \cdot m\_{21} \cdot E[W^2] \cdot E[W] + m\_{111} \cdot (E[W])^3) / (1 - m_3)$.

Following Morris et al.'s moment-matching (MM) approach, we fit the first three moments of $W | W > 0$ to a **generalized gamma distribution** with parameters $(a, d, p)$:

$$f(x; a, d, p) = \frac{p}{a} \left(\frac{x}{a}\right)^{d-1} \exp\left(-\left(\frac{x}{a}\right)^p\right) \Big/ \Gamma(d/p)$$

The generalized gamma moments $E[X^k] = a^k \cdot \Gamma((d+k)/p) / \Gamma(d/p)$ provide a system of equations solvable by numerical optimisation. This gives a parsimonious three-parameter description of the $W$-distribution for each $\psi$.

### Results

**CV²(W) is monotonically decreasing in $\psi$.** This confirms the qualitative ordering: spike ($\psi \to 0$) produces the most variable $W$, smooth ($\psi \to \alpha\_{\text{total}}$) the least. The synchronisation penalty term provides the quantitative explanation.

**Expected time delays of 1–3 days** accumulate for punctuated profiles relative to smooth. For $\psi = 0.5$ (near-spike), the expected delay is approximately 2–3 days; this narrows monotonically as $\psi \to \alpha\_{\text{total}}$.

**Generalized gamma fits** match the first three moments of $W | W > 0$ well for all $\psi$ values, producing overlaid densities that track the empirical simulation histograms closely. The fitted $(a, d, p)$ trace a smooth curve as $\psi$ varies.

**Simulation validation.** Branching-process simulations (10,000 replicates per $\psi$, $N = 50{,}000$ to avoid depletion) confirm the analytical $\text{CV}^2(W)$ within $\pm 5\%$ relative error and the time-delay formula within $\pm 10\%$.

### Figures

- **`fig_W_variance.pdf`**: Three panels. (A) CV²(W) vs $\psi$, analytical curve with simulation points. (B) Var(S) decomposition into constant Poisson+jitter and decreasing synchronisation penalty. (C) Expected time delay $\Delta E[T_n]$ vs $\psi$.
- **`fig_W_gengamma.pdf`**: Four panels ($\psi = 0.5, 2, 6, 9.5$) showing empirical histograms of $W | W > 0$ with overlaid generalized gamma densities from moment matching.
- **`fig_W_timeshift.pdf`**: Time-shift distributions $\tau = (\log W - \log E[W])/\alpha$ for four $\psi$ values, showing how the spread narrows as $\psi \to \alpha\_{\text{total}}$.

### Connection to Morris et al. (2024)

This section implements the "moment-matching route" discussed in Section 18. The closed-form expressions for Var(S) and the moment recursions exploit the specific structure of our Gamma convolutional family — in particular, the decomposition $\tau_j = s + \epsilon_j$ and the resulting conditional independence of $T$ and $s$. Morris et al.'s PE method (Laplace-Stieltjes transform inversion) would give the exact distribution to arbitrary precision, but the moment-matching approach produces an excellent approximation with far less computational cost and provides interpretable parameters.

**Script**: `code/analytical_W_moments.R`

---

## 20. Mean-field ODE vs stochastic expectation

The ODE solution is often loosely described as "the expected epidemic trajectory," but this is not exactly correct. The mean-field ODE and the true expectation $E[Z(t)]$ of the stochastic process are related but distinct quantities, and the distinction matters in finite populations with nonlinear dynamics.

### The general principle

- **Mean-field (ODE):** Replace all random variables with their means *inside* the dynamical equations, then solve the resulting deterministic system. This answers: "what trajectory would the system follow if every quantity were always exactly at its expected value?"
- **True expectation:** Solve the full stochastic process, then average over realisations. This gives $E[X(t)]$ — the quantity that converges as we average more and more simulations at each time point.

These two operations commute for *linear* systems: if $dx/dt = Ax$, then $E[x(t)]$ satisfies the same ODE. But they diverge whenever the dynamics contain *nonlinearities*, because $E[f(X)] \neq f(E[X])$ for nonlinear $f$ when $X$ has variance (Jensen's inequality).

### Application to SEIR

The only nonlinearity in the SEIR system is the transmission term $\beta S I / N$. The true expected force of infection is:

$$E\left[\frac{\beta S I}{N}\right] = \frac{\beta}{N}\bigl(E[S] \cdot E[I] + \text{Cov}(S, I)\bigr)$$

The ODE sets $\text{Cov}(S, I) = 0$. But mechanistically, every new infection simultaneously decrements $S$ and increments $I$ (or $E$), so $\text{Cov}(S, I) < 0$ always. This means the true expected infection rate is *lower* than the ODE assumes:

$$E[S \cdot I] < E[S] \cdot E[I]$$

Consequently, the true $E[Z(t)]$ grows slower than the ODE solution at every time point, and this bias accumulates over the course of the epidemic. The mean-field ODE systematically *overestimates* the expected cumulative incidence trajectory.

### When does the bias matter?

- **Large $N$:** The covariance $\text{Cov}(S, I)$ is $O(N)$ while $E[S] \cdot E[I]$ is $O(N^2)$, so the relative bias is $O(1/N)$ and the ODE becomes exact as $N \to \infty$.
- **Small $N$ or fast epidemics:** With $N = 1000$ and $R_0 \geq 5$, a substantial fraction of the population is infected rapidly, the $S$-$I$ correlation is non-negligible relative to $N$, and the ODE visibly overestimates the stochastic mean.
- **Burstier profiles (spike) show larger bias:** The spike profile concentrates all secondary infections at a single moment, creating larger bursts that deplete susceptibles more abruptly. This amplifies the magnitude of $\text{Cov}(S, I)$ relative to smoother profiles where infections are spread over time.

### Empirical confirmation

In `code/1_1_basicmetrics.R`, we overlay the ODE solution (blue) and the empirical mean of established stochastic trajectories (black) on `fig_cuminf_overlay`. The black curve is consistently pulled below the blue curve, with the gap largest for the spike profile. This gap is *not* primarily a conditioning artefact (conditioning on establishment has a comparatively small effect when the establishment probability is high). It reflects the finite-population mean-field bias described above, and would shrink with increasing $N$.

---

## 21. Distribution of the individual reproduction number under time-varying contacts

When contacts are constant, every individual has the same reproduction number $R_0$. With time-varying contacts, the individual reproduction number $R_i$ becomes a random variable whose distribution depends on (i) the contact pattern, (ii) the infectiousness profile (parametrised by $\psi$), and (iii) the infection time. Here we derive the distribution of $R_i$ for sinusoidal contacts.

### Setup

Suppose contacts vary sinusoidally: $c(t) = \bar{c}(1 + \varepsilon \sin(\omega t))$, where $\varepsilon \in [0, 1)$ is the contact amplitude and $\omega = 2\pi / P$ is the angular frequency for a period $P$. A person infected at time $t_0$ has individual reproduction number:

$$R_i(t_0) = R_0 \int_0^\infty g_i(\tau) \cdot \frac{c(t_0 + \tau)}{\bar{c}} \, d\tau$$

where $g_i(\tau)$ describes how person $i$'s infectiousness is distributed over time after infection.

### Smoothing via the characteristic function

Expanding the sinusoidal contact term and integrating against $g_i$:

$$R_i(t_0) = R_0 \left(1 + \varepsilon \cdot \rho \cdot \sin(\omega t_0 + \varphi)\right)$$

where $\rho = |\varphi_g(\omega)|$ is the modulus of the characteristic function of the generation interval distribution $g$, and $\varphi$ is the corresponding phase shift. The quantity $\rho \in [0, 1]$ measures how much the infectiousness profile *smooths out* the contact variation: a more spread-out profile averages over more of the contact cycle, reducing $\rho$.

### Dependence on $\psi$ in the Gamma convolutional family

In our model, each infection attempt has timing $\tau_j = s + \varepsilon_j$, where $s \sim \text{Gamma}(\alpha\_{\text{total}} - \psi, \beta)$ is a shared component and $\varepsilon_j \sim \text{Gamma}(\psi, \beta)$ is independent jitter. Only the jitter provides smoothing — the shared component shifts all attempts together without averaging. Therefore:

$$\rho(\psi) = |\varphi\_{\text{Gamma}(\psi, \beta)}(\omega)| = \left(\frac{r}{\sqrt{r^2 + \omega^2}}\right)^\psi$$

- **Spike ($\psi \to 0$):** $\rho = 1$. No jitter, no smoothing. The full contact amplitude $\varepsilon$ feeds through to $R_i$.
- **Smooth ($\psi = \alpha\_{\text{total}}$):** $\rho = (r / \sqrt{r^2 + \omega^2})^{\alpha\_{\text{total}}} < 1$. Maximum smoothing; the contact variation is strongly attenuated.
- **Intermediate $\psi$:** Smooth interpolation between these extremes.

### Marginal distribution over uniform infection times

If infection times are uniformly distributed over the contact cycle, then $\sin(\omega t_0 + \varphi)$ has the **arcsine distribution** on $[-1, 1]$, with density $f(x) = 1/(\pi \sqrt{1 - x^2})$. Therefore:

$$R_i \sim R_0 (1 + \varepsilon \cdot \rho(\psi) \cdot X), \qquad X \sim \text{Arcsine}[-1, 1]$$

This gives:

- $E[R_i] = R_0$ (same for all $\psi$)
- $\text{Var}(R_i) = R_0^2 \varepsilon^2 \rho(\psi)^2 / 2$
- Support: $[R_0(1 - \varepsilon \rho), \; R_0(1 + \varepsilon \rho)]$

The overdispersion in $R_i$ is **monotonically decreasing in $\psi$**: spike profiles produce the widest distribution of individual reproduction numbers, smooth profiles the narrowest.

### Role of contact frequency

The contact frequency $\omega$ interacts with $\psi$ through $\rho(\psi) = (r / \sqrt{r^2 + \omega^2})^\psi$:

- **Slow variation ($\omega \to 0$):** $\rho \to 1$ for all $\psi$. Contact changes are slow relative to the generation interval, so all profiles "see" the same local contact rate. No differentiation between profiles.
- **Fast variation ($\omega \to \infty$):** $\rho \to 0$ for all $\psi > 0$. Contact changes are fast relative to the generation interval, so all profiles average over many cycles. No overdispersion from contacts.
- **Intermediate $\omega$ (contact period comparable to generation interval):** Maximum differentiation between profiles. This is the regime where the punctuatedness $\psi$ has the strongest effect on $R_i$ variability.

### Extension to non-sinusoidal contacts

For a general periodic contact pattern $c(t) = \bar{c}(1 + \sum_k \varepsilon_k \sin(k\omega t + \psi_k))$, each Fourier harmonic is smoothed independently by $\rho(\psi)$ evaluated at frequency $k\omega$. The marginal distribution of $R_i$ is then the distribution of a sum of weighted arcsine random variables (which is no longer arcsine in general). For complex contact patterns, simulation is more practical than analytical computation.
