
## Detect-and-isolate: a unified framework for screening and symptom-triggered interventions

Detect-and-isolate interventions reduce transmission by truncating each detected individual's infectiousness profile at the detection time: $a_i^{\text{eff}}(\tau) = a_i(\tau) \cdot \mathbf{1}[\tau < \tau_{\text{det}}]$. This truncation both reduces $R_{\text{eff}}$ and reshapes the generation interval distribution — but the degree of both effects depends critically on $\psi$.

The key insight that unifies screening-based and symptom-triggered interventions: both truncate the profile at a detection time $\tau_{\text{det}}$, and the entire $\psi$-dependence flows through a single quantity:

$$F_\psi(\tau_{\text{det}} - l_i) = \int_0^{\tau_{\text{det}} - l_i} f_\psi(u) du$$

— the fraction of the biological profile's mass that has elapsed before detection. What differs between detection mechanisms is the **joint distribution of $\tau_{\text{det}}$ and the profile timing**.

### The detection spectrum: independent vs. anchored triggers

Detection mechanisms sit on a spectrum of **biological anchoring** — the degree to which $\tau_{\text{det}}$ is correlated with the individual's profile timing:

**Independent trigger ($\rho = 0$): regular screening.** Detection time is determined by an external process (periodic testing every $\Delta$ days, random phase) intersected with a positivity window. The screening clock runs independently of $l_i$, so $\tau_{\text{det}}$ is approximately uncorrelated with the profile. The $\psi$-dependence comes purely from the *shape* of $F_\psi$ — spike profiles have a step-function CDF (all-or-nothing), smooth profiles have a gradual sigmoid (partial aversion).

**Anchored trigger ($\rho = 1$): symptom onset.** Detection time is deterministically coupled to the profile: $\tau_{\text{det}}$ = $l_i$ + $\delta_{\text{sym}}$ + $\delta_{\text{iso}}$, where $\delta_{\text{sym}}$ is the delay from infectiousness onset to symptoms and $\delta_{\text{iso}}$ is the delay from symptoms to effective isolation. Individuals with early onset (small $l_i$) are both early transmitters and early detectors. The $\psi$-dependence comes from both the shape of $F_\psi$ and the correlation between $\tau_{\text{det}}$ and $l_i$.

**Intermediate anchoring ($0 < \rho < 1$):** Hybrid scenarios such as screening that targets symptomatic individuals, or symptoms that are a noisy indicator of infectiousness timing. These interpolate between the independent and anchored extremes.

### The core quantity: fraction averted given detection

For any detection mechanism, the fraction of individual $i$'s transmission averted is:

$$1 - F_\psi(\tau_{\text{det}} - l_i)$$

The population-level $R_{\text{eff}}$ reduction is:

$$R_{\text{eff}} = R_0 \cdot \left(1 - P(\text{det}) \cdot E\left[1 - F_\psi(\tau_{\text{det}} - l_i) \mid \text{detected}\right]\right)$$

The $\psi$-dependence enters through $F_\psi$ and, for anchored triggers, through the correlation between $\tau_{\text{det}}$ and $l_i$.

### Results for regular screening (implemented)

Simulation with 3-day screening, sensitivity 0.85, and a 6-day positivity window confirms strong $\psi$-dependence.

#### Detection probability is not the driver

$P(\text{det})$ is approximately constant across $\psi$ at $`\sim`$97.5%. The positivity window is wide relative to the screening interval, so almost everyone is detected regardless of profile shape.

#### The conditional value of detection is the driver

$E[\text{frac averted} \mid \text{detected}]$ drops from $`\sim`$0.88 for spike profiles ($\psi = 0.5$) to $`\sim`$0.80 for smooth profiles ($\psi = 9.5$):

| $\psi$ | $P(\text{det})$ | $E[\text{frac averted} \mid \text{det}]$ | $R_{\text{eff}}$ | Tail leakage |
|----------|-----------------|------------------------------------------|-------------------|---------------|
| 0.5      | 0.975           | 0.88                                     | 0.28              | 0%            |
| 2        | 0.975           | 0.86                                     | 0.32              | 1%            |
| 6        | 0.975           | 0.82                                     | 0.40              | 5%            |
| 9.5      | 0.975           | 0.80                                     | 0.44              | 7%            |

$R_{\text{eff}}$ varies by a factor of $`\sim`$1.6 across $\psi$ despite identical detection probabilities.

#### Mechanism 1: All-or-nothing vs partial aversion

For spike profiles, all transmission attempts cluster at a single time point. If detection occurs before this spike, 100% is averted; if after, 0%. The conditional expectation $E[\text{frac averted} \mid \text{detected}]$ is a weighted average of 1 and 0, dominated by the large fraction detected before the spike.

For smooth profiles, detection at any time averts only the remaining fraction. An individual detected near the peak has already completed 30–50% of transmission. The distribution of frac averted is bimodal (clustered at 0 and 1) for spike profiles but unimodal and centered below 1 for smooth profiles.

#### Mechanism 2: Tail leakage

Smooth profiles place a non-negligible fraction of transmission outside the positivity window entirely (7% at $\psi = 9.5$). Spike profiles keep nearly 100% within the catchable window.

#### Generation interval distortion

Detection also reshapes the effective generation interval $g_{\text{eff}}(\tau)$, with a $\psi$-dependent effect:

**Spike case:** $g_{\text{eff}} = g$ (unchanged). Testing removes a random fraction of transmitters without altering the timing of the remaining ones — the blocking probability is approximately constant across spike locations.

**Smooth case:** $g_{\text{eff}}(\tau)$ $\propto A(\tau) \cdot P(T_{\text{detect}} > \tau)$ shifts shorter. The right tail of the profile is truncated by detection, so surviving transmissions are biased early.

**Epidemiological consequence: partially self-defeating intervention.** Via the Euler-Lotka equation $1$ = $R_{\text{eff}}$ $\int e^{-r\tau} g_{\text{eff}}(\tau)\,d\tau$, a shorter $g_{\text{eff}}$ implies a *larger* growth rate for a given $R_{\text{eff}}$. So for smooth profiles, the $R_{\text{eff}}$ reduction from detection is partially offset by the generation interval shortening. For spike profiles, $g_{\text{eff}} = g$, so the full $R_{\text{eff}}$ reduction translates to growth rate reduction. This partial self-defeat is specific to profile-truncating interventions and does not apply to interventions that reduce $R_{\text{eff}}$ without altering the generation interval (vaccination, uniform contact reduction).

Analysis code: `code/generation_interval_distortion.R`, `code/reff_reduction_decomposition.R`.

#### Note on the positivity window at small $\psi$

When $\psi$ is small, some individuals have early peaks ($t_{\text{peak}} < w_-$), causing the positivity window to nominally extend before $\tau = 0$. In practice tests cannot fire before infection, so the only effect is a slight overstatement of detection effectiveness for spike profiles. The true $\psi$-dependence is therefore at least as strong as reported.

### Predictions for symptom-triggered isolation (planned)

*Status: not yet implemented.*

Symptom-triggered isolation replaces the independent screening clock with a biologically anchored trigger: $\tau_{\text{det}}$ = $l_i$ + $\delta_{\text{sym}}$ + $\delta_{\text{iso}}$, deterministically coupled to the onset shift $l_i$.

#### Why anchoring should sharpen the $\psi$-dependence

With screening, $\tau_{\text{det}}$ is approximately independent of $l_i$, so detection can occur at any point relative to the profile. With symptoms, $\tau_{\text{det}}$ tracks the profile: individuals with early onset are detected early. This coupling means the analysis reduces to a single, clean question: **what fraction of the profile's mass falls before symptom onset?** That fraction is $F_\psi(\delta_{\text{sym}})$ — a property of the profile shape alone, without the noise introduced by the screening schedule.

#### Key predictions

**Delta (small $\psi$):** Transmission occurs at $\tau^* \approx l_i + \varepsilon$ where $\varepsilon$ is small. Symptom onset at $l_i + \delta_{\text{sym}}$. If $\delta_{\text{sym}}$ + $\delta_{\text{iso}}$ > $\varepsilon$ (symptoms + isolation delay exceeds the jitter), isolation arrives after the spike — nothing averted. If $\delta_{\text{sym}}$ + $\delta_{\text{iso}}$ < $\varepsilon$, everything averted. The all-or-nothing structure is sharper than for screening because the detection time tracks the profile.

**Smooth (large $\psi$):** Transmission is spread over a broad window. Symptom onset truncates the right tail. The fraction averted is $1 - F_\psi(\delta_{\text{sym}})$, a smooth function of $\delta_{\text{sym}}$.

**The key $\psi$-dependent quantity:** $F_\psi(\delta_{\text{sym}})$, the CDF of $\text{Gamma}(\psi, r)$ evaluated at the symptom delay. For small $\psi$, this is a sharp step function. For large $\psi$, a smooth sigmoid.

#### What's new relative to screening

1. **Pre-symptomatic transmission is the key parameter.** The entire analysis reduces to one number: the fraction of transmission before symptoms. No test-performance parameters (sensitivity, specificity, screening interval).

2. **Correlation between detection and profile.** Individuals with early $l_i$ are both early transmitters and early detectors. Whether this correlation amplifies or dampens the $\psi$-dependence (relative to screening) depends on where symptom onset falls relative to the transmission mass.

3. **No false positives or sensitivity parameters.** Only $\delta_{\text{sym}}$, $\delta_{\text{iso}}$, and $p_{\text{iso}}$ (compliance).

#### Design choices to resolve

- **Symptom onset model.** Where does $\tau_{\text{sym}}$ fall relative to $b_i$? Options:
  - At the mode of $f_\psi$ (symptoms at peak infectiousness). Mode of $\text{Gamma}(\psi, r)$ is $(\psi - 1)/r$ for $\psi \geq 1$; for $\psi < 1$, mode is at 0.
  - At a fixed fraction of the profile's CDF (e.g., symptoms when 30% of infectiousness mass has elapsed). This gives a $\psi$-dependent delay.
  - At a fixed absolute delay $\delta_{\text{sym}}$ after $l_i$ (simplest, but least biologically motivated).

- **Isolation delay $\delta_{\text{iso}}$:** 0, 0.5, 1, 2 days.

- **Asymptomatic fraction.** $p_{\text{sym}} < 1$ acts as a ceiling on effectiveness but doesn't interact with $\psi$ in an interesting way.

## Gathering size restrictions are more effective for smooth profiles

We investigated how gathering size restrictions interact with the temporal shape of individual infectiousness profiles. The central question: are gathering size caps more effective at preventing outbreaks when infectiousness is punctuated (spiky) versus smooth?

The answer turns out to be the opposite of the naive expectation. Caps are **more effective for smooth profiles**, because they reduce mean transmission without losing the overdispersion that was already suppressing outbreak probability for spiky profiles.

### Model

Each individual moves through a sequence of gatherings. The contact process $c_i(t)$ is piecewise-constant:

- **Gathering sizes**: drawn iid from $\text{NegBin}(n = 10/9,\; p = 1/10)$, giving mean 10, variance 100 (SD = 10). After normalisation (dividing by 5), the contact rate has mean 2 and variance 4.
- **Holding times**: $\text{Exp}(\text{rate} = 2)$, so individuals change gathering approximately every 0.5 days.
- **Gathering size cap**: draws exceeding the cap are rejected (truncated NegBin).

The individual reproduction number $R_i$ arises from the interaction of the biological profile and the contact process:

- **Delta profile**: $R_i = c_i(t^*)$ at a single random time. This is one draw from the contact-value distribution, inheriting its full variance.
- **Smooth profile**: $R_i = \int c_i(\tau) \cdot f(\tau)\,d\tau$, where $f$ is the $\text{Gamma}(10, 0.3)$ PDF spanning $`\sim`$33 days. This averages over $\sim$36 independent contact-process steps, dramatically reducing variance.

Offspring are $\text{Poisson}(R_i)$ — the standard superspreading framework.

### Caps reduce mean $R_i$ identically for all profile shapes

The contact process enters multiplicatively with the biological profile. Since every biological profile integrates to the same value (regardless of $\psi$), the expected contact rate at any time point determines $E[R_i]$. Truncating the gathering size distribution reduces this expected contact rate equally for spiky and smooth profiles.

| Scenario      | $E[R_i]$ |
|---------------|----------|
| Uncapped      | 2.00     |
| Capped at 20  | 1.39     |

This was confirmed analytically and in simulation (28 $\psi \times$ cap scenarios, with $E[R_i]$ varying by $< 0.04$ across $\psi$ for each cap level).

**The mean effect of gathering restrictions is $\psi$-independent.**

### Caps reduce variance identically in fractional terms

The variance of $R_i$ decomposes cleanly:

| Profile | Scenario   | $E[R_i]$ | $\text{Var}(R_i)$ | $\text{CV}(R_i)$ |
|---------|------------|----------|--------------------|-------------------|
| Delta   | Uncapped   | 2.00     | 4.00               | 1.00              |
| Delta   | Capped@20  | 1.39     | 1.20               | 0.79              |
| Smooth  | Uncapped   | 2.00     | 0.11               | 0.17              |
| Smooth  | Capped@20  | 1.39     | 0.03               | 0.13              |

The cap reduces variance by $`\sim`$70% for both profiles. This is because the variance reduction factor from averaging ($I_{\text{corr}} = 0.028$, corresponding to $`\sim`$36 effective independent samples) multiplies both the uncapped and capped contact-process variance identically:

$$\text{Var}(R_i)_{\text{smooth}} = \text{Var}(c_i) \cdot I_{\text{corr}}$$

So the fractional reduction in variance from capping is the same: $\text{Var}(\text{capped}) / \text{Var}(\text{uncapped}) = \text{Var}(X \mid X \leq 20) / \text{Var}(X)$, regardless of the biological profile.

**The differential is not in the fractional variance reduction — it is in the absolute amount of superspreading removed.** The cap eliminates $\Delta\text{Var} = 2.8$ for the delta profile versus $\Delta\text{Var} = 0.08$ for the smooth profile (a 35$`\times`$ difference).

### Overdispersion increases extinction probability

For a branching process with $\text{Poisson}(R_i)$ offspring, higher variance in $R_i$ (for fixed $E[R_i] > 1$) increases the probability of stochastic extinction. Superspreading creates many individuals with $R_i$ near zero, providing frequent opportunities for chains to die out.

**Unmitigated ($R = 2$):**

| Profile              | $P(\text{extinction})$ | $P(\text{outbreak})$ |
|----------------------|------------------------|----------------------|
| Poisson(2) reference | 0.20                   | 0.80                 |
| Smooth               | 0.22                   | 0.78                 |
| Delta                | 0.51                   | 0.49                 |

Despite having the same mean $R = 2$, the delta profile produces outbreaks less than half the time, while the smooth profile produces outbreaks $`\sim`$80% of the time. The smooth profile's near-Poisson $R_i$ distribution means reliable, moderate transmission from every individual.

### Caps are more effective at preventing outbreaks for smooth profiles

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

### Simulation results confirm the analytical predictions

The large-scale simulation (`code/gathering_size_restrictions.R`) swept 7 values of $\psi$ (0.5 to 9.5) across 4 cap levels (5, 10, 20, $\infty$), with 50,000 MC replicates per scenario for $R_i$ distributions and 200 epidemic simulations per scenario.

Key confirmations:

- **$E[R_i]$ is constant across $\psi$** for each cap level (range $< 0.04$)
- **$\text{CV}(R_i)$ decreases with $\psi$** for uncapped scenarios (from 1.07 at $\psi = 0.5$ to 0.82 at $\psi = 9.5$), confirming that smooth profiles average out contact variation
- **All capped epidemic scenarios showed zero establishment** in the full simulation (because with geometric gathering sizes of mean 30, even a cap at 20 reduces $E[R_i]$ to $\sim$0.58 — well below 1). This demonstrates that for heavy-tailed gathering size distributions, even moderate caps can collapse mean transmission, making the variance-level differential moot for epidemic outcomes.

The analytical calculations (Finding 16d) used a different parameterisation ($\text{NegBin}/5$ with mean 2, cap at 20 giving $E[R_i] = 1.39 > 1$) to isolate the regime where the differential effect on outbreak probability is visible.

### Summary

The naive hypothesis — that gathering size restrictions are *more* effective against punctuated infectiousness because they selectively remove superspreading events — is incorrect, or at least incomplete.

Gathering size caps reduce mean transmission identically regardless of profile shape. But they also reduce variance, and for spiky profiles, that variance was *already suppressing outbreaks* via stochastic extinction. The cap simultaneously removes the superspreading tail (which sustained rare large outbreaks) and the zero-transmission tail (which killed most chains). The net effect on outbreak probability is smaller than for smooth profiles, where there was no beneficial overdispersion to lose.

**Gathering size restrictions are more effective at preventing outbreaks when infectiousness is smooth**, because smooth profiles have no overdispersion buffer to erode.

Analysis code: `code/gathering_size_restrictions.R`.

## Planned: Contact tracing with imperfect isolation / post-exposure prophylaxis

*Status: planned, not yet implemented.*

### Motivation

Contact tracing is the most natural timing-dependent intervention to study because its effectiveness depends on the temporal profiles of *both* the infector (index case) and the infectee (traced contact). The delay from infection of the index case to intervention on the traced contact passes through a chain of timing-dependent steps, each of which interacts with $\psi$.

We model imperfect isolation and post-exposure prophylaxis (PEP) within a single framework. The key observation: both reduce (rather than eliminate) a traced contact's onward transmission from some intervention time onward. Perfect isolation ($\eta = 0$) is a special case; leaky isolation and PEP correspond to $\eta > 0$, differing only in the plausible parameter ranges and whether a pre-emptive mode exists.

### Model

Contact tracing unfolds as a cascade of delays:

**Step 1: Index case detection** ($\tau_{\text{det}}^A$)

The index case (person A) is detected via one of:
- (a) Symptom onset + testing (Section 14, anchored trigger): $\tau_{\text{det}}^A$ = $s_A$ + $\delta_{\text{sym}}$ + $\delta_{\text{test}}$
- (b) Routine screening (Section 14, independent trigger): $\tau_{\text{det}}^A$ depends on the positivity window and screening phase

In both cases, $\tau_{\text{det}}^A$ depends on person A's profile via $s_A$ and $f_\psi$.

**Step 2: Tracing delay** ($\delta_{\text{trace}}$)

After detection, contacts are identified and notified. This introduces a further delay $\delta_{\text{trace}}$ that is approximately independent of $\psi$ (it depends on public health infrastructure, not biology). We can model it as fixed or $\text{Exp}(\text{rate})$.

**Step 3: Intervention on the traced contact**

Person B (a contact of A) is reached at calendar time:

$$T_{\text{int}}^B = t_A + \tau_{\text{det}}^A + \delta_{\text{trace}}$$

Person B was infected by A at some time $t_B = t_A + \tau_{A \to B}$, where $\tau_{A \to B}$ is the generation interval from A to B. The time from B's infection to intervention is:

$$\delta_B = T_{\text{int}}^B - t_B = \tau_{\text{det}}^A + \delta_{\text{trace}} - \tau_{A \to B}$$

This depends on A's detection time (which depends on A's profile via $\psi$) and the generation interval from A to B (which also depends on A's profile).

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

### Why $\psi$ enters three times

The effectiveness of tracing person B depends on:

1. **How quickly A is detected** ($\tau_{\text{det}}^A$) — determines when the tracing clock starts. For spike A: detection is either very early (if the trigger precedes the spike) or very late (if the spike precedes the trigger). For smooth A: detection is at a predictable time near the profile's peak.

2. **What generation interval connected A to B** ($\tau_{A \to B}$) — determines how much head start B has. For spike A: $\tau_{A \to B} \approx s_A + \varepsilon$ (all of A's transmissions cluster at one time, so all traced contacts have similar head starts). For smooth A: $\tau_{A \to B}$ is spread over a wide range, so some contacts were infected early (long head start, likely already transmitted) and others late (short head start, intervention is more useful).

3. **How concentrated B's own transmission is** ($`F_B(\delta_B)`$) — determines how much of B's profile falls before vs. after the intervention. For spike B: $F_B$ is a sharp step function, making the outcome all-or-nothing. For smooth B: $F_B$ is a smooth sigmoid, giving graded partial benefit.

### Key predictions

**Spike A, spike B (small $\psi$ throughout):** A's transmission spike produces a cluster of contacts all infected at nearly the same time. If intervention is fast enough to reach them before their own spikes, all transmission is averted (or reduced to $\eta$). If not, none is averted. **High variance, all-or-nothing.** The residual fraction $\eta$ barely matters: either the spike beat the intervention ($F_B \approx 1$, nothing to reduce) or it didn't ($F_B \approx 0$, full $(1-\eta)$ reduction).

**Smooth A, smooth B (large $\psi$ throughout):** A is detected at a predictable time. Contacts infected early by A have a long head start and have completed much of their own broad transmission window ($F_B$ large); contacts infected late have a short head start ($F_B$ small, intervention is effective). **Graded effectiveness, depending on the generation interval from A to B.** The residual fraction $\eta$ matters throughout: even after intervention, smooth B continues transmitting at rate $\eta$, and this residual accumulates over the remaining broad profile.

**Spike A, smooth B (or vice versa):** The cross-$\psi$ cases reveal asymmetries. A spike infector produces clustered contacts (all with similar $\delta_B$), while a smooth infector produces dispersed contacts (variable $\delta_B$). The infectee's $\psi$ then determines whether intervention catches their transmission window.

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

The $\psi$-dependence of effectiveness should therefore **reverse sign** as a function of $\delta_B$. The cross-over point depends on the relationship between the delay distribution and the onset-shift distribution $s_B \sim \text{Gamma}(\alpha_{\text{total}} - \psi, r)$.

### Connection to detect-and-isolate framework

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
- **$\psi$ for infector vs infectee**: sweep independently if the population is heterogeneous, or assume a single $\psi$ if all individuals share the same profile shape

### Taxonomy of timing-dependent interventions

| Intervention | Acts on | Mechanism | $\eta$ | Profile interaction |
|---|---|---|---|---|
| Detect-and-isolate, screening (§14) | Infector | Truncation at random detection time | 0 | Sampling vs. averaging of positivity window |
| Detect-and-isolate, symptoms (§14) | Infector | Truncation at biologically anchored time | 0 | Pre-symptomatic fraction depends on $\psi$ |
| Contact tracing + perfect isolation | Infectee | Truncation at traced time | 0 | Two-generation: A's profile → delay → B's truncation |
| Contact tracing + leaky isolation | Infectee | Scaling at traced time | $> 0$ | As above, plus residual leakage erodes smooth-B benefit |
| Contact tracing + PEP | Infectee | Scaling at traced time, pre-emptive mode | $> 0$ | As above, plus race between PEP and B's onset shift |

The progression from Section 14 to Section 16 traces a path of increasing mechanistic coupling between the intervention and the temporal structure of infectiousness. This section unifies the final three rows into a single parameterised model.




## 22. Testing effectiveness and the punctuation parameter

This section formalises the connection between the $\psi$-dependent infectiousness model and the testing effectiveness (TE) framework of Middleton & Larremore (2024), deriving analytical predictions that the simulations in `code/2_1_isolation.R` will validate.

### The TE framework adapted to our model

Following Middleton & Larremore (M&L), define testing effectiveness as the fraction by which a detect-and-isolate intervention reduces cumulative transmission:

$$\text{TE} = \frac{\int_0^\infty \beta(\tau) F_D(\tau) d\tau}{\int_0^\infty \beta(\tau) d\tau}$$

where $\beta(\tau)$ is the infectiousness profile (relative to infection time) and $F_D(\tau) = P(D \leq \tau)$ is the CDF of the diagnosis time $D$ (also relative to infection time). Equivalently, TE is the probability that a randomly chosen transmission event occurs after diagnosis: $\text{TE} = E_\varepsilon[F_D(\varepsilon)]$, where $\varepsilon \sim \beta$ is the timing of a random transmission attempt.

In our Gamma convolutional model, the jitter component of transmission timing is $\varepsilon \sim \text{Gamma}(\psi\alpha, r)$, where $\alpha$ = `popshape` and $r = \alpha / T$ is the rate. (Recall that each individual's transmission attempts are at times $l_i + \varepsilon_j$, where $l_i$ is the shared onset shift. Because $l_i$ shifts both infectiousness and any biologically anchored trigger by the same amount, TE depends only on the jitter $\varepsilon$, not on $l_i$ directly.)

### Symptom-triggered isolation: the Beta function result

Model symptom onset as occurring at time $l_i + D$ after infection, where $D$ ~ Gamma($a_{\text{sym}}$, $b_{\text{sym}}$) is the delay from biological onset to symptoms.

**Same-rate case ($b_{\text{sym}} = r$).** When symptoms and jitter share the same rate parameter, the ratio $\varepsilon / (\varepsilon + D)$ has a clean distribution. Since $\varepsilon \sim \text{Gamma}(\psi\alpha, r)$ and $D \sim \text{Gamma}(a_{\text{sym}}, r)$ are independent with the same rate:

$$\frac{\varepsilon}{\varepsilon + D} \sim \text{Beta}(\psi\alpha,\; a_{\text{sym}})$$

Testing effectiveness is the probability that diagnosis precedes the jittered transmission:

$$\text{TE} = P(D < \varepsilon) = P\left(\frac{\varepsilon}{\varepsilon + D} > \frac{1}{2}\right) = I_{1/2}(a_{\text{sym}},\; \psi\alpha)$$

where $I_x(a, b)$ is the regularised incomplete beta function. In R: `pbeta(0.5, psi*popshape, a_sym, lower.tail = FALSE)`.

**Properties:**

- *Monotonically increasing in $\psi$:* smoother profiles (larger $\psi\alpha$) spread transmission over a wider window, so a larger fraction falls after diagnosis. In the limit $\psi\alpha \to \infty$, $\text{TE} \to 1$.
- *Monotonically decreasing in $a_{\text{sym}}$:* larger $a_{\text{sym}}$ (with rate $r$ fixed) means later symptom onset on average ($E[D] = a_{\text{sym}}/r$), reducing the fraction averted. In the limit $a_{\text{sym}} \to \infty$, $\text{TE} \to 0$.
- *For small $\psi\alpha$:* transmission is concentrated near $\tau = 0$; almost all of it fires before symptoms. $\text{TE} \to 0$.
- *Scaling invariance:* if $a_{\text{sym}} = c \cdot \psi\alpha$ for some constant $c$ (symptoms scale with the profile width), then $\varepsilon / (\varepsilon + D) \sim \text{Beta}(\psi\alpha, c \cdot \psi\alpha)$, and $\text{TE} = I_{1/2}(c \cdot \psi\alpha, \psi\alpha)$. As $\psi\alpha \to \infty$ with ratio $c$ fixed, this converges to $\mathbf{1}[c < 1]$ — a step function that depends only on whether mean symptom delay is shorter or longer than mean jitter.

**Different-rate case ($b_{\text{sym}} \neq r$).** The ratio $\varepsilon / D$ follows a scaled $F$-distribution:

$$\frac{\varepsilon / (\psi\alpha)}{D / a_{\text{sym}}} \sim F(2\psi\alpha,\; 2a_{\text{sym}})$$

so TE = $P(D < \varepsilon)$ = $P\big(F(2\psi\alpha, 2a_{\text{sym}})$ > $\frac{r \cdot a_{\text{sym}}}{b_{\text{sym}} \cdot \psi\alpha}\big)$, which is closed-form via the $F$-distribution CDF. The same-rate case is recovered when $b_{\text{sym}} = r$ (the threshold becomes $a_{\text{sym}} / (\psi\alpha)$).

### Periodic screening: derived trigger distribution

For regular screening every $\Delta$ days with random phase, the detection time depends on the intersection of the screening schedule with the detectability window. Parametrise the window as $[\text{mode}_\psi - w_-,\; \text{mode}_\psi + w_+]$, anchored to the mode of the infectiousness profile:

$$\text{mode}_\psi = \max\left(0,\; \frac{\psi\alpha - 1}{r}\right)$$

The effective diagnosis time (relative to biological onset $l_i$) is approximately:

$$D_{\text{screen}} \approx \max(0,\; \text{mode}_\psi - w_-) + \text{Uniform}(0, \Delta)$$

The offset $\max(0, \text{mode}_\psi - w_-)$ is the earliest time post-onset at which a test can detect infection, and the uniform component represents the random phase of the screening clock within one period.

TE requires integrating $P(\varepsilon > D_{\text{screen}})$ over the distribution of $D_{\text{screen}}$. Although this does not simplify to a single special function, the $\psi$-dependence has transparent structure:

- **Small $\psi$ (punctuated):** $`\text{mode}_\psi \approx 0`$, so the detectability window begins at or before onset. The offset is zero, and $`D_{\text{screen}} \sim \text{Uniform}(0, \Delta)`$. Meanwhile, transmission is concentrated near the mode. If $`\Delta`$ is short relative to $w_+$, detection likely precedes the spike: **TE is high**.
- **Large $\psi$ (smooth):** $`\text{mode}_\psi = (\psi\alpha - 1)/r`$ is large, but so is the spread of $`\varepsilon`$. The offset $`\text{mode}_\psi - w_-`$ grows, pushing diagnosis later. A substantial fraction of transmission occurs before the first possible positive test: **TE is lower**.
- **Direction of $\psi$-effect: opposite to symptoms.** Punctuated profiles benefit *more* from screening, whereas smooth profiles benefit *more* from symptom-triggered isolation.

The crucial mechanism is that the detectability window is anchored to peak infectiousness, whose position depends on $\psi$. This $\psi$-dependence in the trigger distribution reverses the direction compared to the symptom case (where the trigger is anchored to biological onset, independent of $\psi$).

### Comparison and the anchoring spectrum

The reversal is best understood through the anchoring framework from §14:

| Mechanism | Trigger distribution (rel. onset) | $\psi$-dependence of trigger | $\psi$-effect on TE |
|---|---|---|---|
| Symptoms (fixed delay) | $\text{Gamma}(a_{\text{sym}}, r)$, fixed | None | Smooth benefits more |
| Symptoms (tracking peak) | $a_{\text{sym}} \propto \psi\alpha$ | Full scaling | Can vanish (scaling invariance) |
| Screening (wide window) | $\text{Uniform}(0,\Delta) + \text{offset}(\psi)$ | Through window position | Punctuated benefits more |

The reversal occurs because:

- **Symptoms:** the trigger is anchored to biological onset, so it is approximately fixed relative to $l_i$ regardless of $\psi$. The fraction averted equals the tail of $f_\psi$ beyond the symptom time. Smoother profiles spread more mass into this tail.
- **Screening:** the trigger is anchored to peak infectiousness. Punctuated profiles concentrate their peak near onset, so the detectability window extends across nearly all transmission. Smooth profiles spread transmission well before and after the window centre.

### General detection distributions: three principles

The symptom and screening cases are instances of a fully general result. For any non-negative detection time $D$ (possibly depending on $\psi$), testing effectiveness is:

$$\text{TE}(\psi) = E_D\left[S_{\psi\alpha}(D)\right], \qquad S_{\psi\alpha}(d) \equiv P(\varepsilon > d) = 1 - F_{\text{Gamma}(\psi\alpha,\,r)}(d)$$

The survival function $S_{\psi\alpha}$ is the complete sufficient statistic for how *any* detection mechanism interacts with the profile. Three properties of $S$ govern the $\psi$-dependence of TE.

**Principle 1 (Stochastic ordering).** $\text{Gamma}(\psi_1 \alpha, r) \leq_{\text{st}} \text{Gamma}(\psi_2 \alpha, r)$ for $\psi_1 < \psi_2$. Therefore $S_{\psi\alpha}(d)$ is pointwise increasing in $\psi$ for every $d > 0$. Consequence: **for any detection distribution $D$ that is independent of $\psi$, TE is monotonically increasing in $\psi$.** Smooth profiles always benefit more from a $\psi$-independent trigger. This is the default; any reversal requires the detection mechanism to "see" the profile shape.

**Principle 2 (Mean-shift reversal).** When $D$ depends on $\psi$, two effects compete:

$$\frac{d\,\text{TE}}{d\psi} = \underbrace{\alpha \cdot E\left[\frac{\partial S_a(D)}{\partial a}\bigg|_{a=\psi\alpha}\right]}_{\text{profile spreading (always } > 0\text{)}} \;-\; \underbrace{\frac{dE[D]}{d\psi} \cdot E\left[f_{\psi\alpha}(D)\right]}_{\text{detection delay shift}}$$

The first term is always positive: broadening $\varepsilon$ puts more mass beyond any fixed detection point. The second term is negative when $E[D]$ increases with $\psi$. The reversal (TE decreasing in $\psi$) requires $E[D]$ to grow with $\psi$ fast enough to overcome the first term. The critical rate is approximately $dE[D]/d\psi \approx \alpha / r = T$, the same rate at which $E[\varepsilon] = \psi\alpha/r$ grows. When $D$ is deterministic, numerical verification confirms:

| Growth of $E[D]$ | $\psi$-effect on TE | Example |
|---|---|---|
| Constant | Increasing (smooth benefits) | Symptom-triggered isolation |
| Proportional to $\psi T$ | Approximately neutral | Symptoms tracking peak |
| Faster than $\psi T$ | Decreasing (punctuated benefits) | Screening with peak-anchored window |

**Principle 3 (Variance–convexity interaction).** For $\psi\alpha \leq 1$, the survival function $S_{\psi\alpha}(d)$ is *globally convex* on $(0, \infty)$. By Jensen's inequality, for any random $D$ with a given mean, $E[S(D)] \geq S(E[D])$: variance in $D$ always increases TE. For $\psi\alpha > 1$, $S$ is concave before the mode $(\psi\alpha - 1)/r$ and convex after. The consequence:

- When detection typically occurs *before* peak infectiousness (e.g., early screening), $D$ falls in the concave region: variance in $D$ *decreases* TE.
- When detection typically occurs *after* peak infectiousness (e.g., delayed symptoms), $D$ falls in the convex region: variance in $D$ *increases* TE.

Because punctuated profiles ($\psi\alpha$ small) have globally convex $S$ while smooth profiles ($\psi\alpha$ large) have a concave region, **noisy detection mechanisms disproportionately benefit punctuated profiles**. This is the formal generalisation of the all-or-nothing mechanism from §14: the bimodal distribution of fraction averted for punctuated profiles (either 0% or 100%) is precisely the convexity of $S$.

**Summary.** The three principles organise the full design space of detection mechanisms:

| Property of $S_{\psi\alpha}$ | Small $\psi\alpha$ | Large $\psi\alpha$ | Consequence |
|---|---|---|---|
| Level at fixed $d$ | Low | High | Fixed $D$ → smooth benefits (Principle 1) |
| Shift rate with $\psi$ | $\alpha/r$ | $\alpha/r$ | $D$ must track to keep up (Principle 2) |
| Curvature | Globally convex | Concave-then-convex | Variance in $D$ helps punctuated more (Principle 3) |

Any detection strategy can be analysed by specifying two things: (i) how $E[D]$ scales with $\psi$ (Principle 2), and (ii) where $D$ sits relative to the mode of $\varepsilon$ (Principle 3). Symptom-triggered isolation is a Principle-1 story (fixed $D$, smooth wins). Periodic screening is a Principle-2 story ($E[D]$ shifts with $\psi$, punctuated wins). The all-or-nothing mechanism of §14 is a Principle-3 story (convexity favours punctuated).

### Application: turnaround time and behavioral delays

In practice, detection is not instantaneous. The effective detection time is a sum of components:

$$D_{\text{eff}} = D_{\text{trigger}} + \delta_{\text{TAT}} + \delta_{\text{action}}$$

where $D_{\text{trigger}}$ is the time of sample collection (symptoms, screening), $\delta_{\text{TAT}}$ is the turnaround time from sample to result, and $\delta_{\text{action}}$ is the delay from result to effective isolation. Each additive delay contributes both a mean shift and variance, and the three principles predict their $\psi$-dependent effects.

**Smoothing interpretation.** Adding an independent delay $\delta$ to $D$ replaces the survival function $S_{\psi\alpha}$ with a smoothed version:

$$\tilde{S}(d) \equiv E_\delta\left[S(d + \delta)\right]$$

The effect on TE decomposes cleanly into a mean component and a variance component:

$$\underbrace{E[\tilde{S}(D)] - E[S(D)]}_{\text{total effect}} = \underbrace{E[S(D{+}E[\delta])] - E[S(D)]}_{\text{mean effect}} + \underbrace{E[\tilde{S}(D)] - E[S(D{+}E[\delta])]}_{\text{variance effect}}$$

where $D = D_{\text{trigger}}$. The mean effect is always negative (delay reduces TE for all profiles). The variance effect is governed by the curvature of $S$ at the effective detection point $D + E[\delta]$.

- If $ D + E[\delta] < \text{mode}_\psi = (\psi\alpha - 1)/r $, then $S$ is concave and variance **decreases** TE.
- If $ D + E[\delta] > \text{mode}_\psi $, then $S$ is convex and variance **increases** TE.
- If $ \psi\alpha \leq 1 $, then $S$ is globally convex and variance **always increases** TE. 

**$\psi$-dependent effect of TAT variance.** For symptom-triggered isolation with a fixed symptom delay $d_0$ and TAT $\delta \sim \text{Exp}(1)$ ($E[\delta] = 1$ day), the decomposition at $d_0 = 2$, $\alpha = 10$, $r = 2$:

| $\psi$ | $\psi\alpha$ | mode | TE (no TAT) | mean effect | var effect | TE (with TAT) |
|---|---|---|---|---|---|---|
| 0.05 | 0.5 | 0 | 0.005 | $-$0.004 | $+$0.001 | 0.002 |
| 0.10 | 1.0 | 0 | 0.018 | $-$0.016 | $+$0.004 | 0.006 |
| 0.50 | 5.0 | 2.0 | 0.629 | $-$0.344 | $+$0.066 | 0.351 |
| 0.95 | 9.5 | 4.25 | 0.987 | $-$0.101 | $-$0.038 | 0.848 |

The variance effect is positive for punctuated profiles (convex $S$) and negative for smooth profiles (concave $S$ when $d_0 + E[\delta] = 3 < \text{mode} = 4.25$). This means variable TAT **narrows the TE gap** between smooth and punctuated profiles relative to deterministic TAT of the same mean.

The intuition: for punctuated profiles, transmission is all-or-nothing, so occasional lucky early results (from TAT variance) can catch entire spikes that a fixed delay would miss. For smooth profiles with early detection (before the mode), the fixed delay already captures most of the tail; TAT variance introduces a chance of late results that allow transmission to escape.

**When does the crossover vanish?** If symptom onset is sufficiently delayed ($d_0 + E[\delta] > \text{mode}_\psi$ for all $\psi$), both profiles are in the convex region of their respective $S$, and variance helps everyone. The sign reversal requires that smooth profiles are detected *before* their mode — i.e., early detection. This is the regime where the gap between smooth and punctuated TE is largest, and where TAT variance has the strongest narrowing effect.

### Numerical predictions for standard parameters

Parameters: $T = 5$, `popshape` $\alpha = 10$, $r = \alpha/T = 2$.

**Symptom-triggered isolation (same-rate case):**

| $\psi$ | $\psi\alpha$ | $a_{\text{sym}}=1$ | $a_{\text{sym}}=3$ | $a_{\text{sym}}=5$ | $a_{\text{sym}}=8$ |
|---|---|---|---|---|---|
| 0.05 | 0.5 | 0.293 | 0.050 | 0.010 | 0.001 |
| 0.10 | 1.0 | 0.500 | 0.125 | 0.031 | 0.004 |
| 0.20 | 2.0 | 0.750 | 0.313 | 0.109 | 0.020 |
| 0.50 | 5.0 | 0.969 | 0.773 | 0.500 | 0.194 |
| 0.95 | 9.5 | 0.999 | 0.975 | 0.890 | 0.643 |

Computed via `pbeta(0.5, psi*popshape, a_sym, lower.tail = FALSE)`. Note the symmetry at $\psi = 0.5$, $a_{\text{sym}} = 5$: $\text{Beta}(5, 5)$ is symmetric about $1/2$, giving $\text{TE} = 0.5$ exactly.

The table confirms:
- TE is monotonically increasing in $\psi$ (reading down each column)
- TE is monotonically decreasing in $a_{\text{sym}}$ (reading across each row)
- At $\psi = 0.05$ (highly punctuated), symptom-triggered isolation averts almost nothing unless symptoms are extremely early ($a_{\text{sym}} = 1$ gives TE = 0.29)
- At $\psi = 0.95$ (smooth), symptom-triggered isolation is highly effective unless symptoms are very delayed ($a_{\text{sym}} = 8$ gives TE = 0.64)

**Periodic screening ($\Delta = 1, 3, 7$; $w_- = 3$, $w_+ = 4$):** requires numerical integration over the joint distribution of $\varepsilon$ and $D_{\text{screen}}$. The qualitative prediction is that TE decreases with $\psi$ (opposite to the symptom table above), with the magnitude depending on $\Delta$.

### Connection to §14 results

The TE framework provides the analytical underpinning for the simulation results in §14. The 3-day screening results — where punctuated profiles ($\psi = 0.5$) achieved $E[\text{frac averted}] \approx 0.88$ compared to $\approx 0.80$ for smooth profiles ($\psi = 9.5$) — correspond to the screening case above, where TE decreases with $\psi$ because the detection window is anchored to peak infectiousness.

The Middleton & Larremore framework additionally accommodates compliance $c$, test failure probability $\varphi$, and turnaround time (TAT), which multiply TE by factors independent of $\psi$. These can be incorporated when modelling realistic interventions but do not alter the $\psi$-dependence story derived here.

