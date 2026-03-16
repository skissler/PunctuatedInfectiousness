## Isolation triggered by detection: a unified framework

We study how three detection-and-isolation mechanisms — regular screening, symptom-triggered isolation, and contact tracing — interact with the punctuation parameter $\psi$. For each mechanism, we derive the fraction of transmission averted and identify where $\psi$ enters.

### Setup and notation

Recall from the theory notes: individual $i$ is infected at calendar time $t_i$ and has

- onset shift $l_i \sim \text{Gamma}((1-\psi)\alpha, \beta)$,
- biological profile $b_i(\tau) = f_\psi(\tau - l_i)$, where $f_\psi = \text{Gamma}(\psi\alpha, \beta)$ is the individual profile density,
- infectiousness $a_i(\tau) = R_0 \cdot b_i(\tau)$ (constant contacts),
- infection attempts $\chi_i \sim \text{Poisson}(R_0)$ at times $\tau_j = l_i + \varepsilon_j$, where $\varepsilon_j \overset{\mathrm{iid}}{\sim} \text{Gamma}(\psi\alpha, \beta)$.

We define the **mode of the individual biological profile** relative to biological onset:

$$m_\psi = \max\left(0, \frac{\psi\alpha - 1}{\beta}\right)$$

The mode of $a_i(\tau)$ (the time since infection at which infectiousness peaks) is $\tau_i^* = l_i + m_\psi$. This is individual-specific through $l_i$ but has the same offset $m_\psi$ from onset for everyone.

We also define the CDF and survival function of the individual biological profile:

$$F_\psi(x) = P(\text{Gamma}(\psi\alpha, \beta) \leq x), \qquad S_\psi(x) = 1 - F_\psi(x)$$

with $F_\psi(x) = 0$ and $S_\psi(x) = 1$ for $x \leq 0$.

### The general isolation model

Suppose individual $i$ is detected at time $\tau_{\text{det}}$ after their infection. After detection, isolation begins with a delay $\delta_{\text{act}}$ (capturing test turnaround time and imperfect temporal adherence):

$$\tau_{\text{iso}} = \tau_{\text{det}} + \delta_{\text{act}}$$

We adopt the convention $\tau_{\text{det}} = \infty$ (and hence $\tau_{\text{iso}} = \infty$) when detection never occurs.

Perfect temporal adherence corresponds to $\delta_{\text{act}} = 0$. Imperfect adherence can be modelled as $\delta_{\text{act}} \sim$ some non-negative distribution (e.g., $\text{Exp}(\lambda_{\text{act}})$ or a point mass).

Isolation reduces subsequent infectiousness by a factor $\eta \in [0, 1]$:

```math
a_i^{\text{eff}}(\tau) = \begin{cases} 
a_i(\tau) & \tau < \tau_{\text{iso}} \\ 
(1 - \eta) \cdot a_i(\tau) & \tau \geq \tau_{\text{iso}} 
\end{cases}
```

where $\eta = 1$ is perfect isolation and $\eta = 0$ is no effect.

The **fraction of individual $i$'s transmission averted** is

$$\theta_i = \eta \cdot S_\psi(\tau_{\text{iso}} - l_i)$$

Since $S_\psi(\infty) = 0$, undetected individuals (for whom $\tau_{\text{iso}} = \infty$) automatically contribute $\theta_i = 0$. For detected individuals, $\theta_i$ is the product of the isolation quality $\eta$ and the fraction of the biological profile remaining after isolation begins.

**Testing effectiveness** (population-level fraction of transmission averted):

$$\text{TE} = E[\theta_i] = P(\text{det}) \cdot E[\theta_i \mid \text{det}] + P(\text{not det}) \cdot 0 = P(\text{det}) \cdot E[\theta_i \mid \text{det}]$$

**Effective reproduction number:**

$$R_{\text{eff}} = R_0 (1 - \text{TE})$$

### Connection to Middleton & Larremore (2024)

Middleton & Larremore (M&L) define testing effectiveness as

$$\text{TE} = \frac{\int E[\beta(\tau) \cdot F_D(\tau)] d\tau}{\int E[\beta(\tau)] d\tau}$$

where $\beta(\tau)$ is the infectiousness profile and $F_D(\tau)$ is the CDF of diagnosis time. In our notation, $\beta(\tau) = R_0 \cdot b_i(\tau)$ and $F_D(\tau) = P(\tau_{\text{det}} \leq \tau)$. With perfect isolation ($\eta = 1$) and no action delay ($\delta_{\text{act}} = 0$), our TE reduces to M&L's expression. The key adaptations from their framework:

1. **Fixed-interval screening** rather than Poisson-rate screening. M&L model tests at rate $1/\Delta$; we model tests at exact intervals $\Delta$. This matters when $\Delta$ is comparable to the detectability window width.
2. **Detectability window** rather than explicit viral kinetics. M&L link detectability to a viral load model; we parametrise the window directly as $[d_{\text{pre}}, d_{\text{post}}]$ relative to the mode of the infectiousness profile.
3. **Explicit $\psi$-dependence.** M&L take the infectiousness profile as given; we vary its punctuation while keeping $A(\tau)$ fixed, which is the source of all our results.

---

## Trigger 1: Regular screening

### Detectability window

A person is testable (detectable by the assay) during a window anchored to the mode of their individual infectiousness profile:

$$\mathcal{W}_i = [\tau_i^* - d_{\text{pre}},\; \tau_i^* + d_{\text{post}}] = [l_i + m_\psi - d_{\text{pre}},\; l_i + m_\psi + d_{\text{post}}]$$

where $d_{\text{pre}} \geq 0$ and $d_{\text{post}} \geq 0$ are the half-widths of the detectability window before and after the mode, and we set $\max(0, \cdot)$ on the left endpoint to ensure the window doesn't extend before infection. The total window width is $w = d_{\text{pre}} + d_{\text{post}}$.

### Screening schedule

Individual $i$ takes tests at times $\phi_i, \phi_i + \Delta, \phi_i + 2\Delta, \ldots$ where:

- $\Delta$ is the screening interval (e.g., 3 days),
- $\phi_i \sim \text{Uniform}(0, \Delta)$ is the random phase of the screening clock relative to infection time.

This models a person who tests regularly (every $\Delta$ days) with a phase that is effectively random relative to their infection time.

### Detection

Each test that falls within $`\mathcal{W}_i`$ independently detects the person with probability $p_{\text{sens}}$ (test sensitivity). The detection time is the first positive test:

$$\tau_{\text{det}} = \min\\{\phi_i + n\Delta : \phi_i + n\Delta \in \mathcal{W}_i, U_n < p_{\text{sens}}, n = 0, 1, 2, \ldots\\}$$

where $U_n \overset{\mathrm{iid}}{\sim} \text{Uniform}(0,1)$. If no test yields a positive, the person is undetected.

**Number of tests in the window.** The number of screening times that fall in $\mathcal{W}_i$ is $\lfloor w / \Delta \rfloor$ or $\lceil w / \Delta \rceil$, depending on alignment. For $w > \Delta$, at least one test is guaranteed to fall in the window.

**Detection probability.** With $k$ tests in the window:

$$P(\text{det}) = 1 - (1 - p_{\text{sens}})^k$$

For $p_{\text{sens}} = 1$, detection is certain whenever $w \geq \Delta$ (at least one test falls in the window, and it's always positive). When $w < \Delta$, there is a probability $w / \Delta$ that any test falls in the window at all.

### Fraction averted conditional on detection

Given detection at time $`\tau_{\text{det}} \in \mathcal{W}_i`$, isolation begins at $`\tau_{\text{iso}} = \tau_{\text{det}} + \delta_{\text{act}}`$. The fraction averted is

$$\theta_i = \eta \cdot S_\psi(\tau_{\text{iso}} - l_i)$$

The random quantity $\tau_{\text{det}} - l_i$ (time from biological onset to detection) is approximately

$$\tau_{\text{det}} - l_i \approx m_\psi - d_{\text{pre}} + U$$

where $U$ is the offset of the first positive test within the window. For $p_{\text{sens}} = 1$ and $\Delta \leq w$, $U$ is approximately $\text{Uniform}(0, \Delta)$, so

$$\tau_{\text{det}} - l_i \approx m_\psi - d_{\text{pre}} + \text{Uniform}(0, \Delta)$$

### How $\psi$ enters

The $\psi$-dependence flows through $S_\psi$, the survival function of $\text{Gamma}(\psi\alpha, \beta)$:

- **Spike ($\psi \to 0$):** $S_\psi$ is approximately a step function at 0. If $\tau_{\text{iso}} - l_i < 0$ (isolation before onset), $\theta = \eta$; if $\tau_{\text{iso}} - l_i > 0$, $\theta = 0$. Outcome is **all-or-nothing**: either the test catches the person before their spike fires, averting everything, or it doesn't, averting nothing.

- **Smooth ($\psi \to 1$):** $S_\psi$ is a gradual sigmoid. Detection at any time averts a partial fraction of the remaining transmission. The outcome is **graded**: even detection near the peak averts the tail.

- **Direction of $\psi$-effect for screening.** The key observation is that the detectability window is anchored to the mode $\tau_i^* = l_i + m_\psi$. For spike profiles, $m_\psi \approx 0$ and all transmission is near $l_i$, so the window starts near the beginning of infectiousness — the person is detectable (and catchable) from the start. For smooth profiles, $m_\psi$ is large, and a substantial fraction of transmission occurs before the mode (and hence before or during the window). So **screening tends to benefit spike profiles more**.

---

## Trigger 2: Symptom-triggered isolation

### Symptom timing model

A symptomatic individual (probability $p_{\text{sym}}$) develops symptoms at time

$$\tau_{\text{sym}} = l_i + m_\psi + \delta_{\text{sym}}$$

where $\delta_{\text{sym}}$ is a random variable describing the timing of symptom onset relative to the mode of the individual's infectiousness profile. Here:

- $\delta_{\text{sym}} > 0$: symptoms arise after peak infectiousness,
- $\delta_{\text{sym}} < 0$: symptoms arise before peak infectiousness,
- $\delta_{\text{sym}} = 0$: symptoms coincide with peak infectiousness.

We impose $\tau_{\text{sym}} > 0$ (symptoms cannot precede infection), which is automatically satisfied for reasonable parameters since $l_i + m_\psi \geq 0$.

Isolation follows at $\tau_{\text{iso}} = \tau_{\text{sym}} + \delta_{\text{act}}$. The fraction of transmission averted is

$$\theta_i = \eta \cdot S_\psi(m_\psi + \delta_{\text{sym}} + \delta_{\text{act}})$$

Note that $l_i$ drops out: the onset shift $l_i$ shifts both the profile and the symptom trigger by the same amount, so the fraction averted depends only on the position of symptoms relative to the profile shape (through $m_\psi$, $\delta_{\text{sym}}$, and $\delta_{\text{act}}$), not on the individual's onset shift.

### Choice of distribution for $\delta_{\text{sym}}$

The framework is agnostic about the distribution of $\delta_{\text{sym}}$. Three natural choices, in order of complexity:

**Option A: Fixed delay (deterministic).** $\delta_{\text{sym}} = \mu_{\text{sym}}$ for all individuals. This is the simplest case and isolates the core $\psi$-dependence cleanly. The fraction averted reduces to $\theta_i = \eta \cdot S_\psi(m_\psi + \mu_{\text{sym}} + \delta_{\text{act}})$, which is non-random (given $\delta_{\text{act}}$). This is the right starting point.

**Option B: Normal jitter.** $\delta_{\text{sym}} \sim \mathcal{N}(\mu_{\text{sym}}, \sigma_{\text{sym}}^2)$, truncated to ensure $\tau_{\text{sym}} > 0$. The two parameters have clean interpretations:
- $\mu_{\text{sym}}$: mean symptom timing relative to peak infectiousness (positive = symptoms lag behind peak; negative = symptoms precede peak),
- $\sigma_{\text{sym}}$: between-individual variability in symptom timing.

The truncation is negligible when $l_i + m_\psi + \mu_{\text{sym}} \gg \sigma_{\text{sym}}$ (which holds for all but the most extreme parameters). This is the recommended default for stochastic analysis.

**Option C: Gamma offset.** $\delta_{\text{sym}} = -d_0 + G$ where $G \sim \text{Gamma}(a_{\text{sym}}, b_{\text{sym}})$ and $d_0 \geq 0$. This is always $\geq -d_0$, avoiding the need for truncation, and is right-skewed (long tail of late symptom onset). The three parameters $(d_0, a_{\text{sym}}, b_{\text{sym}})$ allow flexible asymmetric distributions but are less parsimonious than the Normal.

**Recommendation.** Use Option A (fixed $\mu_{\text{sym}}$) for the main analysis presented in the paper. This isolates the pure $\psi$-dependence and yields the cleanest results. Use Option B (Normal) for robustness checks and to show that between-individual variability in symptom timing doesn't qualitatively change the findings. The Normal is preferred over the Gamma offset because (i) it has fewer parameters, (ii) it allows symptoms before the mode without an extra shift parameter, and (iii) symmetry is a reasonable default when we have no strong prior on skewness.

### Pre-symptomatic transmission fraction

For fixed delay $\delta_{\text{sym}} = \mu_{\text{sym}}$ and perfect immediate isolation ($\eta = 1$, $\delta_{\text{act}} = 0$), the fraction of individual $i$'s transmission that occurs before symptoms is

$$F_\psi(m_\psi + \mu_{\text{sym}})$$

and the fraction averted (= fraction occurring after symptoms) is

$$\theta = S_\psi(m_\psi + \mu_{\text{sym}}) = 1 - F_\psi(m_\psi + \mu_{\text{sym}})$$

**Behaviour in the limits:**

- **Spike ($\psi \to 0$).** $m_\psi \to 0$ and $f_\psi \to \delta(0)$. So $S_\psi(\mu_{\text{sym}}) \to \mathbf{1}[\mu_{\text{sym}} < 0]$. If symptoms come before the mode ($\mu_{\text{sym}} < 0$), all transmission is averted. If after ($\mu_{\text{sym}} > 0$), nothing is averted. **All-or-nothing**, determined by the sign of $\mu_{\text{sym}}$.

- **Smooth ($\psi \to 1$).** $m_\psi = (\alpha - 1)/\beta$ and $S_\psi(m_\psi + \mu_{\text{sym}})$ is a smooth function. For $\mu_{\text{sym}} = 0$ (symptoms at the mode), $F_\psi(m_\psi) \approx 0.35$ for large $\psi\alpha$ (the mode of a Gamma is always left of the mean), so roughly 65% of transmission is averted. For $\mu_{\text{sym}} > 0$, the fraction averted decreases smoothly.

**The direction of the $\psi$-effect depends on $\mu_{\text{sym}}$:**

| $\mu_{\text{sym}}$ | Spike ($\psi \to 0$) | Smooth ($\psi \to 1$) | Which benefits more? |
|---|---|---|---|
| $< 0$ (symptoms before peak) | $\theta = 1$ | $\theta \approx 0.65\text{--}0.95$ | Spike (slightly) |
| $= 0$ (symptoms at peak) | $\theta = 0$ (knife-edge) | $\theta \approx 0.65$ | **Smooth** |
| $> 0$ (symptoms after peak) | $\theta = 0$ | $\theta > 0$ (decreasing) | **Smooth** |

When symptoms arrive after peak infectiousness — the most common scenario for many pathogens — **smooth profiles benefit more** from symptom-triggered isolation. The smooth profile still has transmission in its tail; the spike profile has already fired.

### TE for symptom-triggered isolation

$$\text{TE}_{\text{sym}} = p_{\text{sym}} \cdot E_{\delta_{\text{sym}}}[\eta \cdot S_\psi(m_\psi + \delta_{\text{sym}} + \delta_{\text{act}})]$$

For Option A (deterministic $\delta_{\text{sym}} = \mu_{\text{sym}}$, deterministic $\delta_{\text{act}}$):

$$\text{TE}_{\text{sym}} = p_{\text{sym}} \cdot \eta \cdot S_\psi(m_\psi + \mu_{\text{sym}} + \delta_{\text{act}})$$

This is a single evaluation of the Gamma survival function — trivially computable for any $\psi$.

---

## Trigger 3: Contact tracing

### Setup: the cascade of delays

An index case $A$, infected at time $t_A$, is detected via symptom onset at

$$\tau_{\text{id}}^A = l_A + m_\psi + \delta_{\text{sym}}^A$$

(We use symptom-triggered detection as the default for index case identification; other triggers can be substituted.)

Person $A$ infected $\chi_A \sim \text{Poisson}(R_0)$ contacts, with infection attempt $j$ occurring at time $\tau_j = l_A + \varepsilon_j$, where $\varepsilon_j \sim \text{Gamma}(\psi\alpha, \beta)$.

Contact tracers, dispatched upon $A$'s identification, reach contact $j$ after an independent delay $\delta_{\text{trace}}^{(j)}$. The calendar time at which contact $j$ is reached is

$$T_{\text{reach}}^{(j)} = t_A + \tau_{\text{id}}^A + \delta_{\text{trace}}^{(j)}$$

Contact $j$ was infected at calendar time $t_A + \tau_j$, so the time elapsed since contact $j$'s infection when tracing reaches them is

$$\tau_{\text{reach}}^{(j)} = \tau_{\text{id}}^A + \delta_{\text{trace}}^{(j)} - \tau_j = (l_A + m_\psi + \delta_{\text{sym}}^A + \delta_{\text{trace}}^{(j)}) - (l_A + \varepsilon_j)$$

$$= m_\psi + \delta_{\text{sym}}^A + \delta_{\text{trace}}^{(j)} - \varepsilon_j$$

**The index case's onset shift $l_A$ cancels.** This is because $l_A$ shifts both $A$'s symptom time and $A$'s infection attempts by the same amount. What matters is the position of $A$'s jitter $\varepsilon_j$ relative to $A$'s symptom time — and both are measured from $A$'s onset.

### Tracing delay distribution

The delay $\delta_{\text{trace}}^{(j)}$ for each contact is drawn independently from a distribution that captures the logistics of tracing. Natural choices:

- **Exponential:** $\delta_{\text{trace}} \sim \text{Exp}(\lambda_{\text{trace}})$. Memoryless; mean delay $1/\lambda_{\text{trace}}$.
- **Gamma:** $\delta_{\text{trace}} \sim \text{Gamma}(a_{\text{trace}}, b_{\text{trace}})$. Allows a minimum effective delay (when $a_{\text{trace}} > 1$, the mode is positive).
- **Fixed:** $\delta_{\text{trace}} = d_{\text{trace}}$ (deterministic). Simplest case.

The exponential is a reasonable default: some contacts are found quickly, others take longer, and there's no natural minimum delay.

### Intervention on the traced contact

Contact $j$ has their own onset shift $l_j \sim \text{Gamma}((1-\psi)\alpha, \beta)$, independent of $A$'s profile. Isolation begins at

$$\tau_{\text{iso}}^{(j)} = \max(0,\, \tau_{\text{reach}}^{(j)}) + \delta_{\text{act}}$$

The $\max(0, \cdot)$ handles the case where tracing reaches contact $j$ before they were infected ($\tau_{\text{reach}}^{(j)} < 0$), in which case quarantine begins at the time of infection. The fraction of contact $j$'s transmission averted is

$$\theta_j = \eta \cdot S_\psi\left(\tau_{\text{iso}}^{(j)} - l_j\right)$$

Since $S_\psi(x) = 1$ for $x \leq 0$, this correctly gives $\theta_j = \eta$ whenever isolation precedes $j$'s biological onset.

### TE for contact tracing

The fraction of a traced contact's transmission averted, averaged over all sources of randomness, is

$$\text{TE}_{\text{CT}} = p_{\text{trace}} \cdot E\left[\eta \cdot S_\psi\left(\max\left(0,\, m_\psi + \delta_{\text{sym}}^A + \delta_{\text{trace}} - \varepsilon_j\right) + \delta_{\text{act}} - l_j\right)\right]$$

where the expectation is over $\delta_{\text{sym}}^A$, $\delta_{\text{trace}}$, $\varepsilon_j \sim \text{Gamma}(\psi\alpha, \beta)$, $l_j \sim \text{Gamma}((1-\psi)\alpha, \beta)$, and $\delta_{\text{act}}$; and $p_{\text{trace}} \in [0, 1]$ is the probability that a contact is successfully traced.

**Simplification for perfect immediate isolation** ($\eta = 1$, $\delta_{\text{act}} = 0$) **and deterministic** $\delta_{\text{sym}}^A = \mu_{\text{sym}}$:

$$\text{TE}_{\text{CT}} = p_{\text{trace}} \cdot E\left[S_\psi\left(\max\left(0,\, m_\psi + \mu_{\text{sym}} + \delta_{\text{trace}} - \varepsilon_j\right) - l_j\right)\right]$$

The expectation is over three independent random variables: $\delta_{\text{trace}}$, $\varepsilon_j$, and $l_j$.

### How $\psi$ enters: three pathways

**Pathway 1: the "head start" $m_\psi + \mu_{\text{sym}} - \varepsilon_j$.** This is the time from when contact $j$ was infected to when $A$ develops symptoms, minus the tracing delay. It determines how much of a head start contact $j$ has. The distribution of $\varepsilon_j$ depends on $\psi$:

- **Spike ($\psi \to 0$):** $\varepsilon_j \approx 0$ for all $j$. All of $A$'s contacts were infected at nearly the same time ($\approx l_A$), so they all have the same head start $\approx m_\psi + \mu_{\text{sym}} + \delta_{\text{trace}} = \mu_{\text{sym}} + \delta_{\text{trace}}$ (since $m_\psi \to 0$). Tracing either catches them all or misses them all (modulo variation in $\delta_{\text{trace}}$).

- **Smooth ($\psi \to 1$):** $\varepsilon_j$ is broadly distributed. Some contacts were infected early (large $\varepsilon_j$ means late infection relative to onset — wait, $\varepsilon_j$ is the jitter, which is the time from onset to infection attempt; large $\varepsilon_j$ means the contact was infected late by $A$). Late-infected contacts have small or negative head starts and can be caught; early-infected contacts have large head starts and may have already transmitted.

**Pathway 2: contact $j$'s own profile shape.** The survival function $S_\psi$ determines how much of $j$'s transmission falls after $\tau_{\text{iso}}^{(j)}$. Spike profiles give all-or-nothing outcomes (either $j$'s spike hasn't fired, and everything is averted, or it has, and nothing is). Smooth profiles give graded partial benefit.

**Pathway 3: the interaction.** For spike A → spike B: all contacts are infected at similar times, and each contact's outcome is all-or-nothing. The combined effect is **high variance, correlated outcomes**. For smooth A → smooth B: contacts are infected at dispersed times, and each contact's outcome is graded. The combined effect is **lower variance, more predictable** overall reduction.

---

## Summary: how $\psi$ affects each trigger

| Trigger | Key mechanism | Spike ($\psi \to 0$) | Smooth ($\psi \to 1$) |
|---|---|---|---|
| Screening | Window is anchored to mode; all-or-nothing vs graded | Benefits: window covers all transmission if $w > \Delta$ | Costs: substantial pre-window transmission |
| Symptoms ($\mu_{\text{sym}} \geq 0$) | Symptoms lag behind spike | Costs: spike fires before symptoms | Benefits: tail transmission still avertable |
| Symptoms ($\mu_{\text{sym}} < 0$) | Symptoms precede spike | Benefits: all transmission averted | Benefits (less): most transmission averted |
| Contact tracing | Clustered vs dispersed infection times | All contacts have similar head starts | Contacts have variable head starts |

The key insight is that screening and symptom-triggered isolation have **opposite $\psi$-dependence** (when $\mu_{\text{sym}} \geq 0$): screening benefits spike profiles more, while symptom-triggered isolation benefits smooth profiles more. This reversal arises because screening is anchored to the detectability window (which tracks peak infectiousness), while symptoms are anchored to a biological clock that doesn't speed up when the profile narrows.

Contact tracing inherits the $\psi$-dependence of the index case's detection mechanism (pathway 1) and adds the contact's own profile shape (pathway 2). The net effect depends on both.

---

## Extensions and design choices

### Imperfect isolation

Imperfect isolation ($\eta < 1$) scales all fractions averted by $\eta$. It does not change the $\psi$-dependence of the *relative* effectiveness across profiles — but it does change the *absolute* effectiveness, and this interacts with $\psi$ in one important way.

For spike profiles, the outcome is bimodal: either $S_\psi \approx 1$ (spike not yet fired) or $S_\psi \approx 0$ (spike already fired). Multiplying by $\eta < 1$ reduces the "caught before spike" case from $\theta = 1$ to $\theta = \eta$, but the "missed the spike" case remains $\theta = 0$. The impact of $\eta$ is linear.

For smooth profiles, $S_\psi$ takes intermediate values, and the residual transmission $(1 - \eta) \cdot S_\psi$ persists. This means imperfect isolation erodes the benefit more noticeably for smooth profiles (where there is always remaining profile mass to leak through) than for spike profiles (where the remaining mass is typically near 0 or near 1).

### Imperfect temporal adherence

A positive $\delta_{\text{act}}$ shifts the effective isolation time later, reducing $\theta$ for everyone. But the *marginal cost* of delay differs by $\psi$:

$$\frac{\partial \theta}{\partial \delta_{\text{act}}} = -\eta \cdot f_\psi(\tau_{\text{iso}} - l_i)$$

This is the density of the individual profile at the isolation time. For spike profiles, $f_\psi$ is either very large (if isolation is near the spike) or zero (if not). For smooth profiles, $f_\psi$ is moderate over a wide range. So:

- **Spike profiles** are insensitive to small delays (the spike is either caught or not), but very sensitive to delays that push isolation past the spike.
- **Smooth profiles** are uniformly sensitive to delays (every extra hour of delay costs a small but consistent fraction of transmission).

### Generation interval distortion under isolation

Isolation truncates the effective generation interval distribution, biasing surviving transmissions early. This is ψ-dependent:

- **Spike:** If the spike is caught, no transmissions survive (the generation interval distribution is empty). If the spike is missed, the generation interval distribution is unchanged.
- **Smooth:** Transmissions in the tail are preferentially removed, shifting $g_{\text{eff}}(\tau)$ earlier. Via the Euler–Lotka equation, a shorter effective generation interval implies a higher growth rate for a given $R_{\text{eff}}$. This is a **partially self-defeating** effect: the $R_{\text{eff}}$ reduction from isolation is partly offset by the generation interval shortening for smooth profiles.

This effect is specific to smooth profiles and does not arise for spike profiles (where isolation either removes all transmissions or none). It means that the growth-rate impact of isolation is even more favourable for spike profiles than the $R_{\text{eff}}$ comparison alone suggests.

### Combining triggers

In practice, individuals may be detected by whichever trigger fires first. The combined detection time is

$$\tau_{\text{det}} = \min(\tau_{\text{det}}^{\text{screen}}, \tau_{\text{det}}^{\text{sym}})$$

(with $\tau_{\text{det}}^{\text{screen}} = \infty$ if undetected by screening, $\tau_{\text{det}}^{\text{sym}} = \infty$ if asymptomatic). The combined TE is not the sum of the individual TEs (because the second trigger is redundant conditional on the first firing earlier), but is bounded below by $`\max(\text{TE}_{\text{screen}}, \text{TE}_{\text{sym}})`$.

For contact tracing, the index case is detected by one of the direct triggers, and tracing then acts on their contacts. The overall population-level impact combines: (1) direct isolation of index cases, and (2) secondary isolation of traced contacts. Both depend on $\psi$.
