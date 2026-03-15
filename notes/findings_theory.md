## Three canonical individual infectiousness profiles 

In the renewal equation framework for infectious disease modeling, the force of infection is 

$$ F(t) = \int_0^\infty F(t-\tau) S(t-\tau) A(\tau)$$ 

where $F(t)$ is the force of infection at time $t$, $S$ is the density of susceptible individuals, and $A(\tau)$ is the "infectiousness profile", capturing the contribution to the force of infection from an individual who was infected $\tau$ time units ago. $A(\tau)$ is thus the central object that defines the epidemic's deterministic trajectory. 

The infectiousness profile $A(\tau)$ is related to other key epidemiological parameters, like the reproduction number, $R_0$, and the generation interval distribution, $g(\tau)$: 

$$R_0 = \int_0^\infty A(\tau) d\tau$$ 

$$g(\tau) = \frac{A(\tau)}{R_0} $$ 

Importantly, $A(\tau)$ is a population-level object (an expected value across individuals). At the individual level, various "individual infectiousness profiles" $a_i(\tau)$ can yield the same $A(\tau)$. For example, the SEIR model is commonly motivated by assuming a person undergoes an exponentially-distributed latent period (ending at rate $\eta$), followed by another exponentially-distributed latent period (ending at rate $\gamma$), and is infectious at level $\beta$ during the infectious period, i.e., 

$$ a_i(\tau) = \begin{cases}
\beta &\qquad \tau \in [t_i^\text{on}, t_i^\text{on} + t_i^\text{off}] \qquad \text{ where } \tau_i^\text{on} \sim \text{Exp}(\eta) \text{ and } \tau_i^\text{off} \sim \text{Exp}(\gamma) \\ 
0 &\qquad \text{otherwise}
\end{cases}$$ 

It can be shown that this yields the population-level infectiousness profile 

$$A(\tau) = \beta \frac{\eta}{\eta - \gamma}(e^{-\gamma \tau} - e^{-\eta \tau})$$ 



## Shifted Gamma construction: a clean one-parameter interpolation

The original three profiles (smooth, stepwise, spike) are useful theoretical bookends but leave open the question of how to *interpolate* between the smooth and spike extremes with a single punctuation parameter, while preserving:

1. $\int a_i(\tau) d\tau = R_0$ exactly for every individual,
2. $E[a_i(\tau)] = A(\tau)$ exactly (population kernel recovered in expectation),
3. $A(\tau)$ invariant as the punctuation parameter varies,
4. All individual profiles having the same shape, height, and width (differing only in location).

The shifted Gamma construction achieves all four simultaneously.

### Setup

Choose a population-level generation interval distribution

$$A(\tau) = R_0 \cdot \text{Gamma}(\tau; \alpha, \beta), \qquad \beta = \alpha / \mu,$$

where $\mu$ is the mean generation time and $\alpha$ controls the shape (larger = more symmetric/bell-shaped). This does not need to match the SEIR kernel exactly.

### Decomposition

For each person $i$, draw the number of infection attempts: 

$$ \chi_i \sim \text{Poisson}(R_0)$$ 

Now, we need to identify the timing of each infection event. Decompose each attempted infection time as

$$\xi_j = l_i + \varepsilon_j,$$

where:

- $l_i \sim \text{Gamma}(\alpha (1-\psi), \beta)$ is an individual-specific **latent period,** (drawn once per individual), during which the person is not infectious,
- $\varepsilon_j \sim \text{Gamma}(\psi \alpha, \beta)$ is an **additional lag until infection attempt $j$ occurs,** after the end of the latent period (drawn independently per infection attempt: $j \in \\{1, 2, \dots, \chi_i\\}\$),
- $\psi \in (0, 1)$ is the **punctuation parameter**, where $\psi \rightarrow 1$ is a smooth profile identical to the population-level generation interval distribution, and $\psi \rightarrow 0$ is a delta function. 

### Why it works

The key identity is that Gamma distributions with the same rate are closed under convolution:

$$\text{Gamma}(\alpha (1-\psi), \beta) + \text{Gamma}(\alpha \psi, \beta) = \text{Gamma}(\alpha, \beta).$$

So the marginal distribution of $\xi_j$ (integrating over $l_i$) is $\text{Gamma}(\alpha, \beta)$ regardless of $\psi$. This immediately gives $E[a_i(\tau)] = A(\tau)$ and guarantees that $A(\tau)$ is invariant as $\psi$ varies.

### Individual profiles

Each individual's infectiousness profile is

$$a_i(\tau) = R_0 \cdot f_\psi(\tau - l_i), \qquad f_\psi = \text{Gamma}(\alpha \psi, \beta),$$

which is zero for $\tau < l_i$ and follows a $\text{Gamma}(\alpha \psi, \beta)$ density thereafter. Since $f_\psi$ is a proper density, $\int a_i(\tau)\,d\tau = R_0$ exactly. All individuals share the same $f_\psi$ — the profiles are identical in shape, height, and width; only the onset time $l_i$ differs (though, in later sections, we relax this assumption).

### Limiting behaviour

| $\psi$ | Individual profile $f_\psi$ | Shift distribution (of $l_i$) | Interpretation |
|---|---|---|---|
| $\psi \to 0$ | $\delta$-function (infinitely narrow spike) | $\approx \text{Gamma}(\alpha, \beta) = A(\tau)/R_0$ | Maximally punctuated: all infection attempts at one instant |
| Small $\psi$ (e.g. 0.25) | Narrow unimodal bump | Wide shift distribution | Punctuated but not singular |
| $\psi \to 1$ | $\approx \text{Gamma}(\alpha, \beta) = A(\tau)/R_0$ (broad, matches population kernel) | $\delta$-function at 0 (no shift) | Smooth: all individuals have the same profile $\approx A(\tau)$ |

### Variance decomposition

The total variance of a person's attempted infection times decomposes additively:

$$\text{Var}(\xi_j) = \underbrace{\text{Var}(l_i)}_{\text{between-individual}} + \underbrace{\text{Var}(\varepsilon_j)}_{\text{within-individual}} = \frac{\alpha(1-\psi)}{\beta^2} + \frac{\alpha \psi}{\beta^2} = \frac{\alpha}{\beta^2}.$$

The fraction of total variance that is between-individual (i.e., due to punctuation) is

$$\text{punctuation fraction} = \frac{\alpha (1- \psi)}{\alpha} = 1 - \psi.$$

This gives a clean interpretation: $\psi \in (0, 1)$ is the fraction of generation-time variability that is within-individual.

### Design choices

- **$\alpha$**: Controls the shape of $A(\tau)$. Larger values give a more symmetric, bell-shaped kernel.

- **$\mu$**: The mean generation time, equal to $\alpha / \beta$. Changing $\mu$ scales the time axis without affecting the shape of $A(\tau)$ or the punctuation structure.

### Extension: heterogeneous $\psi$ across individuals

The shifted Gamma construction is remarkably robust to heterogeneity in the punctuation parameter. If different individuals draw different $\psi_i$ values from some distribution $\psi_i \sim F$ on $(0, 1)$, the population-level kernel $A(\tau)$ is **exactly unchanged** — and this holds for *any* distribution $F$, with no moment conditions required.

#### Why it works

The argument is simple. For a given individual $i$ with punctuation parameter $\psi_i$, their attempted infection times decompose as:

$$\xi_j = l_i + \varepsilon_j, \qquad l_i \sim \text{Gamma}(\alpha (1 - \psi_i), \beta), \quad \varepsilon_j \sim \text{Gamma}(\alpha \psi_i, \beta)$$

By the Gamma additivity property, the marginal distribution of each attempted infection time is:

$$\xi_j \mid \psi_i \sim \text{Gamma}(\alpha, \beta) \qquad \text{for every } \psi_i \in (0, 1)$$

The crucial point: the right-hand side **does not depend on $\psi_i$**. The identity $\text{Gamma}(\alpha(1 - \psi), \beta) + \text{Gamma}(\alpha \psi, \beta) = \text{Gamma}(\alpha, \beta)$ is exact for every $\psi$, not just in expectation. So each individual's contribution to $A(\tau)$ is $R_0 \cdot \text{Gamma}(\tau; \alpha, \beta)$ regardless of their punctuation parameter.

When we average over the population:

$$A(\tau) = E_{\psi_i}\left[E[a_i(\tau) \mid \psi_i]\right] = E_{\psi_i}\left[R_0 \cdot \text{Gamma}(\tau; \alpha, \beta)\right] = R_0 \cdot \text{Gamma}(\tau; \alpha, \beta)$$

The expectation over $\psi_i$ passes through because the integrand does not depend on $\psi_i$. The distribution $F$ cancels entirely.

#### Example distributions

With $\alpha = 10$, we can define $\psi_i = B_i$ where $B_i \sim \text{Beta}(a, b)$ on $(0, 1)$.

**Example 1: Homogeneous population (degenerate $F$).** $\psi_i = 0.5$ for all $i$ (i.e., $B_i = 0.5$ with probability 1). Everyone has the same moderate punctuation. This is the baseline case.

**Example 2: Bimodal / U-shaped distribution.** $B_i \sim \text{Beta}(0.3, 0.3)$, giving a U-shaped density on $(0, 1)$ with most mass near the extremes. This produces a population split between highly punctuated individuals ($\psi_i \approx 0.05$, narrow spikes) and highly smooth individuals ($\psi_i \approx 0.95$, broad profiles resembling $A(\tau)$), with few individuals in between. Despite this radical heterogeneity in individual profile shapes, $A(\tau)$ is exactly $R_0 \cdot \text{Gamma}(10, \beta)$.

**Example 3: Right-skewed (mostly smooth).** $B_i \sim \text{Beta}(5, 1)$, concentrating $\psi_i$ near 1. Most individuals have smooth, population-like profiles; a small minority are punctuated. $A(\tau)$ is unchanged.

**Example 4: Left-skewed (mostly spiky).** $B_i \sim \text{Beta}(1, 5)$, concentrating $\psi_i$ near 0. Most individuals have narrow spikes; a few are smooth. $A(\tau)$ is still unchanged.

In all four cases — and in any other distribution on $(0, 1)$ — the population-level kernel, the mean generation interval, and the mean-field ODE dynamics are identical. Importantly, heterogeneous $\psi$ alone does not introduce variation in $R_i$: each individual's profile integrates to exactly $R_0$ regardless of their $\psi_i$, because $f_{\psi_i}$ is a proper density for every $\psi_i$. Variation in $R_i$ arises only when there is a time-varying contact process to interact with; what $\psi$ heterogeneity does is create heterogeneous *sensitivity* to that contact variation. The observables that differ across $\psi$ distributions are therefore conditional on the contact environment: the degree of offspring overdispersion given periodic contacts, the distribution of generation intervals from individual infectors, and the effectiveness of timing-dependent interventions.

#### Implication for interventions

This invariance property means that a population with heterogeneous $\psi$ provides a clean setting for studying how interventions perform across a mixture of profile types, without any confounding from changes in $A(\tau)$. For instance, in a U-shaped population (Example 2):

- Detect-and-isolate would be highly effective against the spiky subpopulation (all-or-nothing aversion, biased toward "all") but less effective against the smooth subpopulation (partial aversion, generation interval distortion).
- Gathering size restrictions would primarily affect the spiky subpopulation in absolute variance terms, while having the same mean effect on everyone.
- Contact tracing effectiveness would depend on the $\psi$ of both infector and infectee, creating a $2 \times 2$ structure of spiky/smooth pairings within the same population.

The population-level $R_{\text{eff}}$ under any intervention is a mixture over individual-level $\psi_i$ values, and the mixture weights are set by $F$ — but the uncontrolled $A(\tau)$ is invariant, so any differences in epidemic outcomes under intervention are attributable purely to the $\psi$ heterogeneity, not to differences in the baseline transmission kernel.

## Decomposing infectiousness into biological and contact components

The individual infectiousness profile $a_i(\tau)$ conflates two distinct processes: the biological trajectory of infectiousness (viral load, shedding) and the rate at which the individual makes potentially transmissive contacts. Separating these gives a cleaner framework for studying how temporal structure in contact patterns interacts with punctuated biological infectiousness.

### Definitions

Decompose

$$a_i(\tau) = b_i(\tau) \cdot c_i(t_i + \tau),$$

where $t_i$ is the calendar time at which person $i$ was infected, and:

- **$b_i(\tau)$** is the *biological infectiousness profile*: a probability density on $[0, \infty)$ satisfying $\int_0^\infty b_i(\tau)\,d\tau = 1$. It describes *when* the individual is capable of transmitting, relative to their time of infection. It contains no information about how many people they will infect — only the temporal shape.

- **$c_i(t)$** is the *instantaneous transmission potential* at calendar time $t$: the expected number of secondary infections the individual would produce at time $t$ if all their biological infectiousness were concentrated at that instant (i.e., if $b_i$ were a delta function). This absorbs contact rate, transmission probability per contact, and susceptibility of contacts. When $c_i$ is constant in time, $c_i = \nu_i$ is exactly the individual reproduction number in the sense of Lloyd-Smith *et al.* (2005).

The individual's realized reproduction number is then

$$\nu_i = \int_0^\infty a_i(\tau) d\tau = \int_0^\infty b_i(\tau) c_i(t_i + \tau) d\tau,$$

a weighted average of their contact potential over the times they are biologically infectious.

### Population-level factorization

Define $B(\tau) = E[b_i(\tau)]$ and $C(t) = E[c_i(t)]$. Under the assumption that $b_i$ and $c_i$ are **independent** — that an individual's biological trajectory is independent of their contact pattern — the population-level infectiousness kernel factorizes:

$$A_t(\tau) = E[a_i(\tau)] = B(\tau) \cdot C(t_i + \tau).$$

In the renewal equation, all contributions to the force of infection at calendar time $t$ share the same calendar-time contact factor $C(t)$, so $C$ factors out:

$$J(t) = S(t) \cdot C(t) \int_0^\infty J(t - \tau)\,B(\tau)\,d\tau.$$

The time-varying contact function acts as a multiplicative modulator on the effective reproduction number. When $C(t) = C$ (constant), we recover the standard time-homogeneous renewal equation with $A(\tau) = C \cdot B(\tau)$.

**Note:** If $b_i$ and $c_i$ are not independent — for instance, if symptomatic individuals reduce contacts late in their infection — the factorization breaks and the joint distribution $E[b_i(\tau) \cdot c_i(t)]$ is needed. We assume independence throughout this project.

### Connection to the shifted Gamma construction

In the shifted Gamma framework (Section 9), with constant contacts:

- $b_i(\tau) = f_\psi(\tau - l_i)$ — the $\text{Gamma}(\alpha \psi, \beta)$ density shifted to onset time $l_i$
- $c_i(t) = R_0$ for all $i, t$

This is the simplest case: all the punctuation lives in $b_i$, contacts are homogeneous and time-invariant, and $\nu_i = R_0$ for every individual.

### Punctuated infectiousness as a source of superspreading

The decomposition reveals a mechanism for generating overdispersion in the individual reproduction number $\nu_i$ that, to our knowledge, has not been previously described.

**Setup.** Suppose:

- Biological profiles are shifted Gammas: $b_i(\tau) = f_\psi(\tau - l_i)$ with punctuation parameter $\psi$
- The contact function is the *same for everyone* and varies periodically in calendar time, e.g.: $c_i(t) = c(t) = R_0 \cdot z(t)$, where $z(t)$ oscillates around 1 (e.g., weekday/weekend cycles, seasonal forcing)

There is no individual-level variation in contact rates — the only heterogeneity is in the biological timing $l_i$. Yet the individual reproduction number becomes:

$$\nu_i = R_0 \int_0^\infty f_\psi(\tau - l_i) \cdot z(t_i + \tau)\,d\tau.$$

The variance of $\nu_i$ depends on how much the biological profile "averages over" the oscillations in $z$:

- **Smooth limit** ($\psi \to 1$): $b_i \approx B$ for all $i$, so $\nu_i \approx R_0 \int B(\tau) z(t_i + \tau) d\tau$. The broad biological window averages over $z$, and $\nu_i \approx R_0$ for everyone. Minimal overdispersion.

- **Delta limit** ($\psi \to 0$): $b_i \to \delta(\tau - l_i)$, so $\nu_i \to R_0 \cdot z(t_i + l_i)$. The individual's reproduction number directly *samples* the contact function at a single point. If $z$ oscillates between 0.5 and 1.5, $\nu_i$ varies between $0.5 R_0$ and $1.5 R_0$. Maximal overdispersion.

- **Intermediate $\psi$**: the narrower $b_i$ is, the less it averages over $z$, and the more variable $\nu_i$ becomes. The punctuation parameter $\psi$ continuously controls the degree of overdispersion.

**The mechanism in words:** when biological infectiousness is concentrated in a narrow window, the individual's total transmission depends on *when* that window falls relative to the contact landscape. Some individuals' narrow windows align with high-contact periods (weekdays, social events, crowded settings); others align with low-contact periods (nights, weekends, holidays). This "sampling" of the contact function creates between-individual variation in $\tilde{\nu}_i$ — i.e., superspreading — even without any intrinsic heterogeneity in contact rates or biological susceptibility.

This is distinct from the standard superspreading mechanisms (heterogeneous $\nu_i$ across individuals, heterogeneous susceptibility, stochastic contact networks). It requires only three ingredients:

1. Punctuated biological infectiousness (narrow $b_i$)
2. Stochastic variation in the timing of infectiousness (random $l_i$)
3. Time-varying contact rates at the population level ($z(t) \neq \text{const}$)

All three are empirically plausible for many infectious diseases. The shifted Gamma construction provides a single parameter ($\psi$) that controls the strength of this effect.