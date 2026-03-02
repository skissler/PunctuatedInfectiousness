# The impact of punctuated infectiousness on infectious disease dynamics 

## Introduction

The trajectory of an epidemic is governed by two key factors: the reproduction number, or the expected number of secondary infections that an infectious persion will produce; and the generation interval distribution, or the distribution of times at which those secondary infections are expected to occur, relative to the index case's infection. 

Both the reproduction number and the generation interval distribution are expected values, averaged over a large population. At the individual level, the individual reproduction number [[x](https://www.nature.com/articles/nature04153)] is the expected number of new infections that a single infectious person will produce, such that $R_0 = \sum_i \nu_i$. Similarly, the generation interval distribution is the average of the times a single person is expected to transmit disease, such that $g(\tau) = \sum_i a_i(\tau)$, where $a_i(\tau)$ is the "individual infectiousness kernel" [[x](https://pubmed.ncbi.nlm.nih.gov/26674948/), [x](https://www.sciencedirect.com/science/article/pii/S0025556406002094)]. 





## Epidemic dynamics with renewal equations

For a generic disease transmission process, the force of infection at time $t$ can be written in terms of a renewal equation: 

$$ F(t) = \int_0^\infty F(t-\tau) S(t-\tau) A(\tau) d\tau $$ 

Here, $S(t)$ is the density of susceptible individuals in the population at time $t$ and $A(\tau)$ is the (population-level) "infectiousness profile", a curve describing the expected contribution to the force of infection from an individual who was infected $\tau$ time units ago. The incidence of disease is 

$$ J(t) = -\dot{S}(t) = S(t) F(t) $$ 

The infectiousness profile $A(\tau)$ is thus the fundamental object that determines how the epidemic unfolds. 

The infectiousness profile is related to the basic reproduction number, $R_0$, and the generation interval, $g(\tau)$: 

$$ R_0 = \int_0^\infty A(\tau) d\tau $$ 

$$ g(\tau) = A(\tau) / R_0 $$ 

Furthermore, various choices of $A(\tau)$ generate disease transmission dynamics that are equivalent to more familiar ordinary differential equation-based models. For example, the SIR model 

$$ 
\begin{align}
&\frac{dS}{dt} = -\beta S I \\
&\frac{dI}{dt} = \beta S I - \alpha I  \\
&\frac{dR}{dt} = \alpha I 
\end{align} 
$$ 

is obtained when 

$$
A(\tau) = \beta e^{-\alpha \tau}
$$

Likewise, the SEIR model 

$$
\begin{align}
&\frac{dS}{dt} = -\beta S I  \\ 
&\frac{dE}{dt} = \beta S I - \gamma E  \\ 
&\frac{dI}{dt} = \gamma E - \alpha I  \\ 
&\frac{dR}{dt} = \alpha I  
\end{align}
$$

is obtained when 

$$ A(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$

## The individual-level infectiousness profile 

The generation interval distribution, $g(\tau)$, describes how secondary infections are distributed over time. Like $R_0$, $g(\tau)$ is a population-level object: it describes the temporal distribution of infection times for a "typical" infectious individual. As with the individual reproduction number, $\nu_i$, where $E_i[\nu_i] = R_0$, we can define an individual generation interval distribution, $\xi_i$, where $E_i[\xi_i(\tau)] = g_i(\tau)$. Furthermore, we can define an individual infectiousness profile, $a_i(\tau) = \nu_i \xi_i(\tau)$, such that $E_i[a_i(\tau)] = A(\tau)$. 

An immediate observation is that it is possible to obtain the same population-level $A(\tau)$ from vastly different forms of $a_i(\tau)$. For example, the following $a_i(\tau)$ all yield the SEIR-equivalent $A(\tau)$: 

$$ a_i(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$ 

$$ a_i(\tau) = \begin{cases} 
\beta & \qquad \tau \in [\tau_i^{\text{on}}, \tau_i^{\text{on}} + \tau_i^{\text{off}}] \\
0 & \qquad \text{otherwise }
\end{cases} \qquad \text{ where } \tau_i^{\text{on}} \sim \text{Exp}(\gamma) \text{ and } \tau_i^{\text{off}} \sim \text{Exp}(\alpha)$$ 

$$ a_i(\tau) = \frac{\beta}{\gamma} \cdot \delta_{t_i}(\tau) \text{ where } t_i \text{ is distributed according to the generation interval density, } A(\tau)/R_0$$ 

Eq XX is simply $A(\tau)$ replicated exactly for each person. Eq XX is the underlying individual-level model often conceptualized when using the SEIR model, where a person spends an exponentially-distributed amount of time in the latent/exposed state and another exponentially-distirbuted amount of time in the infectious state. Eq XX is an extreme case where a person's infectiousness is completely concentrated at a single moment. All three have the same expectation, and thus yield the same average population-level dynamics. 

## A Gamma convolutional model for the individual infectiousness profile 

To examine the impact of the punctuatedness of $a_i(\tau)$, holding all else equal, we introduce a one-parameter family of functions that allow us to reconstruct a Gamma-distributed $A(\tau)$ using Gamma-distributed $a_i(\tau)$ that are identical except for a time shift. 

Let 

$$
A(\tau) = R_0 f_A(\tau) \qquad \text{ where } f_A(\tau) \sim \text{Gamma}(\alpha, r), \qquad r = \alpha / \mu
$$

where $\mu$ is the mean generation interval and $\alpha$ controls the shape of the population-level infectiousness profile. 

Now, let 

$$
a_i(\tau) = R_0 f_a(\tau - s_i) \qquad \text{ where } f_a(\tau) \sim \text{Gamma}(\kappa, r) \qquad \text{ and } s_i \sim \text{Gamma}(\alpha - \kappa, r)
$$

Here, $\kappa \in (0, \alpha)$ is a parameter governing the punctuatedness of $a_i(\tau)$, with smaller $\kappa$ yielding a more concentrated infectiousness profile. Since the sum of two Gamma-distributed random variables with the same rate is also Gamma distributed (with the same rate and with shape equal to the sum of the two component shapes), this formulation is guaranteed to converge to $A(\tau)$ in expectation. 

## Splitting the infectiousness profile into biological infectiousness and contacts 

Previously, we saw that we could decompose the individual infectiousness profile into the individual reproduction number and the individual generation interval (Eq XX). Here, we consider a more general, and more biologically-motivated, decomposition of $a_i(\tau)$ into "biological infectiousness" and "contacts". 

Let 

$$
a_i(\tau) = b_i(\tau) \cdot c_i(t_i + \tau)
$$

Here, $b_i(\tau)$ represents a person's "biological infectiousness profile": a probability deensity on $[0, \infty)$ that describes when the individual is capable of transmitting, relative to their time of infection; and $c_i(t)$ is the instantaneous transmission potential at time $t$, i.e., the expected number of secondary infections that an individual would produce at time $t$ if all their biological infectiousness were concentrated at that instant (i.e., if $b_i$ were a delta function). Note that $c_i(t)$ absorbs the contact rate, the trasmission probability per ocntact, and the susceptibility of contacts. 

This is a general case of Eq XX: when $c_i$ is constant in time, $c_i = \nu_i$ (the individual reproduction number in the sense of Lloyd-Smith *et al.* (2005)) and $b_i(\tau)$ is the individual generation interval distribution, $\xi_i(\tau)$. 

This decomposition feeds back into the population-level model. If we assume that $b_i$ and $c_i$ are independent (*i.e.*, an individual's biological infectiousness is independent of their contact pattern), then 

$$
A_i(\tau) = E[a_i(\tau)] = B(\tau) \cdot C(t_i + \tau)
$$

where, like before, if C is contstant in time, this reduces to the familiar 

$$
A_i(\tau) = g(\tau) R_0
$$

Here, $B(\tau)$ is the population-average biological infectiousness over time and $C(t)$ is the population-average contact pattern, again absorbing contact rates and probability of infection per contacct so that the biological infectiousness profile is a density that integrates to 1. 

In the Gamma construction, we have 

$$ 
b_i(\tau) = f_a(\tau - s_i)
$$ 

and 

$$
c_i(t) = R_0 \text{ for all } i,t
$$

## The impact of the individual infectiousness profile on uncontrolled epidemic dynamics 

Despite yielding the same average population-level dynamics, epidemics that are governed by different forms of $a_i(\tau)$ may yield different stochastic dynamics. 

Superspreading, or overdispersion in the secondary infection distribution, is normally thought to arise from overdispersion in individual infectiousness ($\nu_i$) or overdispersion on contacts ($c_i$). Here, we show that overdispersion can also arise from a third mechanism, the interaction of non-constant contacts and highly punctuated infectiousness profiles. 

## The impact of punctuated infectiousness on outbeak interventions

### Screening and isolation 

### Contact tracing 

## Identifiability of punctuated infectiousness 

While the concentration of the individual infectiousness profile can impact how an epidemic unfolds, we lack frameworks for estimating it. Many approaches exist for estimating the generation interval, but there has not been a systematic approach developed to assess how individual infectiousness varies and to specify where on the continuum we sit -- i.e., how much of the variation in the generation interval distribution is attributable to variation between individuals, and how much to variation in infection timing within individuals. 

# Discussion 

# References 

# Acknowledgements 

# Funding 

# Author contributions 

# Supplementary Information 

## Supplementary Methods


## Supplementary Figures 


## Supplementary Tables 

