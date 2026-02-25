# The impact of punctuated infectiousness on infectious disease dynamics 

## Introduction


## Methods

### Epidemic dynamics with renewal equations

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

### The individual-level infectiousness profile 

The generation interval distribution, $g(\tau)$, describes how secondary infections are distributed over time. Like $R_0$, $g(\tau)$ is a population-level object: it describes the temporal distribution of infection times for a "typical" infectious individual. As with the individual reproduction number, $\nu_i$, where $E_i[\nu_i] = R_0$, we can define an individual generation interval distribution, $\xi_i$, where $E_i[\xi_i(\tau)] = g_i(\tau)$. Furthermore, we can define an individual infectiousness profile, $a_i(\tau) = \nu_i \xi_i(\tau)$, such that $E_i[a_i(\tau)] = A(\tau)$. 

An immediate observation is that it is possible to obtain the same population-level $A(\tau)$ from vastly different forms of $a_i(\tau)$. For example, the following $a_i(\tau)$ all yield the SEIR-equivalent $A(\tau)$: 

$$ a_i(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$ 

$$ a_i(\tau) = \begin{cases} 
\beta & \qquad \tau \in [\tau_i^{\text{on}}, \tau_i^{\text{on}} + \tau_i^{\text{off}}] \\
0 & \qquad \text{otherwise }
\end{cases} \qquad \text{ where } \tau_i^{\text{on}} \sim \text{Exp}(\gamma) \text{ and } \tau_i^{\text{off}} \sim \text{Exp}(\alpha)$$ 

$$ a_i(\tau) = \frac{\beta}{\gamma} \cdot \delta_{t_i}(\tau) \text{ where } t_i \text{ is distributed according to the generation interval density, } A(\tau)/R_0$$ 

Eq XX is simply $A(\tau)$ replicated exactly for each person. Eq XX is the underlying individual-level model often conceptualized when using the SEIR model, where a person spends an exponentially-distributed amount of time in the latent/exposed state and another exponentially-distirbuted amount of time in the infectious state. Eq XX is an extreme case where a person's infectiousness is completely concentrated at a single moment. All three have the same expectation, and thus yield the same average population-level dynamics. 

### The impact of the individual infectiousness profile on uncontrolled epidemic dynamics 

Despite yielding the same average population-level dynamics, epidemics that are governed by different forms of $a_i(\tau)$ may yield different stochastic dynamics. 

### A Gamma convolutional model for the individual infectiousness profile 

To examine the impact of the punctuatedness of $a_i(\tau)$, holding all else equal, we introduce a one-parameter family of functions that allow us to reconstruct a Gamma-distributed $A(\tau)$ using Gamma-distributed $a_i(\tau)$ that are identical except for a time shift. 

Let 

$$
A(\tau) = R_0 \text{Gamma}(\tau; \alpha, r), \qquad r = \alpha / \mu
$$

where $\mu$ is the mean generation interval and $\alpha$ controls the shape of the population-level infectiousness profile. 

### The impact of punctuated infectiousness on test-based screening and isolation countermeasures 


### Splitting the infectiousness profile into biological infectiousness and contacts 


### Assessing overdispersion resulting from variable contact rates and punctuated biological infectiousness 


## Results


## Conclusions