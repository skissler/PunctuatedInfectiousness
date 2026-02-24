# The impact of punctuated infectiousness on infectious disease dynamics 

## Introduction


## Methods

### Population-level disease transmission as a renewal equation

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

The generation interval, $g(\tau)$, is the expected distribution of secondary infection times from an infected individual. Importantly, this is a population average -- the times at which an individual person is expected to transmit may differ significantly from this. Much like the individual reproduction number, $\nu_i$, gives an individual person's expected number of secondary infections (and $E_i[\nu_i] = R_0$), we can define an individual generation interval, $\xi_i(\tau)$, such that $E_i[\xi_i(\tau)] = g(\tau)$. By analog to the population-level model, we can define the *individual infectiousness profile*: 

$$ a_i(\tau) = \nu_i \xi_i(\tau) $$ 

where 

$$ A(\tau) = E_i[a_i(\tau)] $$

It is possible to obtain the same population-level infectiousness profile from vastly different individual infectiousness profiles. For example, we can obtain the same $A(\tau)$ for the SEIR using any of the following formulations: 

$$ a_i(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$ 

$$ a_i(\tau) = \begin{cases} 
\beta & \tau \in [\tau_i^{\text{on}}, \tau_i^{\text{on}} + \tau_i^{\text{off}}] \\
0 & \text{otherwise }
\end{cases} \text{ where } \tau_i^{\text{on}} \sim \text{Exp}(\gamma) \text{ and } \tau_i^{\text{off}} \sim \text{Exp}(\alpha)$$ 

$$ a_i(\tau) = R_0 \cdot \delta_{t_i}(\tau) \text{ where } t_i \text{ is distributed according to the (normalized) generation interval density, } A(\tau)/R_0$$ 





### The impact of the individual infectiousness profile on uncontrolled epidemic dynamics 


### A Gamma convolutional model for the individual infectiousness profile 


### The impact of punctuated infectiousness on test-based screening and isolation countermeasures 


### Splitting the infectiousness profile into biological infectiousness and contacts 


### Assessing overdispersion resulting from variable contact rates and punctuated biological infectiousness 


## Results


## Conclusions