# The impact of punctuated infectiousness on infectious disease dynamics 

## Introduction


## Methods

### The individual-level infectiousness profile 

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




### The impact of the individual infectiousness profile on uncontrolled epidemic dynamics 


### A Gamma convolutional model for the individual infectiousness profile 


### The impact of punctuated infectiousness on test-based screening and isolation countermeasures 


### Splitting the infectiousness profile into biological infectiousness and contacts 


### Assessing overdispersion resulting from variable contact rates and punctuated biological infectiousness 


## Results


## Conclusions