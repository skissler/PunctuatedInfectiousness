# Project plan for "Assessing the impact of individual variation in infectiousness on epidemic dynamics"

## Summary 

Two central quantities - the basic reproduction number ($R_0$) and the generation interval distribution - are critical for defining an epidemic's trajectory. In the renewal equation framework for epidemic modeling, both can be combined into a single "infectiousness trajectory", $A(\tau)$, where $\tau$ represents the time elapsed since infection. Thus, $A(\tau)$ is the "expected contribution to the force of infection by an individual that was itself infected $\tau$ units of time ago" [(Breda *et al.*, 2012)](https://www.tandfonline.com/doi/pdf/10.1080/17513758.2012.716454). Here, $\int_0^\infty A(t) dt = R_0$ and $A(\tau) / R_0$ is equal to the generation interval distribution. In this framework, the force of infection at time $t$ is 

$$ F(t) = \int_0^\infty F(t - \tau) S(t - \tau) A(\tau) d\tau$$

where $S(t)$ is the proportion (density) of susceptible individuals at time $t$. This gives an expression for the incidence: 

$$ J(t) = -\dot{S}(t) = S(t) F(t) = S(t) \int_0^\infty F(t - \tau) S(t - \tau) A(\tau) d\tau = S(t) \int_0^\infty J(t - \tau) A(\tau) d\tau$$ 

For some choices of $A(\tau)$, the renewal equation framework is equivalent to more commonly known systems of ordinary differential equations: for example, we obtain the SIR model with 

$$ A(\tau) = \beta e^{-\alpha \tau}$$ 

and we obtain the SEIR model with 

$$ A(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$

where $\alpha$ is the recovery rate (governing I $\rightarrow$ R transitions) and $\gamma$ is the latency rate (governing E $\rightarrow$ I transitions) [(Breda *et al.*, 2012)](https://www.tandfonline.com/doi/pdf/10.1080/17513758.2012.716454). 

An individual-level analog of the basic reproduction number has been defined and analyzed extensively, *e.g.*, in "Superspreading and the effect of individual variation on disease emergence" [(Lloyd-Smith *et al.*, 2005)]((https://www.nature.com/articles/nature04153)). The individual-level generation interval distribution has received less attention. The goal of this project is to examine how the shape of the generation interval distribution impacts disease transmission dynamics, *holding the population-level generation interval distribution and basic reproduction number constant*. 

We will: 

- Develop a mathematical framework for describing $A(\tau)$ in terms of the individual generation interval distribution, $a_i(\tau)$, where $a_i(\tau)$ converges to $A(\tau)$ in expectation in the limit of large population size, *i.e.*, 

$$\lim_{n \rightarrow \infty} \frac{1}{n} \sum_{i=1}^n a_i(\tau) = A(\tau)$$

- Examine how three ways of parameterizing $a_i(\tau)$ yield different stochastic epidemic dynamics (establishment probabilities, final size distributions, etc), even when all converge to the same $A(\tau)$. We will focus on the SEIR model, and the three cases will be: 
	
$$ a_i(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau}) $$ 

$$ a_i(\tau) = \begin{cases} 
\beta & \tau \in [\tau_i^{\text{on}}, \tau_i^{\text{on}} + \tau_i^{\text{off}}] \\
0 & \text{otherwise }
\end{cases} \text{ where } \tau_i^{\text{on}} \sim \text{Exp}(\gamma) \text{ and } \tau_i^{\text{off}} \sim \text{Exp}(\alpha)$$ 

$$ a_i(\tau) = \delta_{t_i}(\tau) \text{ where } t_i \text{ is distributed according to the unnormalized density } f(\tau) = \beta \frac{\gamma}{\gamma - \alpha} (e^{-\alpha \tau} - e^{-\gamma \tau})$$ 

In other words, our three main cases will be (a) the "all-equal case", where all individuals have the same individual generation interval distribution as the population, so there's no between-individual variation; (b) the "stepwise case", where people are infectious at the same level ($\beta$) for different amounts of time, and the timing of the onset and offset of infectiousness are exponentially distributed with rate $\gamma$ and $\alpha$, respectively; and (c) the "delta function", case, where each person's infectiousness is all concentrated at time $t_i$, and that time is distributed according to the population-level generation interval distribution. To my knowledge, it's the "stepwise" case that people usually have in mind when writing down the SEIR model. 

- Develop a function that acts as a continuum between the all-equal case (a) and the delta function case (c), ensuring that everyone has the same total infectiousness (i.e., $$\int_\tau a_i(\tau) = R_0$$ for all individuals $i$). This continuum won't include the stepwise case (b), since people have different amounts of infectiousness in that case. A parameter $k$ will govern how punctuated (i.e., how close to the delta case) the distribution is. 

- Determine how different choices of $k$ impact the success of test-based control measures (where we assume some relationship between detectability and infectiousness) and the phylogenetics of outbreaks. 

- Split $a_i(\tau)$ further into a biological and "contact" component, $a_i(\tau) = b_i(\tau) c(t + \tau)$, where the normal assumption is that $c(t) = 1$. We'll assess how periodic $c(t)$ impacts overdispersion in the number of secondary infections -- i.e., all we need is punctuated infectiousness, stochastic variation in the timing of infectiousness, and variation in contact rate (at the population level, not at the individual level) to generate superspreading. To my knowledge, this is an unexamine source of superspreading. 













