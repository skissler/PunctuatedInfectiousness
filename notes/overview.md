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

- Develop a mathematical framework for describing $A(\tau)$ in terms of the individual generation interavl distribution, $a_i(\tau)$, where $a_i(\tau)$ converges to $A(\tau)$ in expectation in the limit of large population size (i.e., $\limit_{n \rightarrow \infty} \frac{1}{n} \sum_{i=1}^n a_i(\tau) = A(\tau)$)

