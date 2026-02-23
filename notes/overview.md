# Project plan for "Assessing the impact of individual variation in infectiousness on epidemic dynamics"

## Summary 

Two central quantities - the basic reproduction number ($R_0$) and the generation interval distribution - are critical for defining an epidemic's trajectory. In the renewal equation framework for epidemic modeling, both can be combined into a single "infectiousness trajectory", $A(\tau)$, where $\tau$ represents the time elapsed since infection. Thus, $A(\tau)$ is the "expected contribution to the force of infection by an individual that was itself infected $\tau$ units of time ago" [(Breda *et al.*, 2012)](https://www.tandfonline.com/doi/pdf/10.1080/17513758.2012.716454). Here, $\int_0^\infty A(t) dt = R_0$ and $A(\tau) / R_0$ is equal to the generation interval distribution. In this framework, the force of infection at time $t$ is 

$$ F(t) = \int_0^\infty F(t - \tau) S(t - \tau) A(\tau) d\tau$$

where $S(t)$ is the proportion (density) of susceptible individuals at time $t$. This gives an expression for the incidence: 

$$ J(t) = -\dot{S}(t) = S(t) F(t) = S(t) \int_0^\infty F(t - \tau) S(t - \tau) A(\tau) d\tau = S(t) \int_0^\infty J(t - \tau) d\tau$$ 


An individual-level analog of the basic reproduction number has been defined and analyzed extensively, *e.g.*, in ["Superspreading and the effect of individual variation on disease emergence"](https://www.nature.com/articles/nature04153) (Lloyd-Smith *et al.*, 2005). 