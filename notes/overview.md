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

$$ a_i(\tau) = R_0 \cdot \delta_{t_i}(\tau) \text{ where } t_i \text{ is distributed according to the (normalized) generation interval density } A(\tau)/R_0$$ 

In other words, our three main cases will be (a) the "all-equal case", where all individuals have the same individual generation interval distribution as the population, so there's no between-individual variation; (b) the "stepwise case", where people are infectious at the same level ($\beta$) for different amounts of time, and the timing of the onset and offset of infectiousness are exponentially distributed with rate $\gamma$ and $\alpha$, respectively; and (c) the "delta function", case, where each person's infectiousness is all concentrated at time $t_i$, and that time is distributed according to the population-level generation interval distribution. To my knowledge, it's the "stepwise" case that people usually have in mind when writing down the SEIR model. 

- Develop a function that acts as a continuum between the all-equal case (a) and the delta function case (c), ensuring that everyone has the same total infectiousness (i.e., $$\int_\tau a_i(\tau) = R_0$$ for all individuals $i$). This continuum won't include the stepwise case (b), since people have different amounts of infectiousness in that case. A parameter $k$ will govern how punctuated (i.e., how close to the delta case) the distribution is. 

- Determine how different choices of $k$ impact the success of test-based control measures (where we assume some relationship between detectability and infectiousness) and the phylogenetics of outbreaks. 

- Split $a_i(\tau)$ further into a biological and "contact" component, $a_i(\tau) = b_i(\tau) c(t + \tau)$, where the normal assumption is that $c(t) = 1$. We'll assess how periodic $c(t)$ impacts overdispersion in the number of secondary infections -- i.e., all we need is punctuated infectiousness, stochastic variation in the timing of infectiousness, and variation in contact rate (at the population level, not at the individual level) to generate superspreading. To my knowledge, this is an unexamine source of superspreading. 

## Work so far

- The function `sim_stochastic()` generates stochastic simulations of epidemics under the three canonical individual infectiousness profiles 

- The script `run_analysis.R` uses `sim_stochastic()` and other helper functions to simulate epidemics with each of the canonical individual infectiousness profiles and compare summary statistics 

- We've shown that punctuated individual infectiousness yields dispersion in the secondary infection distribution (see `contacts.nb` and associated slides)

## Plan 

### Task 1: examine the impact of punctuated infectiousness on test-based control 

Diagnostic tests support outbreak control by helping identify a person's infectiousness status and take appropriate action. Screening interventions are a good example: a person takes a rapid test before work (regular screening), or beore an event (e.g., traveler screening or screening to attend a sporting event). If the test is positive, the person isolates. Another example is test-based *removal* from isolation: a person in isolation can leave when they test negative. Of course, there are many nuances and extensions of these, but these cover a widely used set of test-based interventions. 

For this part of the project, we need to assume some relationship between viral load and infectiousness. For viral load, we can assume a log-linear viral load increase and decrease, which corresponds to an exponential increase/decrease in actual virus concentrations. We've already done a lot of work with this sort of model with SARS-CoV-2 viral kinetics from the NBA. 

In general, the relationship between viral load and (biological) infectiousness is unknown. Recently, Dylan Morris and colleagues came up with a compelling viral load model that uses a random time shift to vary viral load curves and enable Bayesian fitting [(Morris *et al.*, 2025)](https://arxiv.org/pdf/2507.02884). I'd like to use a similar model for our viral load curves: everyone has the same viral load beyond some threshold, but the timing at which viral load crosses that threshold is variable. That takes care of the timing, I think; and we can assume that infectiousness is a function of when the peak happens, and we can model the delta function through the smooth case... maybe? 

This needs a bit more thought. In the fully smooth deterministic case of infectiousness, I think everybody needs the same viral load; as we move towards the delta case, I think we need some amount of variation in the viral load. 

One way to maybe use Dylan's approach: in the fully deterministic/smooth case, this is accompanied by identical detection profiles. Everybody has the same profile, starting at the same time. This then varies as we change how infectiousness relates to the profile: the peak infectiousness always aligns with the peak of the viral load profile (maybe), but as we move towards the delta case, we get more variation in the timing of the viral load curve, too. That's necessary for our setup, I think, and makes sense. So, the idea would be to keep the shape of the viral kinetic curve fixed, but allow its timing to vary, so that it gives us the appropriate amount of variation we need in the infectiousness curves. So, the idea would be to (a) generate infectiousness profiles, (b) line up (identical) viral kinetic curves with the peak of those profiles, and (c) use those viral kinetic curves to describe detectability. I think this is nice, because it assumes that we have the same viral kinetic curves across all profiles -- all that varies is their timing, and so all we're really changing here is the relationship between infectiousness and viral kinetics, but not the kinetics themselves. 

Then, we'll want to vary the kinetics curves, too, to see how narrower windows of detectability impact what happens as we also vary the length of the window of infectiousness. 

Rather than having a triangle viral kinetics curve -- what if we say that the viral kinetics are identical to the infectiousness curve in the smooth case? It's not generally the case that the viral kinetic will equal the population-level generation interval, so this might be a stretch -- and it'll also impact how easy it is to vary the sharpness of the kinetics. We know what kinetics curves look like. We should just use them. So we can ignore this paragraph. 

More concisely: how do different shapes of viral kinetics triangles (sharper, faster/slower up/down) interplay with different levels of punctuation in infectiousness? Let's imagine that peak VL aligns with peak infectiousness, and all that changes is the relationship between them. I need to think carefully about what the right parameterization is here of both the viral kinetics curves and the infectiousness profiles -- but of course the infectiousness profiles, for now, need to be parameterized in a way that they collapse into $A(\tau)$ in the limit. 


### Task 2: separate the infectiousness profile into biological and contact-based components 

A person's infectiousness is a function of both their biological propensity to spread and the rate/density of their contacts. Separating these cleanly may not be straightforward, but we'll try. 

Let's imagine separating both the population-level infectiousness profile, $A(\tau)$, and the individual infectiousness profiles, $a_i(\tau)$, into biological and contact-based components. For example, 

$$A(\tau) = B(\tau) \cdot C(t + \tau)$$

$$a_i(\tau) = b_i(\tau) \cdot c_i(t + \tau)$$

There are lots of pitfalls in how one describes these components, but I think I'm converging on a decent option. 

For the sake of this model, we're assuming an *infinite* population (for the stochastic simulations, of course, we'll look at how finite populations impact the dynamics -- but for the theory, the infinite population assumption is critical). In this case, $c_i$ might represent something like a person's "reach" into that infinite population. Specifically: if a person were somehow "infinitely infectious" at time $t$, then $c_i(t)$ is the expected number of people that that person would (attempt to) infect at that moment. It's kind of nice that the contact function $c$ is defined with respect to the delta function biological infectiousness profile. 

In this framework, we'd actually have $b_i$ integrate to 1 -- or, possibly, to the proportion of "effective contacts" that a person is expected to infect. That proportionality constant can be absorbed into either c or b without loss of generality I think; but given my intuition around epidemic models, I think it might actually make sense to define $c_i$ in terms of some baseline contact rate, and $b$ to integrate to some probability of infection per effective contact. There are various choices that could be made here though, and we might adjust this as we go. 

Does this give us what we want? Let's imagine other biological infectiousness profiles -- all of which now integrate to 1. $b_i \cdot c_i$ should give us the right $a_i$ for the delta function case. For the continuous case, keeping that same $c_i$, we should still get a curve that integrates to the right $R_0$ -- right? 

And then, we need to think about how the individual $c_i$ combine into the larger C. This, I think, is just another population-wide, point-by-point average. We need to confirm this, but that should be the case. This should allow us to account for people with different individual contact patterns, and for temporally varying population-wide mean contact patterns. Of course, there's some weirdness if we really try to account for the fact that peoples' contact patterns aren't independent, but I think in the limit of an infinte population, we can side-step this, and for the sake of a clean model and good theory, that seems worthwhile. 

I think another way of thinking about it is that $c_i$ represents a person's individual $R_0$ ($\nu_i$, in the notation of Lloyd-Smith, I think) *at a given point in time*, i.e., if all their infectiousness were concentrated at a single point, then that's what the individual $R_0$ would be. That then absorbs the probability of transmission per contact into the contact function, which is fine; maybe that's what we want, for better consistency with Lloyd-Smith. 

I thnk that all checks out. 

































