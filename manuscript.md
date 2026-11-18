---
title: "Aerosol Transmission Governs Mitigation Efficacy: Model Misspecification, Identifiability, and the Limits of SIR-Based Pandemic Policy"
---

**Aleck H. Alexopoulos**

CERTH, 6th km Harilaou-Thermi road, Box 60361, Thermi, Thessaloniki 57001, Greece

E-mail: aleck@certh.gr

## Abstract

**Background:** During the COVID-19 pandemic, compartmental epidemic models such as the Susceptible-Infected-Recovered (SIR) model provided critical projections and were central to evaluating mitigation strategies. However, the SIR family does not account for transmission via aerosols accumulating in indoor spaces. In a companion paper, a first-principles Aerosol-SIR (A-SIR) model was derived from population balance equations, introducing infected indoor spaces as a distinct transmission compartment.

**Methods:** The A-SIR model is extended to segregated spaces (outdoors, home, workplace) within a consistent model hierarchy: a single-compartment SIR (SIR-1) used as the literature reference, a three-space SIR (SIR-3), and a three-space A-SIR (A-SIR-3). Aerosol scenarios (0–90% of transmission) are calibrated to an identical early epidemic growth rate, so that all scenarios present the same observable epidemic. Five mitigation measures are evaluated under poor (winter) baseline ventilation. SIR-1 is then fitted to synthetic A-SIR-3 epidemic data to quantify the error in its predicted mitigation effectiveness.

**Results:** SIR-3 reproduces SIR-1 to within 0.07 percentage points, confirming the consistency of the hierarchy. The effectiveness ranking of mitigation measures inverts with the aerosol fraction: at 90% aerosol transmission (R₀ = 3), barriers and distancing reduce the epidemic peak by only 5–6%, whereas restoring proper ventilation suppresses the epidemic. Work-from-home raises the peak by 23% at 50% and by 68% at 90% aerosol transmission, because longer home occupancy increases household aerosol exposure. SIR-1 fitted to aerosol-driven data reproduces the epidemic curve to within the observation noise and recovers the correct R₀, yet overestimates the effectiveness of barriers by up to 85 percentage points and assigns no benefit to ventilation, which in truth suppresses the epidemic. Under realistic ventilation, SIR and A-SIR fits are statistically indistinguishable.

**Conclusions:** The aerosol fraction of transmission can be a governing parameter for mitigation success, and it is not identifiable from epidemic curves. A good fit of an SIR model is therefore no evidence that its intervention predictions are correct. The transmission route must be established from environmental and physical evidence, and ventilation should be treated as a primary intervention for airborne pathogens.

**Keywords:** aerosol transmission, SIR model, A-SIR model, mitigation, ventilation, model misspecification, identifiability, pandemic preparedness

## 1. Introduction

During the COVID-19 pandemic, compartmental population models, notably the Susceptible-Infected-Recovered (SIR) model and its variants, provided critical projections and were used to evaluate mitigation strategies, including measures with major social and economic costs such as lockdowns. The pandemic also revealed fundamental limitations in these models (Holmdahl and Buckee, 2020; Ioannidis et al., 2022; Jit et al., 2023): they were repeatedly re-tuned as their predictions diverged from observations, which compromised their long-term predictive value and generated uncertainty in policy guidance.

A critical and largely overlooked source of this divergence was the aerosol transmission of SARS-CoV-2. Despite early evidence that aerosol transmission was a major route of COVID-19 spread (Morawska and Milton, 2020; Lewis, 2020), the dominant modelling framework, the SIR model, has no mechanism to represent it. It is now recognised that respiratory viruses are carried by a continuum of particle sizes, and that the historical dichotomy between "droplet" and "airborne" transmission is not physically meaningful; transmission is better described in terms of a near-field component, which depends on proximity to the infected person, and a far-field component arising from the accumulation of aerosols in shared indoor air (Wang et al., 2021). The SIR model describes transmission through a single term βSI, which captures near-field transmission between infected and susceptible individuals but cannot represent the accumulation of viral aerosols in shared indoor spaces.

The consequence is not merely qualitative. If a disease transmits partly through aerosols and an SIR model is used to evaluate interventions, the model assigns to its single transmission parameter β a value that absorbs both near-field and aerosol contributions. The resulting β is physically misidentified, and any mitigation that acts only on the near-field pathway will be predicted to be more effective than it actually is.

The most direct prior attempt to couple airborne transmission with compartmental epidemic dynamics is that of Noakes et al. (2006), who embedded the Wells–Riley airborne infection equation in an SEIR-type framework for a single, well-mixed enclosed space with a fixed occupant population. In that approach the airborne infection risk is evaluated from a quasi-steady quanta concentration that responds instantaneously to the current number of infectors; the contamination state of the space is not itself a dynamic variable, and the analysis is confined to one enclosure. The A-SIR model differs in three respects: (i) indoor spaces are treated as a population with their own infected/uninfected state; (ii) space contamination persists and decays with its own time scale (1/ζ), independently of whether infected occupants are still present; and (iii) the formulation follows from population balance equations, so it extends naturally to many spaces at population scale and to segregated space types.

In a companion paper (Alexopoulos, 2025), a general first-principles formulation based on multivariate population balance equations (PBEs) was presented, leading to the A-SIR model. The present paper extends this work in three directions: (i) a consistent hierarchy of segregated models (outdoors, home, workplace) that remains directly comparable with the single-compartment SIR models of the literature; (ii) the evaluation of common mitigation measures as a function of the aerosol fraction of transmission; and (iii) the consequences of model misspecification, i.e. of fitting an SIR model to an aerosol-driven epidemic. The third direction is motivated by the COVID-19 experience. Despite mounting evidence for aerosol transmission from 2020, the WHO was slow to acknowledge airborne transmission, notwithstanding sustained advocacy by an international group of researchers led by Lidia Morawska (Morawska et al., 2023). The present analysis shows why epidemic curves alone could not have settled this question, and why the answer nevertheless determines which interventions work.

## 2. Model and Methods

### 2.1 The A-SIR model

The A-SIR model, derived in Alexopoulos (2025), is given by:

$$\frac{dS}{dt} = -\beta S I - \delta S W \qquad (1)$$

$$\frac{dI}{dt} = \beta S I - \gamma I + \delta S W \qquad (2)$$

$$\frac{dZ}{dt} = -\varepsilon Z I + \zeta W \qquad (3)$$

where S, I and R = 1 − S − I are the fractions of susceptible, infected and recovered humans, and W and Z = 1 − W are the fractions of infected and uninfected indoor spaces. β (day⁻¹) is the near-field transmission rate (dominated by large droplets, with short-range exposure in the exhaled plume) and γ the recovery rate; δ, ε and ζ are the aerosol-pathway parameters. The term δSW represents infection of susceptible humans by contaminated indoor spaces, εZI the contamination of clean spaces by infected humans, and ζW the removal of contamination from indoor spaces by ventilation, deposition and natural inactivation. The model is shown schematically in Fig. 1.

![Figure 1. Aerosol transmission in the reduced population compartment model. S = susceptible, I = infected, R = recovered humans; Z = uninfected indoor space, W = infected indoor space.](paper_out/fig1_schematic.png){width=55%}

A key physical distinction underlies the model: near-field, acute transmission (large droplets and short-range plume exposure, captured by βSI) versus far-field, chronic transmission (aerosol accumulation in shared indoor air, captured by δSW) (Wang et al., 2021). Barriers and distancing interrupt the near-field pathway but have no effect on the room-level aerosol concentration; they therefore do not affect δ. In Table 1, ventilation is assigned to ζ only; its additional effect of lowering the probability that an occupied space becomes contaminated (i.e. on ε) is neglected, so the ventilation benefits reported here are conservative (Annex-2).

### 2.2 Model hierarchy and segregated spaces

Three models are distinguished. **SIR-1** is the classical single-compartment SIR model with effective parameters β~eff~ and γ; it is the model used throughout the literature and serves as the reference. **SIR-3** divides the population between three space types, outdoors (i = 0), home (i = 1) and workplace (i = 2), with near-field transmission only. **A-SIR-3** adds the infected-space compartments W₁ and W₂ for the two indoor space types; outdoors there is no W, since aerosols disperse rapidly. The full equations are given in Annex-1.

Humans move between the outdoors and each indoor space type with exchange rates q~i,in~ and q~i,out~ that are back-calculated from a baseline time budget of 12 h at home (f₁ = 12/24), 10 h at work and other indoor activities (f₂ = 10/24) and 2 h outdoors (f₀ = 2/24). Assuming one visit per day, the mean residence time in space i is τ~i~ = f~i~ (days), so that q~i,out~ = 1/τ~i~, and the steady-state condition q~i,in~S₀ = q~i,out~S~i~ gives q~i,in~ = 1/f₀. Mitigation measures that change the time budget (work-from-home, outdoor time shifts) change the exchange rates and hence the occupancy of each space.

Within each space, near-field transmission is frequency-dependent, β~i~ S~i~ I~i~ / N~i~, where N~i~ = S~i~ + I~i~ + R~i~ is the current occupancy: the infection rate of a susceptible person depends on the fraction of infected people present, not on the share of the population that happens to be in that space type. Indoor spaces share a common rate β~in~, while outdoors the rate is reduced, β₀ = κβ~in~ with κ = 0.2, because convective air currents remove exhaled air outdoors far faster than indoor air movement does.

Because human exchange is fast compared with the epidemic (q~i,in~ = 12 day⁻¹ and q~i,out~ = 2.0–2.4 day⁻¹, versus γ = 0.1 day⁻¹), the occupancy of each space remains close to its quasi-steady value, S~i~ ≈ f~i~ S and I~i~ ≈ f~i~ I. Summing the SIR-3 equations then yields SIR-1 with

$$\beta_{eff} = \sum_{i=0}^{2} f_i \beta_i = \beta_{in}\,(f_1 + f_2 + \kappa f_0) \qquad (4)$$

up to corrections of order γ/q. The indoor rate β~in~ is therefore set from Eq. (4) so that SIR-3 reproduces SIR-1 for a given β~eff~; this is verified numerically in Section 3.1. This consistency is what allows the A-SIR-3 results to be compared directly with literature SIR results expressed in terms of β~eff~ and R₀ = β~eff~/γ.

### 2.3 Ventilation states

The space clearance rate ζ is estimated from aerosol physics in Annex-2. For proper ventilation (0.5 air changes per hour, ACH, at home and 2 ACH at work), clearance of a contaminated space takes hours, giving ζ ≈ 14 day⁻¹ (home) and ζ ≈ 29 day⁻¹ (work). In winter, windows are kept closed and removal relies largely on deposition and on viral inactivation, which is slower in cold, dry air. The **baseline** used in this study is therefore poor (winter) ventilation, with clearance rates ten times lower: ζ₁ = 1.4 day⁻¹ and ζ₂ = 2.9 day⁻¹. The ventilation measure (ζ × 10, Table 1) then corresponds to restoring proper ventilation. A proper-ventilation baseline is examined as a sensitivity case. The aerosol emission rates are ε₁ = 0.5 day⁻¹ (home) and ε₂ = 0.8 day⁻¹ (work).

### 2.4 Mitigation measures

Five COVID-19 mitigation measures are modelled as multiplicative modifiers on the base parameters (Table 1). The modifiers are approximate estimates; they primarily affect the magnitudes of the results, whereas the sign of the work-from-home effect depends on the balance between contact reduction and increased home occupancy (Section 4.3). Vaccination is not examined in this study.

**Table 1.** Mitigation parameter multipliers. γ is unaffected by all interventions.

| Measure | β | δ | ε | ζ | Hours (out/home/work) | Physical basis |
|-------------|-----|-----|-----|-----|----------------|-----------------------------------|
| Barriers | ×0.4 | ×1.0 | ×1.0 | ×1.0 | 2/12/10 | Blocks near-field droplets and plume; no effect on room aerosol |
| Distancing | ×0.5 | ×1.0 | ×1.0 | ×1.0 | 2/12/10 | Reduces near-field exposure; room aerosol is distance-independent |
| Ventilation | ×1.0 | ×1.0 | ×1.0 | ×10 | 2/12/10 | Restores proper aerosol removal; no effect on near-field transmission |
| Masking | ×0.1 | ×0.4 | ×0.6 | ×1.0 | 2/12/10 | Strong source control of droplets; inhalation protection (δ) exceeds exhalation source control (ε) |
| Work-from-home | ×0.2 | ×1.0 | ×1.0 | ×1.0 | 1/23/0 | 80% fewer person-to-person contacts; occupancy shifts to the home |

Work-from-home is implemented through the time budget rather than through a separate emission multiplier. Moving from 12 to 23 hours at home raises the fraction of infected people present at home, I₁/I, from 12/24 to 23/24, i.e. by a factor of 1.92, and the contamination of homes (ε₁Z₁I₁) increases accordingly. The increased household aerosol load is therefore a consequence of the occupancy shift and is not imposed separately.

### 2.5 Aerosol scenarios and calibration

Four scenarios are defined by the fraction f~a~ of transmission attributable to aerosols: 0%, 20%, 50% and 90%. For a given R₀ = β~eff~/γ (3 or 5, with γ = 0.1 day⁻¹), the near-field rate is set to β~eff~(1 − f~a~), and δ is determined by root-finding so that the early exponential growth rate r of A-SIR-3 equals that of the SIR-1 reference, r = β~eff~ − γ. The growth rate is the dominant eigenvalue of the system linearised about the disease-free state. With this calibration, the aerosol share of the force of infection during the exponential phase equals f~a~ exactly, and all scenarios present the same observable early epidemic. Matching the growth rate, rather than a nominal transmission rate, is essential: an earlier scaling based on a nominal rate produced baseline epidemics of very different size across scenarios.

The basic reproduction number of A-SIR-3 is computed with the next-generation method (van den Driessche and Watmough, 2002), treating new human infections (by either pathway) as new infections and space contamination as a transition. For a single compartment this gives

$$R_0 = \frac{\beta}{\gamma} + \frac{\delta\,\varepsilon}{\gamma\,\zeta} \qquad (5)$$

The calibrated values are listed in Table 2. With realistic clearance rates, matching the growth rate also keeps R₀ close to its nominal value (Table 2), because contaminated spaces clear within days.

**Table 2.** Calibrated aerosol infection rate δ and next-generation R₀ of A-SIR-3 (winter baseline ventilation).

| Nominal R₀ | Aerosol fraction | δ (day⁻¹) | r (day⁻¹) | R₀ (A-SIR-3) |
|---|---|---|---|---|
| 3 | 0% | 0 | 0.20 | 3.00 |
| 3 | 20% | 0.486 | 0.20 | 3.07 |
| 3 | 50% | 1.215 | 0.20 | 3.17 |
| 3 | 90% | 2.184 | 0.20 | 3.30 |
| 5 | 0% | 0 | 0.40 | 5.00 |
| 5 | 20% | 0.890 | 0.40 | 5.22 |
| 5 | 50% | 2.228 | 0.40 | 5.56 |
| 5 | 90% | 4.004 | 0.40 | 6.01 |

All simulations start from an infected fraction of 1%, distributed between spaces and between the I and W compartments along the dominant growth eigenvector of the unmitigated system, so that no artificial transient is introduced. Mitigation measures apply from t = 0.

### 2.6 Fitting experiment

To quantify the consequences of model misspecification, synthetic epidemic data are generated with A-SIR-3 (R₀ = 3, winter baseline, each aerosol scenario): the total infected fraction is sampled daily over the first 60 days, which span the epidemic peak, and Gaussian observation noise with σ = 0.2% of the population is added. SIR-1 is fitted to these data by least squares (β and γ free, with several starting points), and the fitted model is then used to predict the peak reduction achieved by each mitigation measure, applying the β multipliers of Table 1 (SIR-1 has no aerosol parameters, so ventilation has no effect in it). These predictions are compared with the true peak reductions from A-SIR-3.

To test whether aerosol transmission can be detected from such data, single-compartment SIR and A-SIR fits (ε and ζ fixed at their true values, δ free) are compared using the Akaike Information Criterion, AIC = n ln(RSS/n) + 2k (Akaike, 1974), where n is the number of observations, RSS the residual sum of squares and k the number of free parameters; a difference ΔAIC = AIC~SIR~ − AIC~A-SIR~ > 10 is taken as strong evidence for A-SIR (Burnham and Anderson, 2002). Because A-SIR with δ = 0 reduces to SIR, the models are nested and ΔAIC ≥ −2 when both fits are converged.

### 2.7 Outcome metrics

Four metrics are reported: (i) peak infected fraction, I~max~, and its reduction relative to the unmitigated epidemic; (ii) final attack rate, the cumulative fraction of the population infected; (iii) epidemic duration, the time during which the infected fraction exceeds the initial prevalence of 1%; and (iv) the share of all infections that occurred through the aerosol pathway. A peak reduction of about 97–98% means that the infected fraction never rises above its initial value of 1%, i.e. complete suppression.

### 2.8 Numerical methods

All ODEs were solved with the Radau implicit solver (scipy.integrate.solve~ivp~; Python 3.13, SciPy 1.18) with a maximum step size of 0.5 days and a relative tolerance of 10⁻⁸. Simulations ran to 360 days. Calibration used Brent's method (scipy.optimize.brentq) on the dominant eigenvalue of the linearised system; fits used scipy.optimize.least~squares~. The code is available from the author.

## 3. Results

### 3.1 Validation

Figure 2 reproduces the single-compartment comparison of the companion paper (Alexopoulos, 2025). With β = 0.4, γ = 0.2, δ = 0.2, ε = 0.8 and ζ = 0.05 day⁻¹, A-SIR gives a peak infected fraction of 37.0% at day 11.6, matching the published values. At β/γ = 0.5, below the SIR epidemic threshold, SIR predicts monotonic decay, whereas A-SIR still produces an epidemic peaking at 27.0% because of the aerosol pathway.

![Figure 2. Single-compartment SIR vs A-SIR (δ = 0.2, ε = 0.8, ζ = 0.05 day⁻¹; I(0) = W(0) = 1%). Left: β/γ = 2.0. Right: β/γ = 0.5, below the SIR epidemic threshold.](paper_out/fig2_validation.png){width=100%}

Figure 3 verifies the model hierarchy of Section 2.2. With δ = 0, SIR-3 reproduces SIR-1 for both R₀ = 3 and R₀ = 5: the peaks agree (30.39% vs 30.38% and 48.01% vs 48.01%), and the infected fractions never differ by more than 0.07 percentage points. All differences between A-SIR-3 and SIR-1 reported below are therefore due to the aerosol pathway and not to the spatial segregation itself.

![Figure 3. SIR-3 (three spaces, frequency-dependent mixing) vs SIR-1 with β~eff~ from Eq. (4). Left: R₀ = 3. Right: R₀ = 5.](paper_out/fig3_sir3_vs_sir1.png){width=100%}

### 3.2 Mitigation effectiveness across aerosol scenarios

Figure 4 shows the infected fraction under each mitigation measure, and Table 3 lists the peak reductions. Because all scenarios are calibrated to the same growth rate, the unmitigated epidemics are nearly identical (peaks of 30.4–31.2% for R₀ = 3 and 48.0–49.7% for R₀ = 5), so differences between scenarios reflect the transmission route alone. The mitigation hierarchy inverts as the aerosol fraction increases.

![Figure 4. Infected fraction under five mitigation measures for R₀ = 3 (top) and R₀ = 5 (bottom), at 0%, 50% and 90% aerosol transmission (winter baseline ventilation).](paper_out/fig4_mitigation_grid.png){width=100%}

**Table 3.** Reduction of the epidemic peak (%) relative to no mitigation; negative values indicate a higher peak. The unmitigated peak (%) is given in the third column.

| R₀ | Aerosol | Peak, none | Barriers | Distancing | Ventilation | Masking | Work-from-home |
|---|---|---|---|---|---|---|---|
| 3 | 0% | 30.4 | 92.4 | 77.0 | 0.0 | 96.7 | 96.7 |
| 3 | 20% | 30.6 | 70.8 | 57.1 | 23.9 | 96.7 | 48.9 |
| 3 | 50% | 30.9 | 38.5 | 31.3 | 67.7 | 96.8 | −23.2 |
| 3 | 90% | 31.2 | 6.4 | 5.3 | 96.8 | 96.8 | −67.9 |
| 5 | 0% | 48.0 | 67.0 | 50.5 | 0.0 | 97.9 | 97.9 |
| 5 | 20% | 48.4 | 44.1 | 34.2 | 14.1 | 97.9 | 28.8 |
| 5 | 50% | 49.0 | 21.5 | 17.3 | 43.0 | 97.9 | −11.8 |
| 5 | 90% | 49.7 | 3.3 | 2.8 | 97.8 | 90.9 | −32.9 |

**Barriers and distancing.** At 0% aerosol transmission, barriers reduce the peak by 92.4% and distancing by 77.0% (R₀ = 3). At 90% aerosol, the same measures reduce it by only 6.4% and 5.3%. Both act only on β and are therefore structurally equivalent in the model; they differ only through the magnitude of their multipliers.

**Ventilation.** Ventilation has no effect at 0% aerosol, reduces the peak by 67.7% at 50% aerosol and suppresses the epidemic completely at 90% aerosol (R₀ = 3). At R₀ = 5 it remains the most effective measure at 90% aerosol (97.8%), ahead of masking (90.9%).

**Masking.** Masking is the most robust measure across all aerosol fractions because it acts simultaneously on β (×0.1), δ (×0.4) and ε (×0.6), i.e. on both transmission pathways. With the multipliers of Table 1 it suppresses the epidemic in almost all scenarios. Although the population-level effectiveness of masking remains debated, substantial evidence supports it (Chu et al., 2020; Crespo et al., 2025).

**Work-from-home backfire.** At 0% aerosol, work-from-home suppresses the epidemic, as an SIR analysis would predict. At 20% aerosol its benefit is already halved, and at 50% and 90% aerosol it raises the peak by 23.2% and 67.9% (R₀ = 3); at 90% aerosol the peak reaches 52.4% at day 14, against 31.2% at day 27 without mitigation. The mechanism is the occupancy shift of Section 2.4: infected individuals spending 23 rather than 12 hours at home contaminate the home environment more, and household members are exposed for longer. This prediction is consistent with the high household secondary attack rates reported for SARS-CoV-2, which exceeded those in other close-contact settings (Madewell et al., 2020). It cannot be obtained from a homogeneous SIR model with a single transmission term.

**Sensitivity to baseline ventilation.** With proper ventilation as the baseline (and the ventilation measure improving it a further tenfold), the results are qualitatively and largely quantitatively unchanged: at 90% aerosol and R₀ = 3, barriers reduce the peak by 7.3%, distancing by 6.1%, ventilation and masking by 96.7%, and work-from-home raises it by 82.6%. The reason is that, when ζ is much larger than γ, the infected-space fraction remains in quasi-steady state, W ≈ εI/ζ, and once δ is calibrated to the observed growth rate, the absolute value of ζ drops out of the baseline dynamics. The conclusions are therefore robust to the uncertain absolute clearance rate; what matters is the relative change in ventilation produced by an intervention.

### 3.3 Outcome metrics

Figure 5 compares the peak, the final attack rate and the epidemic duration for R₀ = 3. At 0% aerosol, barriers and distancing not only lower the peak but also reduce the attack rate (to 33.7% and 58.4%), at the cost of a longer epidemic (166 and 143 days above 1% prevalence, against 77 days without mitigation). At 90% aerosol they leave both the peak and the attack rate (93.3–93.5%) essentially unchanged. At 50% aerosol, ventilation lowers the peak from 30.9% to 10.0% and the attack rate from 93.9% to 67.3%, while extending the epidemic from 77 to 127 days: the epidemic is flattened, spreading the healthcare load over a longer period.

![Figure 5. Peak infected fraction, final attack rate and epidemic duration (time above 1% prevalence) by mitigation measure and aerosol fraction (R₀ = 3, winter baseline). Suppressed epidemics have zero duration.](paper_out/fig5_metrics.png){width=100%}

The share of infections that occur through the aerosol pathway is larger than the calibrated early-phase share. In the 50% scenario, 51.0% of all infections are aerosol-borne without mitigation, but 72.9% under barriers and 91.8% under work-from-home, as measures that act only on near-field transmission shift the epidemic towards the aerosol route.

### 3.4 Aerosol-specific strategies

Figure 6 evaluates strategies aimed specifically at the aerosol pathway (R₀ = 3): a summer outdoor shift (+2 h outdoors, taken equally from home and work), improved home ventilation (ζ₁ × 5), restored workplace ventilation (ζ₂ × 10), both, and the outdoor shift combined with both. At 50% aerosol, these reduce the peak by 14.3%, 33.8%, 19.4%, 61.5% and 70.3%, respectively; at 90% aerosol by 18.4%, 64.3%, 33.8%, 96.8% and 96.8%.

![Figure 6. Aerosol-specific strategies at 50% (left) and 90% (right) aerosol transmission (R₀ = 3, winter baseline). Outdoor shift: f₀ from 2 h to 4 h; home ventilation: ζ₁ × 5; work ventilation: ζ₂ × 10.](paper_out/fig6_strategies.png){width=100%}

Under poor winter ventilation, improving home ventilation is more effective than improving workplace ventilation, even though the workplace has the higher emission rate (ε₂ = 0.8 vs ε₁ = 0.5 day⁻¹). The home is the higher-risk space at baseline because people spend more time there and its clearance rate is lower; the aerosol infection rate in space i scales with ε~i~ f~i~²/ζ~i~ (contamination grows with occupancy, and so does exposure), which at baseline is almost twice as large at home as at work. Targeting the space with the highest aerosol load is the more efficient intervention, and in winter this is the home.

### 3.5 Model misspecification: the identifiability gap

Figure 7 shows the result of fitting SIR-1 to synthetic A-SIR-3 data. For every aerosol fraction, the SIR-1 fit reproduces the data to within the observation noise (root-mean-square residuals of 0.19–0.24%, against σ = 0.2%), and the fitted parameters are close to their effective values (R̂₀ = 3.00–3.05; β̂ = 0.300, γ̂ = 0.098–0.100 day⁻¹). Nothing in the fit signals that the model is misspecified.

![Figure 7. SIR-1 fitted to synthetic A-SIR-3 data (R₀ = 3, 60-day window, σ = 0.2%). Top: data, true trajectory and SIR-1 fit for each aerosol fraction. Bottom: SIR-predicted minus true peak reduction for each mitigation measure; positive values mean that SIR overestimates effectiveness.](paper_out/fig7_identifiability.png){width=100%}

The mitigation predictions of the fitted model are nevertheless badly wrong (Table 4). At 90% aerosol transmission, SIR-1 predicts that barriers reduce the peak by 91.7% and distancing by 76.0%, whereas the true reductions are 6.4% and 5.3%, overestimates by a factor of about 14. SIR-1 assigns no benefit to ventilation, which in truth suppresses the epidemic. For work-from-home, SIR-1 predicts suppression, whereas in truth the peak rises by 67.9%. Already at 50% aerosol, SIR-1 overestimates the effectiveness of barriers and distancing by a factor of 2.4. Only for masking, which acts on both pathways, are the SIR-1 predictions correct.

**Table 4.** Peak reduction (%) predicted by the fitted SIR-1 versus the true A-SIR-3 value (R₀ = 3). Gap in percentage points (pp).

| Aerosol | Measure | SIR-1 predicted | True (A-SIR-3) | Gap (pp) |
|---|---|---|---|---|
| 50% | Barriers | 92.0 | 38.5 | +53.5 |
| 50% | Distancing | 76.4 | 31.3 | +45.1 |
| 50% | Ventilation | 0.0 | 67.7 | −67.7 |
| 50% | Masking | 96.8 | 96.8 | 0.0 |
| 50% | Work-from-home | 96.8 | −23.2 | +119.9 |
| 90% | Barriers | 91.7 | 6.4 | +85.3 |
| 90% | Distancing | 76.0 | 5.3 | +70.7 |
| 90% | Ventilation | 0.0 | 96.8 | −96.8 |
| 90% | Masking | 96.8 | 96.8 | 0.0 |
| 90% | Work-from-home | 96.8 | −67.9 | +164.7 |

The near-perfect fit is not a coincidence. With realistic clearance rates the infected-space fraction tracks the infected fraction, W ≈ εI/ζ, so the aerosol term δSW takes the same mathematical form as the near-field term βSI, and the two pathways produce indistinguishable epidemic curves. Accordingly, AIC-based comparison of SIR and A-SIR fits found no evidence for aerosol transmission at clearance rates of 1.4 or 14 day⁻¹: ΔAIC remained between −2.0 and −1.7 for 50% and 90% aerosol transmission, with both 25-day and 60-day windows. Only for unrealistically persistent contamination (ζ = 0.05 day⁻¹, a contamination lifetime of 20 days) did a signal appear, and then only after the epidemic peak (ΔAIC ≈ 96–120 with a 60-day window). The aerosol fraction is therefore not identifiable from epidemic curves, although it governs the effect of every intervention except masking. This reflects a general property of compartmental fitting: aggregate quantities such as R₀ are estimated robustly even when the underlying structure and individual parameters are not identifiable (Tuncer and Le, 2018; Roosa and Chowell, 2019).

## 4. Discussion

### 4.1 Implications for COVID-19 policy

The identifiability gap provides a quantitative explanation for the overoptimism of contact-reduction measures in SIR-based policy analysis. An SIR model trained on aerosol-driven data absorbs the aerosol contribution into its single transmission parameter β, which is inflated relative to the true near-field rate. Interventions that act only on near-field transmission are then credited with reducing transmission that they cannot touch. In the present analysis, the effectiveness of barriers and distancing is overestimated by a factor of about 2.4 at 50% and about 14 at 90% aerosol transmission, while the fitted model reproduces the epidemic curve perfectly.

The work-from-home result is particularly relevant. Many governments implemented stay-at-home orders early in the pandemic, projecting large reductions in transmission. The A-SIR-3 model predicts that, once aerosol transmission is substantial, concentrating occupancy in poorly ventilated homes increases household exposure enough to outweigh the reduction in contacts. Households were among the most important transmission settings during the pandemic, with pooled secondary attack rates of approximately 16.6% (Madewell et al., 2020). This is consistent with the predicted mechanism, although household data alone cannot separate the near-field and aerosol contributions, since prolonged co-residence increases exposure through both pathways.

The underestimation of ventilation is the most consequential institutional consequence. SIR-based models assign zero effectiveness to ventilation by construction, making it invisible as a policy lever. The present analysis shows that restoring proper ventilation is among the most effective individual interventions once aerosol transmission is significant, and that it suppresses the epidemic at high aerosol fractions. Investment in building ventilation infrastructure, which is durable, scalable and pathogen-agnostic, was systematically underweighted during COVID-19 (Morawska et al., 2021). Practical ventilation targets for occupied indoor spaces, expressed as air changes per hour, were proposed early in the pandemic (Allen and Ibrahim, 2021), but such targets had no entry point into SIR-based policy models. The seasonal dimension deserves emphasis: in winter, when windows are closed, the home becomes the dominant aerosol environment, and improving home ventilation is more effective than improving workplace ventilation.

### 4.2 Why surveillance data cannot settle the transmission route

A natural question is whether the aerosol fraction could be estimated from early surveillance data by comparing SIR and A-SIR fits. The present analysis indicates that it cannot, under realistic ventilation. Because contaminated indoor air clears within hours to days, the aerosol pathway acts as a near-instantaneous transmission term and produces epidemic curves that an SIR model fits exactly. The transmission route is structurally non-identifiable from case curves, yet it determines the effect of almost every intervention.

This has a direct implication for pandemic preparedness. The question of transmission route must be answered with evidence that does not come from epidemic curves: environmental measurements (indoor CO₂ as a proxy for exhaled air, air sampling for viral RNA), outbreak investigations in settings with known ventilation, and laboratory studies of aerosol emission and viral persistence. Such evidence should feed directly into the parameters of models such as A-SIR, through the emission rate ε, the clearance rate ζ and the aerosol infection rate δ. In this light, the protracted debate over airborne transmission of SARS-CoV-2 is not surprising: the epidemiological data that dominated policy modelling could not, by their nature, resolve it.

### 4.3 Limitations and future directions

Several limitations should be acknowledged. The three-space decomposition (outdoors, home, workplace) is a simplification; real urban environments contain transport vehicles, restaurants, schools, hospitals and large gathering venues, each with distinct aerosol characteristics. Each space type is also treated as a single well-mixed compartment, so that all homes mix as one space; household structure, in which near-field transmission saturates within small groups, is not represented. Extension to a richer space taxonomy and to explicit household structure is possible within the PBE framework.

The mitigation multipliers in Table 1 carry uncertainty. In particular, the masking asymmetry (δ×0.4 vs ε×0.6) and the 80% contact reduction assumed for work-from-home involve assumptions that would benefit from experimental validation. The sign of the work-from-home effect depends on the balance between contact reduction and increased home occupancy.

The fitting experiment uses synthetic prevalence data with a single, low noise level (σ = 0.2%). Real surveillance data are incidence counts subject to under-reporting, reporting delays and overdispersion; these would make the structural difference between SIR and A-SIR even harder to detect, reinforcing the conclusion of Section 4.2. The clearance-rate estimates of Annex-2 rest on a threshold interpretation of the infected-space state and on literature values for deposition and viral inactivation, which vary with temperature and humidity.

Vaccination is intentionally excluded. The interaction between vaccination and aerosol transmission is non-trivial (e.g. effects on viral load and shedding, waning of protection, variant escape) and will be addressed in a subsequent paper.

## 5. Conclusions

This study shows that the aerosol fraction of transmission can be a governing parameter that determines which public health interventions succeed or fail. Three findings have direct policy relevance.

First, the effectiveness ranking of mitigation measures inverts as the aerosol fraction increases. Barriers and distancing lose most of their effectiveness, ventilation gains, and work-from-home becomes counterproductive because it concentrates occupancy in poorly ventilated homes. These results follow from the structure of the A-SIR model and cannot be obtained from a homogeneous SIR model with a single transmission term.

Second, an SIR model fitted to aerosol-driven data reproduces the epidemic curve to within the observation noise and recovers the correct R₀, while overestimating the effectiveness of barriers and distancing by up to 85 percentage points, assigning no value to ventilation, and predicting suppression where work-from-home in fact makes the epidemic worse.

Third, under realistic ventilation the aerosol fraction is not identifiable from epidemic curves. A good fit is no evidence that a model's intervention predictions are correct.

Recommendation: epidemic modelling frameworks for respiratory pathogens should include aerosol transmission as a standard component from the first days of an outbreak; its magnitude should be constrained by environmental and physical evidence rather than inferred from case curves; and ventilation should be treated as a primary rather than supplementary intervention for airborne pathogens.

## References

Akaike, H. (1974). A new look at the statistical model identification. IEEE Trans. Autom. Control 19(6), 716–723. https://doi.org/10.1109/TAC.1974.1100705

Alexopoulos, A.H. (2025). Pandemic preparedness and new model developments for airborne diseases transmitting via aerosols. BioSystems 255, 105532.

Allen, J.G., Ibrahim, A.M. (2021). Indoor air changes and potential implications for SARS-CoV-2 transmission. JAMA 325(20), 2112–2113.

Burnham, K.P., Anderson, D.R. (2002). Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach, 2nd ed. Springer, New York.

Chu, D.K., Akl, E.A., Duda, S., et al. (2020). Physical distancing, face masks, and eye protection to prevent person-to-person transmission of SARS-CoV-2 and COVID-19: a systematic review and meta-analysis. Lancet 395(10242), 1973–1987. https://doi.org/10.1016/S0140-6736(20)31142-9

Crespo, N.C., Shifflett, S., Kosta, K., et al. (2025). Evidence of face masks and masking policies for the prevention of SARS-CoV-2 transmission and COVID-19 in real-world settings: a systematic literature review. Int. J. Environ. Res. Public Health 22(10), 1590. https://doi.org/10.3390/ijerph22101590

Holmdahl, I., Buckee, C. (2020). Wrong but useful — what Covid-19 epidemiologic models can and cannot tell us. N. Engl. J. Med. 383(4), 303–305.

Ioannidis, J.P.A., Cripps, S., Tanner, M.A. (2022). Forecasting for COVID-19 has failed. Int. J. Forecast. 38(2), 423–438.

Jit, M., et al. (2023). Reflections on epidemiological modeling to inform policy during the COVID-19 pandemic in Western Europe, 2020–23. Health Aff. 42(12), 1630–1636.

Lewis, T. (2020). How coronavirus spreads through the air: what we know so far. Scientific American. [AUTHOR: add date and URL]

Madewell, Z.J., Yang, Y., Longini, I.M., Halloran, M.E., Dean, N.E. (2020). Household transmission of SARS-CoV-2: a systematic review and meta-analysis. JAMA Netw. Open 3(12), e2031756.

Morawska, L., Milton, D.K. (2020). It is time to address airborne transmission of coronavirus disease 2019 (COVID-19). Clin. Infect. Dis. 71(9), 2311–2313.

Morawska, L., et al. (2021). A paradigm shift to combat indoor respiratory infection. Science 372(6543), 689–691.

Morawska, L., et al. (2023). Coronavirus disease 2019 and airborne transmission: science rejected, lives lost. Clin. Infect. Dis. 76(10), 1854–1859.

Noakes, C.J., Beggs, C.B., Sleigh, P.A., Kerr, K.G. (2006). Modelling the transmission of airborne infections in enclosed spaces. Epidemiol. Infect. 134(5), 1082–1091.

Roosa, K., Chowell, G. (2019). Assessing parameter identifiability in compartmental dynamic models using a computational approach: application to infectious disease transmission models. Theor. Biol. Med. Model. 16, 1.

Tuncer, N., Le, T.T. (2018). Structural and practical identifiability analysis of outbreak models. Math. Biosci. 299, 1–18.

van den Driessche, P., Watmough, J. (2002). Reproduction numbers and sub-threshold endemic equilibria for compartmental models of disease transmission. Math. Biosci. 180, 29–48.

Wang, C.C., Prather, K.A., Sznitman, J., Jimenez, J.L., Lakdawala, S.S., Tufekci, Z., Marr, L.C. (2021). Airborne transmission of respiratory viruses. Science 373(6558), eabd9149.

## Annex-1: Segregated A-SIR model

The three-space model is depicted in Fig. A1. Humans in each space i = 0, 1, 2 are tracked as S~i~, I~i~ and R~i~ (fractions of the total population), and the infected-space fractions W₁, W₂ are defined for the two indoor space types.

![Figure A1. Segregated spaces in A-SIR-3: outdoors (i = 0), homes (i = 1), workplace (i = 2).](paper_out/figA1_segregated.png){width=80%}

With N~i~ = S~i~ + I~i~ + R~i~ the occupancy of space i and X ∈ {S, I, R}, the exchange of humans between the outdoors and the indoor spaces is

$$\Phi_i^X = q_{i,in} X_0 - q_{i,out} X_i, \qquad i = 1, 2 \qquad (A1.1)$$

and the full equation set is

$$\frac{dS_0}{dt} = -\beta_0 \frac{S_0 I_0}{N_0} - \Phi_1^S - \Phi_2^S \qquad (A1.2)$$

$$\frac{dS_i}{dt} = -\beta_i \frac{S_i I_i}{N_i} - \delta_i S_i W_i + \Phi_i^S, \qquad i = 1, 2 \qquad (A1.3)$$

$$\frac{dI_0}{dt} = \beta_0 \frac{S_0 I_0}{N_0} - \gamma I_0 - \Phi_1^I - \Phi_2^I \qquad (A1.4)$$

$$\frac{dI_i}{dt} = \beta_i \frac{S_i I_i}{N_i} + \delta_i S_i W_i - \gamma I_i + \Phi_i^I, \qquad i = 1, 2 \qquad (A1.5)$$

$$\frac{dR_0}{dt} = \gamma I_0 - \Phi_1^R - \Phi_2^R, \qquad \frac{dR_i}{dt} = \gamma I_i + \Phi_i^R \qquad (A1.6)$$

$$\frac{dW_i}{dt} = \varepsilon_i (1 - W_i) I_i - \zeta_i W_i, \qquad i = 1, 2 \qquad (A1.7)$$

SIR-3 is obtained with δ~i~ = 0. With the exchange rates of Section 2.2, the steady-state occupancy of the three spaces equals the time budget, N~i~ = f~i~. Interventions that change the time budget change the exchange rates and hence the occupancy of each space, in addition to any change in β~i~, δ~i~, ε~i~ and ζ~i~.

## Annex-2: Physical estimation of the ventilation parameter ζ

### A2.1 Aerosol balance in a well-mixed indoor space

For an indoor space of volume V (m³) occupied by an infected person emitting infectious aerosol at rate E (quanta h⁻¹) (Buonanno et al., 2020), the well-mixed concentration C (quanta m⁻³) obeys

$$V\frac{dC}{dt} = E - \lambda V C, \qquad \lambda = \lambda_v + \lambda_f + k_d + k_i \qquad (A2.1)$$

where λ~v~ = Q/V is the outdoor-air change rate (ACH, h⁻¹), λ~f~ = η~f~ Q~f~/V the equivalent rate of filtration (portable HEPA cleaner or recirculation filter of efficiency η~f~), k~d~ the gravitational and turbulent deposition rate (≈ 0.1–0.5 h⁻¹ for 1–5 µm particles; Thatcher et al., 2002) and k~i~ the biological inactivation rate (k~i~ = ln2/t½ ≈ 0.63 h⁻¹ for t½ ≈ 1.1 h; van Doremalen et al., 2020). The steady state is C~ss~ = E/(λV); when the source leaves, C(t) = C~ss~ exp(−λt).

### A2.2 Link to the W compartment

A space is counted as infected (W) while C exceeds an infectious threshold C\*, the concentration at which the Wells–Riley probability of infection, P = 1 − exp(−pCt~exp~), for a typical occupant (breathing rate p, exposure time t~exp~) becomes non-negligible (Riley et al., 1978; Rudnick and Milton, 2003). The clearance time of a contaminated space is t~c~ = ln(C~ss~/C\*)/λ, so that

$$\zeta = \frac{1}{t_c} = \frac{\lambda}{\ln\left(\dfrac{E}{\lambda V C^*}\right)} \qquad (A2.2)$$

Defining the baseline contamination intensity r₀ = C~ss~/C\* at the baseline removal rate λ~b~, an intervention that raises the total removal rate to λ = xλ~b~ gives

$$\frac{\zeta}{\zeta_b} = \frac{x \ln r_0}{\ln(r_0 / x)}, \qquad 1 \le x < r_0 \qquad (A2.3)$$

For x ≥ r₀ the steady concentration never reaches C\* and the space cannot become infected. Eq. (A2.3) shows that ζ grows faster than λ: ventilation both shortens clearance and lowers the contamination level from which clearance starts.

### A2.3 Clearance rates for proper and poor ventilation

With r₀ = 10, k~d~ = 0.2 h⁻¹ and k~i~ = 0.63 h⁻¹ (no baseline filtration):

- Home, 0.5 ACH: λ = 1.33 h⁻¹, t~c~ ≈ 1.7 h, ζ₁ ≈ 14 day⁻¹.
- Workplace, 2 ACH: λ = 2.83 h⁻¹, t~c~ ≈ 0.8 h, ζ₂ ≈ 29 day⁻¹.

These are the proper-ventilation values of Section 2.3. Even with no ventilation at all, deposition and inactivation impose a lower bound on λ; with the values above, and allowing for the higher contamination intensity r₀ of an unventilated room, this gives ζ of about 7 day⁻¹. The poor (winter) baseline used in this study, ζ₁ = 1.4 and ζ₂ = 2.9 day⁻¹, ten times below the proper values, therefore implies both closed windows and slower viral inactivation, as occurs in cold, dry indoor air. Spaces that are never ventilated, such as classrooms with permanently closed windows, may approach this regime.

### A2.4 Ventilation required for the multipliers used in this work

**Table A2.** Total removal-rate ratio x = λ/λ~b~ required to obtain the ζ multipliers of Table 1 and Fig. 6.

| r₀ = C~ss~/C\* | ζ × 5 | ζ × 10 |
|---|---|---|
| 5 | 2.35 | 3.06 |
| 10 | 2.78 | 3.99 |
| 30 | 3.26 | 5.17 |
| 100 | 3.61 | 6.08 |

For reference, 4–6 ACH has been proposed as a target for small indoor spaces (Allen and Ibrahim, 2021), and airborne-infection isolation rooms are designed for at least 12 ACH. Starting from poor winter ventilation, the multipliers of Table 1 and Fig. 6 therefore correspond to restoring ventilation to ordinary good practice.

### A2.5 Remark: ventilation also acts on ε

Because C~ss~ ∝ 1/λ, increased ventilation also lowers the probability that an occupied clean space reaches C\* in the first place, i.e. it reduces the effective contamination rate ε. Table 1 conservatively assigns ε × 1.0 to ventilation, so the ventilation benefits reported in Section 3 are, within the model, lower bounds.

### References for Annex

Buonanno, G., Stabile, L., Morawska, L. (2020). Estimation of airborne viral emission: quanta emission rate of SARS-CoV-2 for infection risk assessment. Environ. Int. 141, 105794.

Riley, E.C., Murphy, G., Riley, R.L. (1978). Airborne spread of measles in a suburban elementary school. Am. J. Epidemiol. 107(5), 421–432.

Rudnick, S.N., Milton, D.K. (2003). Risk of indoor airborne infection transmission estimated from carbon dioxide concentration. Indoor Air 13(3), 237–245.

Thatcher, T.L., Lai, A.C.K., Moreno-Jackson, R., Sextro, R.G., Nazaroff, W.W. (2002). Effects of room furnishings and air speed on particle deposition rates indoors. Atmos. Environ. 36(11), 1811–1819.

van Doremalen, N., et al. (2020). Aerosol and surface stability of SARS-CoV-2 as compared with SARS-CoV-1. N. Engl. J. Med. 382(16), 1564–1567.
