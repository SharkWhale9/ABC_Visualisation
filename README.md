# ABC / ABC-MCMC SIR Dashboard

## Objective

This dashboard demonstrates how **Approximate Bayesian Computation (ABC)** can infer parameters of a stochastic simulator by replacing likelihood evaluation with repeated simulation and comparison. It uses one fixed synthetic SIR epidemic to compare two algorithms—independent ABC rejection sampling and ABC-MCMC—and shows how the resulting inference depends on the summary statistic, tolerance, simulation budget, and MCMC proposal scale.

The main statistical message is not that ABC-MCMC should produce a narrower posterior than rejection sampling. At matched settings, both methods target the same tolerance-defined ABC posterior. The dashboard instead shows that:

- rejection sampling is simple and produces independent accepted draws, but can waste most simulations across the prior;
- ABC-MCMC proposes locally and can produce compatible moves more efficiently, but creates dependent samples that require burn-in, proposal tuning, and mixing checks;
- the summary statistic determines what information ABC can recover;
- `beta` and `gamma` can remain individually uncertain along a ridge while the derived ratio `R0 = beta / gamma` is more strongly constrained;
- small differences between finite-run credible intervals should not automatically be interpreted as differences between the algorithms.

## Suggested walkthrough

1. Start with the fixed observed epidemic and note its peak, peak timing, and final size.
2. In Panel 1, move `epsilon` and watch the accepted rejection-sampling region change without new simulations.
3. Switch between the rich and poor summaries to see how discarded information changes the inferred parameter region.
4. In Panel 2, inspect the ABC-MCMC path and trace plots, then change the proposal scales and rerun the chain.
5. In Panel 3, compare simulation efficiency, joint parameter geometry, and the credible intervals for `R0`.
6. Finish with the validation section, which checks that the stochastic simulator behaves sensibly relative to a deterministic SIR reference.

---

## 1. What this is

This is a self-contained browser dashboard for likelihood-free parameter inference in a stochastic SIR epidemic simulator. It compares independent ABC rejection sampling with ABC-MCMC and lets a visitor see how tolerance, summary statistics, and random-walk proposal scales affect the approximate posterior.

The implementation is intended as an interactive statistics portfolio piece. The stochastic simulator, binomial sampler, deterministic reference solver, rejection sampler, ABC-MCMC transition, validation checks, and comparison metrics are implemented explicitly in vanilla JavaScript rather than delegated to an ABC or MCMC package.

---

## 2. The problem being solved

### The SIR model

The model divides a fixed population into three compartments:

- `S`: susceptible people;
- `I`: currently infected people;
- `R`: recovered people.

The unknown parameters are:

- `beta`: transmission rate. Larger values generally make infections grow faster.
- `gamma`: recovery rate. Larger values remove people from `I` more quickly and shorten the infectious period.

The dashboard also studies the derived quantity:

```text
R0 = beta / gamma
```

`R0` summarizes the balance between transmission and recovery. Different values of `beta` and `gamma` can produce a similar ratio and therefore similar broad epidemic dynamics.

### Why ABC is used

The dashboard performs inference by simulation:

1. propose `beta` and `gamma`;
2. simulate a stochastic epidemic;
3. reduce the simulated trajectory to selected summary statistics;
4. compare those summaries with a fixed observed epidemic;
5. retain parameter values whose distance is below a tolerance `epsilon`.

This is the standard motivation for ABC: a model can generate data, but its likelihood is unavailable, impractical, or intentionally treated as a black box.

### Likelihood caveat for this exact teaching model

ABC is **not mathematically necessary in the strongest sense for this particular simplified SIR model**. Conditional on a fully observed daily `S`, `I`, and `R` path, the implemented transitions are explicit binomial distributions and a complete-path likelihood could be written down.

The dashboard should therefore be understood as a controlled demonstration of:

- likelihood-free workflow;
- summary-statistic information loss;
- tolerance effects;
- rejection versus MCMC simulation efficiency;
- parameter non-identifiability;
- Monte Carlo variability.

It is not proof that this toy SIR model has no tractable likelihood under any observation scheme.

### What “observed data” means here

The observed data are synthetic, not real epidemiological measurements. They are one stochastic epidemic generated from known parameters and then embedded permanently in the HTML file.

The inferential exercise is:

> Pretend `beta` and `gamma` are unknown, then ask which parameter values generate epidemics that resemble this one fixed synthetic outbreak.

---

## 3. The model and ground truth

### Stochastic daily updates

At each day, both random counts are drawn from the **start-of-day** compartments:

```text
infection_probability = 1 - exp(-beta * I / N * dt)
recovery_probability  = 1 - exp(-gamma * dt)

new_infections ~ Binomial(S, infection_probability)
new_recoveries ~ Binomial(I, recovery_probability)

S_next = S - new_infections
I_next = I + new_infections - new_recoveries
R_next = R + new_recoveries
```

Infections and recoveries are conditionally independent given the current state. A person infected during the current daily step cannot recover until a later step because recoveries are drawn from the original `I`, not from `I + new_infections`.

### Fixed simulation settings

| Setting | Value |
|---|---:|
| Population `N` | 1,000 |
| Initial susceptible `S0` | 995 |
| Initial infected `I0` | 5 |
| Initial recovered `R0` | 0 |
| Horizon `T` | 100 days |
| Time step `dt` | 1 day |

### Priors

The priors are independent uniforms:

```text
beta  ~ Uniform(0.05, 1.00)
gamma ~ Uniform(0.02, 0.50)
```

The feasible support of the derived ratio is therefore:

```text
minimum R0 = 0.05 / 0.50 = 0.10
maximum R0 = 1.00 / 0.02 = 50.00
```

This is only the feasible support. Independent uniform priors on `beta` and `gamma` do **not** induce a uniform prior distribution on their ratio.

### Ground truth

```text
beta_true  = 0.4
gamma_true = 0.1
R0_true    = 4.0
```

### Fixed observed epidemic

```text
observed-generation seed = 4041001
```

The observed epidemic is one literal, hardcoded stochastic realization. It was generated once using the same SIR simulator and seed `4041001`, after which its `S`, `I`, `R`, `newInfections`, and `newRecoveries` arrays were embedded in the page.

It is **not regenerated**:

- on page load;
- when validation is rerun;
- when a new inference seed is requested;
- when the summary statistic or tolerance changes.

Its displayed summaries are:

```text
peak infected = 410
peak day      = 24
final size    = 980
```

Final epidemic size is:

```text
N - S(100)
```

It counts everyone infected at least once by day 100. If infections remained after day 100, this would not necessarily equal the eventual post-extinction epidemic size.

### Custom binomial sampler

The page does not use a probability library for binomial draws.

For reflected `p <= 0.5`, it uses:

- inverse-CDF recursion when `n * p < 30`;
- tangent-envelope rejection sampling when `n * p >= 30`;
- the identity `Binomial(n,p) = n - Binomial(n,1-p)` when the original `p > 0.5`.

The large-mean branch uses a Lanczos approximation to `logGamma` when evaluating binomial probability masses. This avoids one Bernoulli draw per individual, but the dashboard’s epidemic-level validation is not a formal exhaustive distributional test of the sampler for every possible `(n,p)`.

---

## 4. Summary statistics, distance, and tolerance

ABC does not compare the full epidemic trajectory directly. It compresses each simulation into a summary vector.

### Rich summary

```text
[peak infection height, peak timing, final epidemic size]
= [max(I), index of first maximum of I, N - S(T)]
```

For the fixed observed epidemic:

```text
[410, 24, 980]
```

If the infected trajectory has a tied maximum, the earliest maximum is used.

### Poor summary

```text
[final epidemic size]
= [N - S(T)]
```

For the fixed observed epidemic:

```text
[980]
```

The poor summary deliberately discards peak height and peak timing. It therefore cannot distinguish a fast, sharp epidemic from a slower, longer epidemic if both infect similar totals.

### Distance metric

The implementation uses normalized Euclidean distance.

```text
peak_difference   = (simulated_peak - 410) / N
timing_difference = (simulated_peak_day - 24) / T
final_difference  = (simulated_final_size - 980) / N
```

For the rich summary:

```text
rich_distance = sqrt(
  peak_difference^2
  + timing_difference^2
  + final_difference^2
)
```

For the poor summary:

```text
poor_distance = abs(final_difference)
```

With `N = 1000` and `T = 100`, the count-valued summaries are divided by 1,000 and peak timing is divided by 100 days. This prevents a statistic from dominating only because it is measured on a numerically larger scale.

This is fixed range-based normalization. It is not standardization by simulated, prior, or posterior standard deviations, so it does not guarantee that the three rich-summary components have equal inferential importance.

### Tolerance `epsilon`

A parameter proposal is ABC-compatible when:

```text
distance < epsilon
```

Default:

```text
epsilon = 0.080
```

Slider range:

```text
0.005 to 0.250
step = 0.005
```

A smaller `epsilon` demands a closer match. It usually produces a more concentrated ABC approximation but lowers acceptance and increases Monte Carlo instability.

A larger `epsilon` accepts weaker matches. It increases acceptance but makes the approximation less informed by the observed epidemic. At sufficiently large values, the ABC result approaches the prior.

The tolerance therefore controls a statistical-accuracy versus computational-cost trade-off.

---

## 5. Panel-by-panel explanation

### Observed epidemic panel

#### What it shows

- The fixed `S(t)`, `I(t)`, and `R(t)` trajectories for days 0 through 100.
- The infection peak at `I = 410` on day 24.
- The rich and poor observed summary vectors.
- Population, time horizon, true `R0`, and observed-data provenance seed.

#### What it is trying to tell you

The fixed outbreak is the common target for every later calculation. The panel also makes visible how aggressively the summaries compress the original data: the rich summary reduces an entire trajectory to three numbers, while the poor summary keeps only one.

#### What not to over-interpret

This is one stochastic realization. A different epidemic generated at the same true parameters would generally have a different peak, timing, and final size.

---

### Panel 1: ABC rejection sampling

#### Algorithm actually run

1. Draw 4,000 independent `(beta, gamma)` pairs from the priors.
2. Simulate one stochastic SIR epidemic for each pair.
3. Cache each proposal’s:
   - `beta`;
   - `gamma`;
   - full infected trajectory `I(t)`;
   - peak height;
   - peak day;
   - final size.
4. Compute the selected normalized summary distance.
5. Accept when `distance < epsilon`.

#### What is cached versus recomputed

The 4,000 epidemics are generated once per inference seed.

- Moving the `epsilon` slider only re-filters cached distances.
- Changing rich versus poor summary recomputes distances from cached raw summaries.
- Neither interaction reruns the SIR simulator.
- The trajectory panel displays the eight accepted simulations with the smallest current distances.

#### Prior-proposal scatter plot

Every point is one proposed `(beta, gamma)` pair.

The plot distinguishes:

- rejected proposals;
- accepted proposals;
- the true point `(0.4, 0.1)`.

The accepted cloud approximates the ABC posterior under the current summary and tolerance.

A compact circular cloud would suggest that `beta` and `gamma` are separately constrained. A diagonal or curved ridge indicates a trade-off: different transmission and recovery rates can produce similar epidemic summaries.

#### Accepted-trajectory plot

The closest accepted infected curves are overlaid on the observed infected trajectory.

This answers a useful question that the parameter scatter alone cannot:

> Do simulations that pass the summary-statistic rule also resemble the observed epidemic as full trajectories?

Under the poor summary, accepted curves may have similar final totals while differing substantially in peak height and timing.

#### Acceptance rate

```text
acceptance rate = accepted proposals / 4000
```

This measures how much of the independent simulation budget passed the current ABC rule. It is not the probability that the model or a parameter value is correct.

#### What to try

- Lower `epsilon` and watch the accepted cloud shrink while the acceptance rate falls.
- Increase `epsilon` and observe the result move toward the prior.
- Switch to the poor summary and examine how the accepted parameter region and epidemic curves broaden.

#### What Panel 1 is trying to teach

Rejection sampling is transparent and its accepted draws are independent:

```text
propose globally -> simulate -> compare -> keep or discard
```

Its weakness is computational. It continues proposing across the full prior even when only a small region produces acceptable epidemics.

#### What not to over-interpret

The rejection result uses only 4,000 proposals. At strict tolerances, the number accepted may be small, making empirical means, quantiles, modes, and credible-interval endpoints unstable.

For the default rich-summary run examined during development, only 120 proposals were accepted. The resulting rejection `R0` credible interval was unusually narrow relative to repeated batches at the same settings. It should not be interpreted as evidence that rejection sampling targets a genuinely narrower posterior than ABC-MCMC.

---

### Panel 2: ABC-MCMC

#### Why MCMC is introduced

ABC-MCMC tries to avoid repeatedly proposing parameters in regions that are clearly incompatible with the observed summaries. Once it finds a compatible state, it proposes nearby values rather than drawing globally from the prior.

#### Initialization

The algorithm repeatedly:

1. draws a parameter pair from the priors;
2. simulates an epidemic;
3. computes its distance;
4. stops at the first draw with `distance < epsilon`.

Initialization fails after 25,000 attempts and asks the visitor to increase `epsilon` if no compatible state is found.

#### Chain configuration

```text
iterations = 6500
burn-in    = 500
retained   = 6000 states
```

Default proposal scales:

```text
sigma_beta  = 0.045
sigma_gamma = 0.018
```

Slider ranges:

```text
sigma_beta:  0.005 to 0.150, step 0.005
sigma_gamma: 0.002 to 0.080, step 0.002
```

#### Proposal distribution

At each transition:

```text
beta_proposed  = beta_current  + sigma_beta  * Normal(0,1)
gamma_proposed = gamma_current + sigma_gamma * Normal(0,1)
```

#### Acceptance rule

A proposal outside either prior range is rejected immediately and does not call the SIR simulator.

An in-bounds proposal receives one fresh epidemic simulation and is accepted exactly when:

```text
distance(proposed simulation, observed summaries) < epsilon
```

If rejected, the next chain state is an exact copy of the current state.

#### Relation to the textbook ABC-MCMC rule

The implementation uses a hard-threshold Marjoram-style ABC-MCMC rule with the Metropolis-Hastings ratio simplified to one for valid proposals.

That simplification is valid here because:

- both priors are uniform inside fixed bounds;
- the untruncated Gaussian random-walk proposal is symmetric;
- the code rejects proposals outside the prior support before simulation;
- within the support, the prior and proposal ratios cancel.

The algorithm does **not** accept a proposal merely because its distance is better than the current state’s distance. It accepts any in-bounds proposal whose newly simulated data pass the fixed threshold.

#### Chain path plot

The path plot contains:

- the matched rejection-compatible region in the background;
- the burn-in path;
- the post-burn-in path;
- the initial state;
- the final/current state;
- the true parameter.

It is intended to show whether the chain explores the compatible region or remains trapped in one part of it.

#### Trace plots

The beta and gamma traces show parameter values against iteration number.

Look for:

- repeated movement across the explored range;
- slow one-directional drift;
- long flat stretches caused by rejection;
- excessive stickiness;
- failure to revisit regions;
- differences between burn-in and post-burn-in behavior.

The red references show the true generating values. The first 500 states are marked as burn-in.

#### Proposal-scale interpretation

If the proposal scales are too small:

- acceptance may be high;
- moves are local;
- autocorrelation can be severe;
- the chain may explore the ridge slowly.

If the proposal scales are too large:

- more proposals leave the prior;
- more simulated proposals fail the tolerance;
- acceptance falls;
- the chain repeats its current state more often.

#### Settings snapshot and stale results

When the chain begins, it snapshots:

```text
summary mode
epsilon
sigma_beta
sigma_gamma
```

Changing the controls later does not retroactively alter the completed chain. Instead, the chain and comparison are marked stale until the chain is rerun.

This differs from rejection sampling because an `epsilon` or summary change would have altered every historical MCMC transition decision.

#### What Panel 2 is trying to teach

ABC-MCMC can concentrate simulation effort in the compatible region, but this efficiency creates dependent output.

The method exchanges rejection sampling’s global inefficiency for:

- burn-in;
- proposal tuning;
- autocorrelation;
- repeated states;
- the need to assess mixing.

#### What not to over-interpret

- Six thousand retained states are not six thousand independent samples.
- A high acceptance rate is not automatically good; it may indicate steps that are too small.
- A low acceptance rate can result from oversized steps, poor simulated matches, or out-of-bounds proposals.
- The dashboard does not compute effective sample size, autocorrelation functions, R-hat, or multi-chain convergence diagnostics.
- Burn-in is fixed at 500 rather than selected from a formal diagnostic.

---

### Panel 3: Matched comparison

#### What “matched settings” means

The comparison uses the completed chain’s frozen:

- summary mode;
- `epsilon`.

It then re-filters the cached rejection proposals using those exact settings. Both methods are therefore compared against the same tolerance-defined ABC target.

The MCMC proposal scales affect how the chain explores that target. They do not change which cached rejection proposals are compatible.

If the visible controls differ from the completed chain settings, the comparison is marked stale and continues to describe the last valid run.

#### Central statistical point

At matched settings, rejection and ABC-MCMC are intended to approximate the **same ABC posterior**.

The comparison is therefore not asking:

> Which method creates the narrower posterior?

It is asking:

> How does each method spend its simulation budget, and what computational limitations accompany its output?

#### Joint parameter plot

The large joint plot overlays:

- accepted rejection points;
- post-burn-in MCMC states;
- the true `(beta, gamma)` point;
- each method’s sample mean.

Separate marginal histograms show `beta` and `gamma` individually.

The plot displays at most 1,800 evenly spaced points from each method for browser performance, but numerical summaries use all accepted rejection draws and all retained MCMC states.

#### Rejection yield

```text
1000 * accepted rejection proposals / 4000
```

This reports how many compatible independent draws rejection produced per 1,000 SIR simulations.

#### MCMC move yield

```text
1000 * accepted MCMC moves
/ (initialization simulations + in-bounds transition simulations)
```

This reports how often a simulator call generated an accepted random-walk move after accounting for initialization.

Out-of-bounds proposals are excluded because they do not call the simulator.

Move yield is a computational measure, not an effective-sample-size measure. Accepted moves and retained states remain dependent.

#### Simulator budget

For rejection:

```text
4000 simulator calls
```

For ABC-MCMC:

```text
initialization attempts
+ in-bounds transition simulations
```

#### Central 80% footprint

The footprint is:

```text
(beta 90th percentile - beta 10th percentile)
*
(gamma 90th percentile - gamma 10th percentile)
```

The result is divided by the area of the full prior rectangle.

This is an axis-aligned descriptive rectangle, not:

- a joint highest-posterior-density region;
- a formal credible region;
- a geometry-aware measure of a diagonal or curved ridge.

A narrow diagonal ridge can still generate a relatively large axis-aligned footprint.

#### What Panel 3 is trying to tell you

- Rejection continues drawing across the entire prior.
- ABC-MCMC proposes locally once it reaches a compatible region.
- ABC-MCMC can therefore produce more compatible moves per simulator call.
- That efficiency does not make its retained states independent.
- A correctly mixing chain should reproduce the same target as rejection, not an artificially narrower one.
- Summary-statistic choice can have a larger effect on the result than algorithm choice.

---

### R0 credible intervals

For every accepted or retained sample, the dashboard computes:

```text
R0 = beta / gamma
```

It calculates separately for rejection and ABC-MCMC:

- mean `R0`;
- median `R0`;
- 5th percentile;
- 95th percentile;
- equal-tailed 90% credible interval;
- credible-interval width;
- width as a fraction of the feasible `R0` support.

No `R0` values are clipped. Valid samples satisfy `gamma >= 0.02`, so division by zero is not expected. A non-finite value or denominator below the prior lower bound is treated as malformed and causes the summary to be omitted rather than silently modified.

#### Why R0 can be better identified than beta and gamma

Consider:

```text
beta = 0.4, gamma = 0.10 -> R0 = 4
beta = 0.6, gamma = 0.15 -> R0 = 4
```

The individual rate values differ, but their ratio is the same.

If epidemic summaries mainly constrain the balance between transmission and recovery, accepted points can form a ridge of approximately constant `beta / gamma`. The rates remain individually dispersed while the ratio is more tightly constrained.

This is the dashboard’s most important inferential finding:

> A scientifically meaningful combination of parameters can be better identified than the individual parameters that define it.

#### Default rich-summary run

The default run examined during development produced approximately:

```text
True R0 = 4.00

Rejection median R0 = 3.71
Rejection 90% credible interval = [3.09, 4.50]

ABC-MCMC median R0 = 3.67
ABC-MCMC 90% credible interval = [3.06, 4.67]
```

Both intervals contain the true value and both are narrow relative to the feasible support `[0.10, 50.00]`.

The rejection interval happens to be narrower in this run. That ordering should **not** be interpreted as a method effect. The rejection interval is based on only 120 accepted points and was found to be unusually narrow relative to repeated independent rejection batches at the same settings. The MCMC interval was not unusually wide.

The defensible conclusion is:

- both methods identify a constrained `R0` region under the rich summary;
- `beta` and `gamma` remain more dispersed individually;
- small differences between the two finite-run interval widths are Monte Carlo variation, not evidence that one method targets a different posterior.

#### Poor-summary behavior

In the tested poor-summary run:

```text
Rejection 90% credible interval = [2.26, 18.14]
ABC-MCMC 90% credible interval  = [2.50, 20.66]
```

These much wider intervals demonstrate that final epidemic size alone does not reliably identify the transmission-recovery balance.

The lesson is:

> ABC cannot recover information that the chosen summary statistic has discarded.

A more sophisticated sampling algorithm cannot repair an uninformative summary.

#### Credible-interval terminology

These are **credible intervals**, not confidence intervals. They summarize the empirical ABC posterior samples produced under the model, prior, summaries, tolerance, and finite simulation run.

They are equal-tailed 5th-to-95th percentile intervals, not highest-posterior-density intervals.

---

### Summary-statistic callout

The summary callout interprets the completed chain’s actual settings.

Under the rich summary:

- peak height constrains epidemic scale;
- peak timing constrains epidemic speed;
- final size constrains total spread;
- the joint information narrows the compatible dynamic behavior.

Under the poor summary:

- only total spread remains;
- fast/sharp and slow/long outbreaks can look equivalent;
- the beta-gamma region widens;
- the `R0` credible intervals widen;
- the approximate posterior center can move away from the generating point.

This section is intended to show that the definition of “similar data” is part of the model-based inference, not a harmless plotting choice.

---

### Simulator validation panels

The validation section is not another inference method. It tests whether the stochastic SIR simulator behaves sensibly.

The page compares 500 stochastic trajectories with a deterministic SIR ODE reference in each regime:

| Regime | beta | gamma | R0 | Validation seed |
|---|---:|---:|---:|---:|
| Below threshold | 0.08 | 0.20 | 0.4 | 20260806 |
| Moderate outbreak | 0.25 | 0.10 | 2.5 | 20261806 |
| Ground-truth regime | 0.40 | 0.10 | 4.0 | 20262806 |

The deterministic reference is solved using fourth-order Runge-Kutta with an internal step of `0.05` days and returned at daily intervals.

The validation checks:

- population conservation;
- non-negative integer stochastic states;
- non-increasing susceptible counts;
- non-decreasing recovered counts;
- normalized RMSE of the stochastic ensemble mean against the deterministic infected trajectory;
- peak height;
- peak timing;
- final epidemic size.

The plots show:

- deterministic ODE trajectory;
- stochastic ensemble mean;
- pointwise 10th-to-90th percentile band;
- one stochastic realization.

#### What validation establishes

It shows that the stochastic simulator has valid compartment behavior and reproduces the broad qualitative scale and timing of deterministic SIR dynamics across several regimes.

#### What validation does not establish

It does not prove:

- epidemiological realism;
- exhaustive correctness of every binomial-sampler case;
- accuracy of the ABC posterior;
- exact equality between the daily stochastic model and the continuous ODE.

The stochastic process and deterministic ODE are not algebraically identical, so modest differences are expected.

---

### Expandable source-code panel

The implementation section displays the exact JavaScript functions and configuration objects used by the page, including:

- stochastic SIR simulation;
- custom binomial sampling;
- deterministic RK4 reference;
- observed-data arrays;
- summary statistics and distance;
- rejection sampling and caching;
- ABC-MCMC initialization and transition rule;
- comparison and `R0` summaries.

It is included for transparency and to show that the inference methods are implemented explicitly.

The displayed source is documentation generated from the functions already defined by the page. Editing only the displayed text in the browser would not change the running algorithms.

---

## 6. Implementation notes

### Client-side architecture

- Everything runs in the browser.
- There is no backend.
- There are no external data files.
- The one external runtime dependency is Plotly.js `2.35.2`, loaded from a CDN.
- The SIR simulator, random-number generator, binomial sampler, RK4 solver, rejection sampler, and ABC-MCMC chain are written in vanilla JavaScript.

Because Plotly is loaded remotely, the HTML needs network access for chart rendering unless Plotly is vendored locally. The numerical JavaScript is embedded in the page, but the charts cannot render if the CDN request fails.

### Caching and performance

- Rejection sampling stores 4,000 parameter pairs.
- It stores raw epidemic summaries and all 4,000 infected trajectories of length 101.
- The rejection batch runs in chunks of 200 and yields to the browser between chunks.
- `epsilon` changes re-filter cached distances.
- Summary changes recompute distances without rerunning the simulator.
- ABC-MCMC must be rerun after summary, `epsilon`, or proposal-scale changes because those settings alter the transition history.
- MCMC displays are refreshed every 250 transitions rather than after every transition.
- The comparison panel performs no new epidemic simulations; it is derived from the rejection and MCMC caches.
- `R0` values and credible intervals are calculated directly from cached `beta` and `gamma` samples.

### Reproducibility

- Simulation randomness uses an explicit `mulberry32` seeded pseudo-random number generator after a seed is selected.
- Default inference base seed: `314159265`.
- Session storage key: `sessionStorage["sir-abc-rejection-seed-v1"]`.
- Reloads in the same browser tab/session reproduce the same rejection batch while the stored seed remains available.
- The MCMC seed is derived as `base_seed XOR 0x9E3779B9`.
- With the default base seed, the derived MCMC seed is `2358167832`.
- “New inference seed” uses `crypto.getRandomValues` when available, with a time/`Math.random` fallback.
- A new inference seed reruns rejection sampling and ABC-MCMC.
- The fixed observed epidemic never changes with the inference seed.
- Validation uses separate fixed seeds and therefore reproduces the same validation ensembles.

### Top-level buttons

#### New inference seed

Creates and stores a new base seed, then reruns:

- the 4,000 rejection proposals;
- the 6,500-transition ABC-MCMC chain.

It does not modify the fixed observed epidemic.

This button demonstrates finite-run Monte Carlo variability: accepted counts, sample clouds, means, and interval estimates can change across seeds.

#### Rerun validation

Reruns the three deterministic-versus-stochastic validation experiments using their fixed seeds.

It does not change:

- the observed epidemic;
- the inference seed;
- the current rejection or MCMC target.

### Status indicators

Status boxes indicate whether results are:

- running;
- complete;
- stale;
- unavailable because of an error.

A stale result remains valid for the settings under which it was generated. It simply no longer corresponds to the currently visible controls.

### Repository files

For a minimal repository or personal-site deployment, keep:

- `abc-sir-dashboard-r0-credible-intervals.html`: the complete dashboard. It may be renamed to `index.html`.
- `README.md`: this documentation.

Optional development artifacts:

- `abc-sir-dashboard-r0-credible-intervals.patch`: unified diff of the `R0` addition against the previous dashboard version.
- `abc-sir-dashboard-r0-test-results.json`: captured browser-test output for the `R0` additions. The dashboard does not load or require this file.
- `abc-sir-dashboard-r0-interval-investigation.json`: detailed diagnostics from the finite-run interval investigation. It is not used by the application.
- screenshots: optional documentation or preview assets.

### Known limitations and simplifications

- **ABC is pedagogical for this exact toy model.** A complete-path binomial transition likelihood could be written for a fully observed SIR trajectory.
- **The observed epidemic is synthetic.** No claim is made about a real disease or population.
- **Finite budgets matter.** Rejection uses 4,000 proposals and MCMC uses one 6,500-transition chain. Apparent means, modes, widths, and tail quantiles can vary with the seed.
- **Rejection tail estimates can be unstable.** With few accepted points, the 5th and 95th percentiles are determined by only a handful of observations.
- **MCMC states are dependent.** The dashboard reports retained-state count but does not calculate effective sample size.
- **Extreme epsilon values have predictable failure modes.** Very small values can yield no rejection points or fail MCMC initialization after 25,000 attempts. Very large values accept weak matches and approach the prior.
- **Summaries lose information.** Even the rich vector compresses a 101-day trajectory to three numbers.
- **No observation-error model is included.** There is no reporting noise, under-ascertainment, delay, missingness, or measurement process.
- **One simulated epidemic is used per proposal.** The code does not average repeated simulations at the same parameter value and does not use a synthetic likelihood.
- **The chain is not adaptively tuned.** Proposal scales are fixed during each run.
- **Burn-in is fixed.** The first 500 states are always discarded without a formal convergence criterion.
- **No formal MCMC convergence diagnostics are computed.** There is one chain and no R-hat, effective sample size, integrated autocorrelation time, or multi-chain comparison.
- **The binomial sampler is custom.** Its two numerical branches are practical implementations, not a call to a separately validated probability library.
- **The deterministic comparison is only a reference.** It uses a continuous ODE while the stochastic simulator uses one-day discrete transitions.
- **Peak timing records the first maximum.** Tied peak plateaus are represented by their earliest day.
- **Final size is measured at day 100.** It is not necessarily the eventual epidemic size if infections remain.
- **The central footprint ignores ridge geometry.** It is based on marginal quantiles and an axis-aligned rectangle.
- **Credible intervals inherit Monte Carlo error.** Rejection intervals use a finite independent accepted sample; MCMC intervals use dependent retained states.
- **Credible intervals are equal-tailed.** They are not highest-posterior-density intervals.
- **`R0` support normalization is descriptive.** The ratio prior is not uniform on `[0.10, 50.00]`.
- **`R0` values are not clipped.** Valid samples are protected by the `gamma >= 0.02` prior bound; malformed denominators cause the summary to be omitted.
- **A narrower displayed interval is not automatically better.** It can reflect Monte Carlo noise, autocorrelation, or incomplete exploration rather than genuine posterior concentration.

---

## What the dashboard establishes

The dashboard supports the following conclusions:

1. **ABC replaces likelihood evaluation with simulation and comparison.**
2. **Tolerance controls a fit-versus-acceptance trade-off.**
3. **Summary-statistic choice controls what can be identified.**
4. **Independent rejection sampling can waste simulations across the prior.**
5. **ABC-MCMC can concentrate simulator calls in the compatible region.**
6. **MCMC efficiency comes with dependence, tuning, and mixing concerns.**
7. **At matched settings, both methods target the same ABC posterior.**
8. **Beta and gamma can trade off along a ridge.**
9. **The derived ratio `R0` can be better identified than either rate separately.**
10. **Finite-run credible-interval differences should be interpreted cautiously.**

## What the dashboard does not establish

It does not show that:

- ABC-MCMC inherently produces a narrower posterior;
- a higher MCMC acceptance rate guarantees better mixing;
- 6,000 retained states are equivalent to 6,000 independent draws;
- the true parameter must be the sample mean or posterior center for one stochastic dataset;
- the rich summary preserves all information in the epidemic trajectory;
- rejection’s current narrower `R0` interval represents a real algorithmic advantage;
- the model is a realistic description of a specific epidemic;
- ABC is the only possible inferential method for this toy simulator.

---

## One-paragraph interpretation

The page creates one fixed stochastic epidemic and asks which transmission and recovery rates can reproduce its important features. ABC rejection sampling searches by drawing independently across the entire prior, while ABC-MCMC searches locally after reaching a compatible region. MCMC can use simulator calls more efficiently, but its samples are dependent and require proposal tuning and mixing assessment. More importantly, the dashboard shows that inference is controlled by the selected summaries: peak height, peak timing, and final size constrain epidemic dynamics much more strongly than final size alone. The compatible `beta` and `gamma` values can form a ridge, meaning the two rates are not separately well identified, while their ratio `R0 = beta / gamma` remains comparatively constrained. Both algorithms should approximate the same ABC posterior at matched settings, so small differences between their displayed interval widths are finite-run Monte Carlo effects rather than the intended scientific conclusion.

---

## 7. How to extend this

- Run multiple ABC-MCMC chains and add R-hat, effective sample size, autocorrelation, and integrated autocorrelation-time diagnostics.
- Add adaptive or pilot-tuned proposal covariance, including correlated beta-gamma proposals aligned with the posterior ridge.
- Compare hard-threshold ABC with kernel-weighted ABC, ABC-SMC, or sequential Monte Carlo.
- Add alternative summaries such as early growth rate, epidemic duration, selected time-point counts, or learned summaries.
- Add a realistic observation model with reporting noise, delay, under-ascertainment, or missing data.
- Compare the ABC approximation with an exact transition likelihood or particle-filter method for this small SIR model.
- Infer additional parameters or initial conditions, or replace SIR with SEIR or another simulator.
- Add posterior predictive checks based on fresh simulations from retained parameter draws.
- Add Monte Carlo uncertainty estimates for credible-interval endpoints.
- Add a dedicated `R0` marginal plot or geometry-aware joint credible-region summary.
- Vendor Plotly locally or bundle it for fully offline deployment.
