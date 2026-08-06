# ABC_Visualisation# ABC / ABC-MCMC SIR Dashboard

## 1. What this is

This is a self-contained browser dashboard for likelihood-free parameter inference in a stochastic SIR epidemic simulator. It compares independent ABC rejection sampling with ABC-MCMC and lets a visitor see how the tolerance, summary statistics, and random-walk proposal scales change the resulting approximate posterior.

The implementation is intended as an interactive statistics portfolio piece: the simulator, binomial sampler, rejection sampler, ABC-MCMC transition, validation code, and comparison metrics are written explicitly in vanilla JavaScript rather than delegated to an ABC or MCMC package.

## 2. The problem being solved

The SIR model divides a fixed population into three compartments:

- `S`: susceptible people
- `I`: currently infected people
- `R`: recovered people

The two inferred parameters are:

- `beta`: transmission rate. Larger values generally produce faster epidemic growth.
- `gamma`: recovery rate. Larger values generally shorten infectious durations and remove people from `I` more quickly.

The dashboard performs inference by simulation: propose `beta` and `gamma`, generate a stochastic epidemic, reduce the trajectory to selected summary statistics, and compare those summaries with a fixed observed outbreak.

### Likelihood caveat

The page deliberately uses an ABC workflow and never evaluates a likelihood. This matches the simulator-as-a-black-box setting for which ABC is commonly useful, especially when realistic epidemic simulators contain latent structure, partial observation, or mechanisms whose likelihood is unavailable or expensive.

For this exact teaching model, however, ABC is **not mathematically necessary** in the strongest sense. Conditional on the complete daily `S`, `I`, and `R` path, the implemented transitions are binomial and a complete-path likelihood could be written down. The dashboard should therefore be understood as a controlled demonstration of likelihood-free inference and summary-statistic loss, not as proof that this small SIR model has no tractable likelihood under any observation scheme.

The "observed data" are also synthetic, not epidemiological measurements. They are one stochastic epidemic generated from known parameters and embedded permanently in the HTML file.

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

Infections and recoveries are therefore sampled independently conditional on the current state. A person infected during the current daily step cannot recover until a later step because recoveries are drawn from the original `I`, not from `I + new_infections`.

### Fixed simulation settings

| Setting | Value |
|---|---:|
| Population `N` | 1,000 |
| Initial infected `I0` | 5 |
| Initial recovered `R0` | 0 |
| Initial susceptible `S0` | 995 |
| Horizon `T` | 100 days |
| Time step `dt` | 1 day |

### Priors

The two priors are independent uniforms:

```text
beta  ~ Uniform(0.05, 1.00)
gamma ~ Uniform(0.02, 0.50)
```

### Ground truth and fixed observed epidemic

```text
beta_true  = 0.4
gamma_true = 0.1
beta_true / gamma_true = 4
observed-generation seed = 4041001
```

The observed epidemic is one literal, hardcoded stochastic realization. It was generated once with the same SIR simulator and seed `4041001`, then its `S`, `I`, `R`, `newInfections`, and `newRecoveries` arrays were embedded in the page.

It is **not regenerated on page load**, when validation is rerun, or when a new inference seed is requested. Its displayed summaries are:

```text
peak infected = 410
peak day      = 24
final size    = 980
```

Here, final epidemic size is `N - S(T)`, so it counts everyone infected at least once by day 100.

## 4. Summary statistics

The visitor can choose between two summary vectors.

### Rich summary

```text
[peak infection height, peak timing, final epidemic size]
= [max(I), index of first maximum of I, N - S(T)]
```

For the fixed observed epidemic, the rich vector is:

```text
[410, 24, 980]
```

### Poor summary

```text
[final epidemic size]
= [N - S(T)]
```

For the fixed observed epidemic, the poor vector is:

```text
[980]
```

### Distance metric

The implementation uses normalized Euclidean distance. For a simulated epidemic:

```text
peak_difference   = (simulated_peak - 410) / N
timing_difference = (simulated_peak_day - 24) / T
final_difference  = (simulated_final_size - 980) / N

rich_distance = sqrt(
  peak_difference^2
  + timing_difference^2
  + final_difference^2
)

poor_distance = abs(final_difference)
```

With `N = 1000` and `T = 100`, peak heights and final sizes are scaled by 1,000 while peak timing is scaled by 100 days. This prevents the count-valued components from dominating the distance only because they have larger numerical units.

This is a fixed range-based normalization, not standardization by simulated or posterior standard deviations. It gives each component a comparable nominal scale but does not guarantee equal inferential importance.

## 5. Panel-by-panel explanation

### Observed epidemic panel

**What it shows**

- The fixed `S(t)`, `I(t)`, and `R(t)` trajectories for days 0 through 100.
- The infection peak marked at `I = 410` on day 24.
- The rich and poor summary vectors used by both ABC methods.
- The offline provenance seed `4041001`.

**What to look at**

The summaries intentionally discard most of the full trajectory. Compare the rich vector with the poor vector: the poor version retains only the final total and cannot tell whether an epidemic was fast and sharp or slow and prolonged.

**Do not over-interpret**

This is one stochastic realization. Even at the true parameters, another simulated epidemic will generally have a different peak, timing, and final size.

### Panel 1: ABC rejection sampling

**Algorithm actually run**

1. Draw 4,000 independent parameter pairs from the uniform priors.
2. Simulate one stochastic SIR trajectory for each pair.
3. Cache each proposal's `beta`, `gamma`, full `I(t)` curve, peak height, peak day, and final size.
4. Compute the selected normalized summary distance.
5. Accept a proposal when `distance < epsilon`.

The default tolerance is `epsilon = 0.080`. The slider spans `0.005` to `0.250` in increments of `0.005`.

**Caching and interaction**

- The 4,000 epidemics are generated once per inference seed.
- Moving the epsilon slider only re-filters cached distances.
- Switching between rich and poor summaries recomputes distances from cached summaries; it does not resimulate epidemics.
- The curve overlay displays the eight accepted proposals with the smallest current distances.

**What to try**

- Lower epsilon and watch the accepted parameter cloud shrink while the acceptance rate falls.
- Switch to the poor summary and observe the broader ridge of parameter pairs that can reproduce a similar final size.
- Compare the eight accepted infected curves with the fixed observed `I(t)` trajectory.

**Do not over-interpret**

Accepted rejection draws are independent conditional on the fixed simulation batch, but there are only 4,000 proposals. At small epsilon, the apparent posterior can be sparse, noisy, or empty because of finite Monte Carlo budget rather than because the compatible region is truly absent.

### Panel 2: ABC-MCMC

**Algorithm actually run**

The chain uses 6,500 random-walk transitions and discards the first 500 states as burn-in, leaving 6,000 retained states.

Initialization repeatedly draws from the prior, simulates an epidemic, and stops at the first draw with `distance < epsilon`. It gives up after 25,000 attempts and asks the visitor to increase epsilon if no valid initial state is found.

At iteration `t`, the proposal is:

```text
beta_proposed  = beta_current  + sigma_beta  * Normal(0, 1)
gamma_proposed = gamma_current + sigma_gamma * Normal(0, 1)
```

Default proposal scales:

```text
sigma_beta  = 0.045
sigma_gamma = 0.018
```

The sliders allow:

```text
sigma_beta:  0.005 to 0.150, step 0.005
sigma_gamma: 0.002 to 0.080, step 0.002
```

A proposal outside either prior range is rejected immediately and does not call the SIR simulator. An in-bounds proposal gets one simulated epidemic and is accepted exactly when:

```text
distance(proposed simulation, observed summaries) < epsilon
```

Otherwise the next chain state is an exact copy of the current state.

**Relation to the textbook ABC-MCMC rule**

The code implements the hard-threshold ABC-MCMC rule associated with Marjoram-style ABC-MCMC, but with its Metropolis-Hastings ratio simplified away. That simplification is valid for this implementation because:

- both priors are uniform within fixed bounds;
- the Gaussian random-walk proposal is symmetric before applying the prior support check;
- therefore the prior and proposal ratios equal one for an in-bounds move.

The code does not accept a proposal because its distance is better than the current distance. It accepts any in-bounds proposal whose newly simulated data pass the fixed tolerance, as the ABC-MCMC rule requires under these assumptions.

**What to try**

- Inspect the path plot to see the chain move within the rejection-compatible region.
- Inspect the beta and gamma traces for long repeated stretches, slow drift, or poor mixing.
- Change the proposal scales and rerun the chain. Parameter, summary, and epsilon changes mark the existing chain as stale rather than silently altering it.

**Do not over-interpret**

- The 6,000 retained states are dependent and include exact repeats after rejected proposals.
- A high move acceptance rate is not automatically good; very small steps can create a highly autocorrelated chain.
- A low move acceptance rate can result from oversized steps, out-of-bounds proposals, or simulated epidemics that miss the tolerance.
- The page provides traces and simple summaries, but no effective sample size, autocorrelation estimate, R-hat, or multi-chain convergence assessment.

### Panel 3: Matched comparison

**What "matched settings" means**

The comparison uses the completed MCMC chain's frozen summary choice and epsilon. It re-filters Panel 1's cached rejection proposals under those same settings, so both methods are compared against the same tolerance-defined ABC target.

If the visitor changes the summary, epsilon, or proposal scales after a chain finishes, the comparison is marked stale and continues to describe the last completed chain until the chain is rerun.

**What the combined figure shows**

- Accepted rejection points and post-burn-in MCMC states in `(beta, gamma)` space.
- The true parameter `(0.4, 0.1)` and each method's sample mean.
- Separate beta and gamma marginal histograms.

The plot displays at most 1,800 evenly spaced points from each sample for browser performance, but the numerical summaries use all accepted rejection draws and all retained MCMC states.

**Efficiency metrics**

- **Rejection yield per 1,000 simulations**: `1000 * accepted rejection proposals / 4000`.
- **MCMC move yield per 1,000 simulator calls**: `1000 * accepted MCMC moves / (initialization attempts + in-bounds proposal simulations)`.
- **Simulator budget**: 4,000 for rejection; for MCMC, initialization simulations plus in-bounds transition simulations. Out-of-bounds proposals are not counted because they do not run the SIR simulator.
- **Central 80% footprint**: the area of the axis-aligned rectangle formed by the marginal 10th-to-90th percentile interval of beta times the marginal 10th-to-90th percentile interval of gamma, divided by the full prior rectangle area.

**What to look at**

ABC-MCMC often produces more accepted moves per simulator call because it proposes locally after reaching the compatible region, whereas rejection sampling continues drawing across the entire prior. Under the poor summary, both methods should admit a much wider, ridge-like set of beta and gamma combinations.

**Do not over-interpret**

- Accepted MCMC moves are not independent posterior samples, so move yield is not an effective-sample-size measure.
- The central 80% footprint is not a highest-posterior-density region. It is an axis-aligned marginal rectangle and can misrepresent a diagonal or curved ridge.
- A narrower MCMC cloud is not automatically evidence of greater accuracy; it may indicate autocorrelation, inadequate mixing, or failure to explore all compatible regions.

### Simulator validation and source panels

The validation section compares 500 stochastic trajectories with a deterministic RK4 SIR reference in each of three regimes:

| Regime | beta | gamma | Validation seed |
|---|---:|---:|---:|
| Below threshold | 0.08 | 0.20 | 20260806 |
| Moderate outbreak | 0.25 | 0.10 | 20261806 |
| Ground-truth regime | 0.40 | 0.10 | 20262806 |

It checks population conservation, non-negative integer states, non-increasing `S`, non-decreasing `R`, and the ensemble mean's normalized RMSE from the deterministic infected trajectory. The page also shows 10th-to-90th percentile bands and one example stochastic path.

The expandable implementation section prints the exact JavaScript functions and configuration objects executed by the page. It is included for transparency; changing the displayed text alone would not alter the already-defined functions unless the HTML source itself is edited.

## 6. Implementation notes

### Client-side architecture

- Everything runs in the browser; there is no backend and no external data file.
- The only external runtime dependency is Plotly.js `2.35.2`, loaded from a CDN.
- The SIR, PRNG, binomial sampler, RK4 solver, ABC rejection sampler, and ABC-MCMC chain are implemented in vanilla JavaScript.
- Because Plotly is loaded remotely, the HTML needs network access for charts unless Plotly is vendored locally. Numerical code can still exist in the page, but chart rendering reports a CDN load failure when Plotly is unavailable.

### Caching and performance

- Rejection sampling stores 4,000 parameter pairs, raw summaries, and all 4,000 infected trajectories of length 101.
- The rejection batch runs in chunks of 200 and yields to the browser between chunks.
- Epsilon changes re-filter the cache; summary changes recompute only distances.
- ABC-MCMC must be rerun when its summary, epsilon, or proposal scales change because those settings alter the transition history.
- MCMC plots and status are refreshed every 250 transitions rather than after every transition.
- The comparison panel is derived from the rejection and MCMC caches and performs no new epidemic simulations.

### Reproducibility

- All simulation randomness uses the explicit `mulberry32` seeded PRNG once a seed has been selected.
- The default inference base seed is `314159265`.
- The base seed is stored under `sessionStorage["sir-abc-rejection-seed-v1"]`, so reloads in the same browser tab/session reproduce the same rejection batch.
- The MCMC seed is derived deterministically as `base_seed XOR 0x9E3779B9`. With the default base seed, this is `2358167832`.
- "New inference seed" uses `crypto.getRandomValues` when available, with a time/`Math.random` fallback, stores the result in `sessionStorage`, and reruns rejection sampling and MCMC.
- The fixed observed epidemic never changes with the inference seed.
- Validation uses its own fixed seeds and therefore reproduces the same validation ensembles on each rerun.

### Known limitations and simplifications

- **ABC is pedagogical for this exact toy model.** The page treats the simulator as likelihood-free, but the complete daily binomial transition likelihood could be written for the fully observed SIR path.
- **Finite budgets matter.** Rejection uses 4,000 proposals; MCMC uses one 6,500-transition chain. Sparse regions, sample means, widths, and apparent modes can change with the seed.
- **Extreme epsilon values behave differently.** Very small epsilon can produce zero rejection draws or make MCMC initialization fail after 25,000 attempts. Very large epsilon accepts weak matches and makes the ABC approximation approach the prior rather than a sharply informed posterior.
- **Summaries lose information.** Even the rich vector reduces a 101-day trajectory to three numbers. The poor vector intentionally creates substantial non-identifiability.
- **No observation-error model is included.** The synthetic trajectory is treated as the target directly; there is no reporting noise, under-ascertainment, missingness, or measurement process.
- **One simulation per proposed parameter is used.** There is no averaging over repeated simulated datasets at the same parameter and no synthetic-likelihood estimate.
- **The chain is not adaptively tuned.** Proposal scales are fixed during a run, burn-in is always 500, and there is no automatic covariance adaptation.
- **No formal MCMC convergence diagnostics are computed.** The dashboard uses one chain, traces, acceptance rates, and sample means only.
- **The custom binomial sampler has two numerical branches.** For reflected `p <= 0.5`, it uses inverse-CDF recursion when `n * p < 30`; otherwise it uses tangent-envelope rejection sampling with a Lanczos approximation to `logGamma`. This avoids a library call and avoids `n` Bernoulli draws, but the dashboard's epidemic-level validation is not a formal exhaustive goodness-of-fit test of the sampler for every possible `(n, p)`.
- **The deterministic comparison is only a reference.** RK4 uses an internal step of `0.05` days for the continuous-time ODE. The stochastic simulator uses one-day discrete transitions, so its ensemble mean is not expected to match the ODE exactly.
- **Peak timing uses the first maximum.** If `I(t)` has a tied plateau, the earliest day is recorded.
- **Final size is evaluated at day 100.** If infections remain at day 100, `N - S(T)` includes all infections to date but is not necessarily the eventual post-extinction epidemic size.
- **The comparison footprint ignores dependence and geometry.** It uses marginal quantiles and an axis-aligned area, not a joint credible region or effective sample size.

## 7. How to extend this

- Run multiple ABC-MCMC chains and add R-hat, effective sample size, and autocorrelation diagnostics.
- Add adaptive or pilot-tuned proposal covariance, including correlated beta/gamma proposals.
- Compare hard-threshold ABC with kernel-weighted ABC, ABC-SMC, or sequential Monte Carlo.
- Add alternative summaries such as early growth rate, epidemic duration, several time-point counts, or learned summaries.
- Add a realistic observation model for noisy, delayed, or partially observed case counts.
- Compare the ABC approximation with an exact or particle-filter likelihood method for this small SIR model.
- Infer additional parameters or initial conditions, or replace SIR with SEIR or another simulator.
- Add posterior predictive checks based on fresh simulations from retained parameter draws.
- Vendor Plotly locally or bundle the JavaScript for fully offline deployment.
