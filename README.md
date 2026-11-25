## ecDNA–ctDNA modelling code

Julia code accompanying the ecDNA/cfDNA tumour modelling project.

## Files

- `model.jl`: Baseline stochastic tumour-growth model with ecDNA copy-number-dependent selection, ecDNA segregation, ctDNA/ecDNA shedding, blood clearance, resampling, and longitudinal 15 ml blood sampling.

- `treatment.jl`: Extends the baseline model to include ecDNA-targeted chemotherapy. During treatment, the death rate of a cell carrying `n` ecDNA copies is `d(n) = d0 + D*n`. Treatment starts after the first ecDNA-positive 15 ml blood sample and lasts for 30 days. This can be easily modified in the code.

- `ABC_selection.jl`: Approximate Bayesian computation for distinguishing neutral from positive ecDNA selectio, estimating tumour burden at the first blood draw (`N_first`) and estimating the net growth rate (`r`). ABC uses longitudinal blood measurements only.

- `ABC_treatment.jl`: ABC inference of treatment strength `D` from early longitudinal blood measurements over the first 48 hours after treatment initiation. Treatment strength is also classified as Weak, Intermediate, or Strong.

- `sequencing.jl`: Sequencing/observation model used to assess the effect of limited molecular information and ecDNA detection thresholds on treatment inference. Particularly, we use effective molecular depth `N_target` and ecDNA detection thresholds `L`.

## Main parameter priors

For the ABC analyses:

- `r ~ Uniform(0.02, 0.06)`
- `log10(N_first) ~ Uniform(7, 9)`
- positive-selection model: `s_max ~ Uniform(1.05, 3.0)`
- treatment model: `D ~ Uniform(0, 0.012)`

Treatment inference fixes `s_max = 1.2`.

## Treatment classes

- Weak: `D < 0.003`
- Intermediate: `0.003 <= D < 0.0075`
- Strong: `D >= 0.0075`

## Running the code

The files depend on one another in the following order:

`model.jl` → `treatment.jl` → `ABC_treatment.jl` → `sequencing.jl`

Each ABC file also includes the model files it requires, so they can be run individually.

Small reference tables can be used for testing. The final treatment analysis used 5000 evaluable reference particles with the closest 50 particles retained. Of course these numbers can be increased for deeper accuracy of the ABC estimate.

## Notes and Requirements

- ABC distances are based only on longitudinal blood measurements (sampled ecDNA and tumour-derived chromosomal DNA). Tumour size, mean ecDNA copy number, and other latent tumour states are not included in the ABC distance.

- The sequencing model is an adaptation of a finite-molecular-information observation model and uses Binomial sampling of realised molecular depth followed by Multinomial sampling of ecDNA/chromosomal composition.

- The Figures in the report are produced from the functions in this repository. The actual codes for reproducing the visualisations are not attached.

- The only external package required is Distributions. It can be installed using:

  `using Pkg`
  
  `Pkg.add("Distributions")`.

  Random, Statistics and Serialization are part of the Julia standard library.
