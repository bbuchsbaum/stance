# Basic example for Continuous Bayesian Decoder
library(stance)

# simulate small dataset
sim <- simulate_fmri_data(V = 200, T = 100, K = 3, algorithm = "CBD")

# initialize decoder
cbd <- ContinuousBayesianDecoder$new(
  Y = sim$Y,
  K = 3,
  r = 10,
  hrf_basis = "canonical"
)

# fit model
cbd$fit(max_iter = 50)

# extract results
maps <- cbd$get_spatial_maps()
states <- cbd$get_state_sequence()

