# PMWGS test script
devtools::load_all()

load(here::here("scratch", "forstmann_long.RData"))
message("Setup")
importance_samples <- 100 # number of importance samples
n_particles <- 10 # number of particles

importance_samples <- is2(sampled, importance_samples, n_particles)

message("Get Maximum Likelihood and Bootstrap for Standard error")

summary_like <- summarise(importance_samples)
print(summary_like)

save.image(here::here("scratch", "IS2_pmwg3.RData"))
