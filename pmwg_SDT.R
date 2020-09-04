require(pmwg)

load("wagenmakers2008.RData")
wagenmakers2008 <- data
data <- as.data.frame(table(data$subject, data$prop, data$freq, data$resp))
names(data) <- c("subject", "prop", "freq", "resp", "n")


SDT_loglike_fast <- function(x, data) {
  out <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    if (data$prop[i] == "w") {
      if (data$freq[i] == "hf") {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["HF.d"], sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["HF.d"], sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      } else if (data$freq[i] == "lf") {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["LF.d"], sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["LF.d"], sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      } else if (data$freq[i] == "vlf") {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["VLF.d"], sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = x["VLF.d"], sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      } else {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = 0, sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.w"], mean = 0, sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      }
    }
    else {
      if (data$freq[i] == "hf") {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["HF.d"], sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["HF.d"], sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      } else if (data$freq[i] == "lf") {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["LF.d"], sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["LF.d"], sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      } else if (data$freq[i] == "vlf") {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["VLF.d"], sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = x["VLF.d"], sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      } else {
        if (data$resp[i] == "W") {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = 0, sd = 1, log.p = TRUE, lower.tail = FALSE)
        } else {
          out[i] <- data$n[i] * pnorm(x["C.nw"], mean = 0, sd = 1, log.p = TRUE, lower.tail = TRUE)
        }
      }
    }
  }
  sum(out)
}


# Specify the parameters and priors -------------------------------------------
pars <- c("C.w", "C.nw", "HF.d", "LF.d", "VLF.d")

# Create the Particle Metropolis within Gibbs sampler object ------------------

sampler <- pmwgs(
  data = data,
  pars = pars,
  ll_func = SDT_loglike_fast
)

sampler <- init(sampler, particles=100)

burned <- run_stage(sampler, stage = "burn", iter = 200, particles = 30)
matplot(t(burned$samples$theta_mu), type = "l")

adapted <- run_stage(burned, stage = "adapt", iter = 200, particles = 30)

sampled <- run_stage(adapted, stage = "sample", iter = 200, particles = 20)

matplot(t(sampled$samples$alpha[,8,]), type = "l")
abline(v=200)
abline(v=sampled$samples$idx - 200)
