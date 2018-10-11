get_design <- function(years_BP, D14_atm){
    # Arguments: years (BP) of measurements and basis
    # Output: Designmatrix constructed by convolving basis with atmospheric 14C-curve
    taumax <- 50
    tau <- 0:taumax
    X <- numeric()
    for (year_BP in years_BP){
        X <- rbind(X, approx(x = c(-100, D14_atm$year_BP, 1000), 
                                y = c(0, D14_atm$D14C, 0),
                                xout = year_BP + tau)$y)
    }
    X
}


write(
    "model {
    mu <- (1 - f_Old) * X %*% f
    f_Old ~ dunif(0, 1)
    f ~ ddirch(alpha)
    b ~ dnorm(0, 0.00001)
    C ~ dunif(-1000,0)
    for (i in 1:n){
        D14C[i] ~ dnorm(mu[i] + f_Old * C + b * (year_BP[i] - 50), tau[i])
    }
    tau_Old <- -log((C / 1000) + 1) * 8267
    }", 
file = "deconvolve_f.jags")

deconvolve_f <- function(data, D14_atm, file = "deconvolve_f.jags"){
    # Arguments: data (containing columns year_BP, D14 and sd)
    # Output: A JAGS mcarray object
    X <- get_design(data$year_BP, D14_atm)
    alpha <- rep(1, ncol(X))
    jags_data <- list(D14C = data$D14C, X = X, 
                      n = length(data$D14C), 
                      alpha = alpha,
                      tau = 1 / data$sd^2,
                      year_BP = data$year_BP)
    jags_model <- jags.model(file,
                             data = jags_data, 
                             n.chains = 2, n.adapt = 10000, quiet = TRUE)
    update(jags_model, 50000)
    jags_par <- c("f_Old", "C", "mu", "b", "f", "tau_Old")
    output <- jags.samples(jags_model, jags_par, n.iter = 1000000, thin = 10)
    output
}

draw_tau_hist <- function(out_jags, tau.grid){
        f <- out_jags$f[, , 1]  # Pick first chain
        tau <- numeric(ncol(f))
        for (i in 1:ncol(f)){
            tau[i] <- sample(tau.grid, size = 1, prob = f[, i] / sum(f[, i]))
        }
        data.frame(tau = tau) %>% 
            ggplot(aes(x = tau)) + stat_bin(aes(y=..density..), bins = 51, fill = "grey") +
            scale_y_continuous(expand = c(0, 0), limits = c(NA, .08)) + theme_bw() +
            geom_hline(yintercept=1/50, color="red")+
            ylab("Posterior predictive density") + xlab(expression(paste(tau, " (years)")))
        
            
}

plot_fit <- function(out_jags, data, D14_atm){
    out.data <- data.frame()
    f_Old <- out_jags$f_Old[1,,1]
    b <- out_jags$b[1,,1]
    C <- out_jags$C[1,,1]
    f <- t(out_jags$f[,,1])
    years_BP <- min(data$year_BP):max(data$year_BP)
    mu <- (f * matrix(rep((1 - f_Old), ncol(f)), ncol = ncol(f))) %*% t(get_design(years_BP, D14_atm)) + outer(b, (years_BP - 50)) + outer(C * f_Old, rep(1,length(years_BP)))
    mu.mean <- apply(mu, 2, mean)
    mu.95 <- apply(mu, 2, function(x){quantile(x, probs = .95)})
    mu.05 <- apply(mu, 2, function(x){quantile(x, probs = .05)})
    out.data <- data.frame(year_cal = (1950 - years_BP),
                        mu.mean = mu.mean,
                        mu.95 = mu.95,
                        mu.05 = mu.05)
    ggplot() + geom_errorbar(data = data, aes(x = year_cal, ymax=D14C+3*sd, ymin=D14C-3*sd), width = 2)  +
        geom_point(data = data, aes(x = year_cal, y = D14C), size = 2, fill = "white") +
        geom_line(data = out.data, aes(x = year_cal, y=mu.mean)) + 
        geom_ribbon(data = out.data, aes(x=year_cal, ymax = mu.95, ymin = mu.05), alpha = .2) +
        labs(x = "Calendar year", y = expression(paste(Delta ^14, "C", " (\u2030)"))) + 
        theme_bw() +
        coord_flip()
}

table_pars <- function(out_jags){
    pars <- data.frame(f_Old = out_jags$f_Old[1, ,1],
                       b = out_jags$b[1, ,1],
                       tau_Old = out_jags$tau_Old[1,,1])
    gather(pars, parameter, value) %>% 
        group_by(parameter) %>% 
        summarise(mean(value), sd(value)) %>% 
        ungroup() %>% 
        rename(mean=`mean(value)`, sd=`sd(value)`) %>% 
        mutate(`mean (sd)`=paste(round(mean, 2), " (", round(sd, 2), ")", sep = "")) %>% 
        select(parameter, `mean (sd)`, ) %>% 
        knitr::kable(caption = "Posterior mean and standard deviaions of parameter posterior distributions")
}

