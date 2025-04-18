library(pwr)
library(tidyverse)
library(patchwork)
###Question 1: Use pwr to find # of obs needed to detect mod to lrg effect
num.obs <- pwr.t.test(n=NULL, d = 0.65, sig.level = 0.05, power = 0.80, type="one.sample")

###Question 2: Pull Figure Data
fig2data <- read.csv("figure2data.csv")|>
  rename(
    closer.vals = X0.275714146,
    further.vals = X.0.191685046
  )|>
mutate(difference = fig2data$closer.vals - fig2data$further.vals)

##Question 3: Summarize the data
#Further Data
further.sum <- fig2data |>
  summarize(
    mean = mean(further.vals, na.rm = TRUE),
    sd = sd(further.vals, na.rm = TRUE)
  )

#closer Data
closer.sum <- fig2data |>
  summarize(
    mean = mean(further.vals, na.rm = TRUE),
    sd = sd(further.vals, na.rm = TRUE)
  )

diff.sum <- fig2data |>
  summarize(
    mean = mean(difference, na.rm = TRUE),
    sd = sd(difference, na.rm = TRUE)
  )

data.long <- fig2data |>
  pivot_longer(
    cols = c(closer.vals, further.vals, difference),
    names_to = "measurement_type",
    values_to = "fluorescence_change"
  )

# Change the labels
data.long <- data.long |>
  mutate(measurement_type = factor(measurement_type, 
                                   levels = c("closer.vals", "further.vals", "difference"),
                                   labels = c("Closer Values", "Further Values", "Difference")))

graphical.sums <- ggplot(data.long, aes(x = measurement_type, y = fluorescence_change, fill = measurement_type)) +
  geom_boxplot() +
  labs(
    title = "Box Plots of Fluorescence Changes",
    x = "Measurement Type",
    y = "Change in Fluorescence",
    fill = "Measurement Type"
  ) +
  theme_minimal() 

#Question 4: Conducting Statistical Inferences

mu0 = 0
x <- fig2data$closer.vals
y <- fig2data$further.vals
z <- fig2data$difference
(xbar <- mean(x))
(s <- sd(x))
(n <- length(x))
any(is.na(x)) # no missing data


library(effectsize)
x.hedges <- hedges_g(x = x, mu = mu0, alternative = "two.sided")
interpret_hedges_g(1.34)

y.hedges <- hedges_g(x = y, mu = mu0, alternative = "two.sided")
interpret_hedges_g(1.34)

z.hedges <- hedges_g(x = z, mu = mu0, alternative = "two.sided")
interpret_hedges_g(1.34)

close.t.stat <- t.test(x=x, mu = mu0, alternative = "two.sided")
far.t.stat <- t.test(x=x, mu = mu0, alternative = "two.sided")
diff.t.stat <- t.test(x=x, mu = mu0, alternative = "two.sided")
(x.p.val <- 2*pt(q=-abs(close.t.stat), df = n-1))
(y.p.val <- 2*pt(q=-abs(far.t.stat), df = n-1))
(z.p.val <- 2*pt(q=-abs(diff.t.stat), df = n-1))
# t.test() function does a lot of the heavy lifting

#Question 5: Hypothesis test plots

plot_t_test <- function(vec, mu0 = 0, title_text) {
  # compute basic quantities
  n      <- length(vec)
  s      <- sd(vec)
  xbar   <- mean(vec)
  t.stat <- (xbar - mu0) / (s / sqrt(n))
  df     <- n - 1
  alpha  <- 0.05
  
  # null distribution curve
  ggdat.t <- tibble(t = seq(-5, 5, length.out = 1000)) %>%
    mutate(pdf.null = dt(t, df = df))
  # observed point
  ggdat.obs <- tibble(t = t.stat, y = 0)
  
  # bootstrap‐resampled t’s
  R <- 1000
  resamples <- tibble(t = replicate(R, {
    samp <- sample(vec, size = n, replace = TRUE)
    (mean(samp) - mu0) / (sd(samp) / sqrt(n))
  }))
  
  # critical values
  lower.q  <- qt(alpha/2, df = df)
  upper.q  <- qt(1 - alpha/2, df = df)
  t.breaks <- c(lower.q, 0, upper.q, t.stat)
  xbar.breaks <- t.breaks * (s / sqrt(n)) + mu0
  
  # build the plot
  ggplot() +
    # null‐curve
    geom_line(data = ggdat.t, aes(x = t, y = pdf.null)) +
    geom_hline(yintercept = 0) +
    
    # rejection regions (two‐tailed)
    geom_ribbon(data = filter(ggdat.t, t <= lower.q),
                aes(x = t, ymin = 0, ymax = pdf.null),
                fill = "grey", alpha = 0.5) +
    geom_ribbon(data = filter(ggdat.t, t >= upper.q),
                aes(x = t, ymin = 0, ymax = pdf.null),
                fill = "grey", alpha = 0.5) +
    
    # p‐value shading
    geom_ribbon(data = filter(ggdat.t, t <= lower.q | t >= upper.q),
                aes(x = t, ymin = 0, ymax = pdf.null),
                fill = "red", alpha = 0.25) +
    
    # observed t
    geom_point(data = ggdat.obs, aes(x = t, y = y),
               color = "red", size = 3) +
    
    # bootstrap density
    stat_density(data = resamples, aes(x = t),
                 geom = "line", linetype = "dashed", size = 1) +
    
    # axes & titles
    theme_bw() +
    scale_x_continuous("t",
                       breaks = t.breaks,
                       sec.axis = sec_axis(
                         ~ . * (s / sqrt(n)) + mu0,
                         name    = bquote(bar(x)),
                         breaks  = t.breaks,
                         labels  = round(xbar.breaks, 2)
                       )
    ) +
    ylab("Density") +
    ggtitle(title_text,
            subtitle = bquote(H[0] == .(mu0) * ";" ~~ H[a] != .(mu0)))
}


p_closer  <- plot_t_test(fig2data$closer.vals, mu0 = 0,
                         title_text = "T‐Test for Closer Values")
p_further <- plot_t_test(fig2data$further.vals, mu0 = 0,
                         title_text = "T‐Test for Further Values")
p_diff    <- plot_t_test(fig2data$difference, mu0 = 0,
                         title_text = "T‐Test for Difference")

# Stack them vertically (or arrange however you like)
combined_plot <- (p_closer / p_further / p_diff) +
  plot_layout(ncol = 1, guides = "collect")

# Print to the device
print(combined_plot)

