library(pwr)
library(tidyverse)
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
(xbar <- mean(x))
(s <- sd(x))
(n <- length(x))
any(is.na(x)) # no missing data


library(effectsize)
hedges_g(x = x, mu = mu0, alternative = "two.sided")
interpret_hedges_g(1.34)

(p.val <- 2*pt(q=-abs(t.stat), df = n-1))

# t.test() function does a lot of the heavy lifting
t.stat <- t.test(x=x, mu = mu0, alternative = "two.sided")

