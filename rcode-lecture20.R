################################################################################
# LECTURE 19 R CODE
# CIPOLLI
# MATH 240 - SPRING 2025
################################################################################
library(tidyverse)
library(patchwork)

##########################################################
# Plot the distribution (Chapter 6)
##########################################################
n <- 50         # sample size
mu0 <- 5        # null mean
mu.lower  <- 4  # lower population mean
mu.higher <- 6  # higher population mean
sigma <- 1.3    # population standard deviation

# For plotting the distribution
ggdat.t.data <- tibble(t.lower  = numeric(1000),
                       t.null   = numeric(1000),
                       t.higher = numeric(1000))
# Simulation
for(i in 1:1000){
  ggdat.t.data[i,] <- tibble(x.lower = rnorm(n, mean = mu.lower, sd=sigma),
                             x.null = rnorm(n, mean = mu0, sd=sigma),
                             x.higher = rnorm(n, mean = mu.higher, sd=sigma)) |>
    summarize(t.lower = (mean(x.lower) - mu0)/(sd(x.lower)/sqrt(n)),
              t.null = (mean(x.null) - mu0)/(sd(x.null)/sqrt(n)),
              t.higher = (mean(x.higher) - mu0)/(sd(x.higher)/sqrt(n)))
}

# For plotting the *actual* sampling distributions
ggdat.t <- tibble(t=seq(-10,10,length.out=1000))|>
  mutate(pdf.null = dt(t, df=n-1),
         pdf.lower = dt(t, df=n-1, ncp = (mu.lower-mu0)/(sqrt(sigma)/sqrt(n))),
         pdf.higher = dt(t, df=n-1, ncp = (mu.higher-mu0)/(sqrt(sigma)/sqrt(n))))

# Create Plot
ggplot()+
  # Plot the actual distributions
  geom_line(data=ggdat.t, aes(x=t, y=pdf.null, 
                              color="T-distribution (Null)"))+
  geom_line(data=ggdat.t, aes(x=t, y=pdf.lower, 
                              color="T-distribution (Lower)"))+
  geom_line(data=ggdat.t, aes(x=t, y=pdf.higher, 
                              color="T-distribution (Higher)"))+
  # Plot the simulated data
  stat_density(data=ggdat.t.data, 
               aes(x=t.null, color="T-distribution (Null)"),
               geom="line", linetype="dotted", linewidth=1)+
  stat_density(data=ggdat.t.data, 
               aes(x=t.lower, color="T-distribution (Lower)"),
               geom="line", linetype="dotted", linewidth=1)+
  stat_density(data=ggdat.t.data, 
               aes(x=t.higher, color="T-distribution (Higher)"),
               geom="line", linetype="dotted", linewidth=1)+
  # Clean up aesthetics
  geom_hline(yintercept=0)+
  theme_bw()+
  xlab("t")+
  ylab("Density")+
  labs(color="",
       caption="Dotted lines represent simulated T")+
  ggtitle("Sampling Distributions for T",
          subtitle=bquote(mu[X]==5*","~n==50*","~mu[lower]==4*","~mu[higher]==6))


##########################################################
# Plot the data (Chapter 5)
##########################################################
dat.SS <- read_csv(file = "socsecRC.csv")

ggplot(data=dat.SS) +                                         # Specify data to plot
  geom_histogram(aes(x=White_SS, y = after_stat(density)),    # Density histogram
                 breaks=seq(1,6, length.out=10),
                 color="grey30", fill="lightgray")+             # Specify bins
  geom_density(aes(x=White_SS), color="darkred")+             # Add density estimate
  geom_hline(yintercept=0)+                                   # Add x-axis
  theme_bw() +                                                # Remove gray background
  xlab("Sample Mean Perceived Whiteness") +                   # x-axis label
  ylab("Density")                                             # y-axis label


library(e1071)
dat.SS |>
  summarize(mean       = mean(White_SS),
            sd         = sd(White_SS),
            median     = median(White_SS),
            IQR        = IQR(White_SS),
            skewness   = skewness(White_SS),
            exkurtosis = kurtosis(White_SS))

##########################################################
# Compute the test statistic
##########################################################
mu0 <- 3.5
x <- dat.SS$White_SS
(xbar <- mean(x))
(s <- sd(x))
(n <- length(x))
any(is.na(x)) # no missing data
(t.stat <- (xbar - mu0)/(s/sqrt(n)))

library(effectsize)
hedges_g(x = x, mu = mu0, alternative = "two.sided")
interpret_hedges_g(1.34)


# t.test() function does a lot of the heavy lifting
t.test(x=x, mu = mu0, alternative = "two.sided")


##########################################################
# Plot it
##########################################################
# For plotting the null distribution
ggdat.t <- tibble(t=seq(-5,5,length.out=1000))|>
  mutate(pdf.null = dt(t, df=n-1))
# For plotting the observed point
ggdat.obs <- tibble(t    = t.stat, 
                    y    = 0) # to plot on x-axis

# Resampling to approximate the sampling distribution 
# on the data
R <- 1000
resamples <- tibble(t=numeric(R))
for(i in 1:R){
  curr.sample <- sample(x=x,
                        size=n,
                        replace=T)
  resamples$t[i] = (mean(curr.sample)-mu0)/(sd(curr.sample)/sqrt(n))
}

t.breaks <- c(-5, qt(0.025, df = n-1), # rejection region (left)
              0, 
              qt(0.975, df = n-1), 5,  # rejection region (right)
              t.stat)                  # t-statistic observed
xbar.breaks <- t.breaks * s/(sqrt(n)) + mu0

# Create Plot
ggplot() +
  # null distribution
  geom_line(data=ggdat.t, 
            aes(x=t, y=pdf.null))+
  geom_hline(yintercept=0)+
  # rejection regions
  geom_ribbon(data=subset(ggdat.t, t<=qt(0.025, df=n-1)), 
              aes(x=t, ymin=0, ymax=pdf.null),
              fill="grey", alpha=0.5)+
  geom_ribbon(data=subset(ggdat.t, t>=qt(0.975, df=n-1)), 
              aes(x=t, ymin=0, ymax=pdf.null),
              fill="grey", alpha=0.5)+
  # plot p-value (not visible)
  geom_ribbon(data=subset(ggdat.t, t>=t.stat), 
              aes(x=t, ymin=0, ymax=pdf.null),
              fill="reg", alpha=0.25)+
  # plot observation point
  geom_point(data=ggdat.obs, aes(x=t, y=y), color="red")+
  # Resampling Distribution
  stat_density(data=resamples, 
               aes(x=t),
               geom="line", color="grey")+
  # clean up aesthetics
  theme_bw()+
  scale_x_continuous("t",
                     breaks = round(t.breaks,2),
                     sec.axis = sec_axis(~.,
                                         name = bquote(bar(x)),
                                         breaks = t.breaks,
                                         labels = round(xbar.breaks,2)))+
  ylab("Density")+
  ggtitle("T-Test for Mean Perceived Whiteness of Social Security Recipients",
          subtitle=bquote(H[0]==3.5*";"~H[a]!=3.5))
