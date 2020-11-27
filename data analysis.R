#
## I Can See Myself Helping: The Effect of Self-Awareness on Prosocial Behavior
## Jerome Lewis, Zachary Himmelberger, and J. Dean Elmore
#

# importing libraries
library(psych)
library(effsize)
library(ggpubr)
library(tidyverse)
library(brms)

# importing data
url.file <- "https://raw.githubusercontent.com/ZacharyHimmelberger/I-can-see-myself-helping/master/Raw%20Data.csv"
df.raw <- read.csv(url.file)

# separating condition into two columns 
# helps with 'readability' but is not strictly necessary
df <- df.raw %>% mutate(experimental = ifelse(df.raw$condition == "E", 1, 0),
                        control = ifelse(df.raw$condition == "C", 1, 0))

# descriptive statistics for age
psych::describe(df$age)

# gender breakdown
table(df$gender)

# Main Analysis

# visualizing priors
sample.priors <- brm(data = df, family = gaussian,
                     ciphers ~ 0 + experimental + control,
                     prior = c(prior(normal(10, 5), lb = 0, ub = 20, class = b),
                               prior(cauchy(0, 2.5), class = sigma)),
                     iter = 10000, warmup = 2000, chains = 4, sample_prior = "only",
                     seed = 19)

# summarizing model
summary(sample.priors)

# sample priors and make a difference score variable
df.priors.only <- as.data.frame(sample.priors)
df.priors.only$diff <- df.priors.only$b_experimental - df.priors.only$b_control

# checking that difference scores are reasonably centered around zero
round(quantile(df.priors.only$diff, probs=c(.05, .25, .5, .75, .95)), 2)
rethinking::dens(df.priors.only$diff, show.HPDI=.50)

# checking prior beliefs about sigma
ggplot(data=df.priors.only, aes(x=sigma)) +
  geom_histogram() +
  scale_x_continuous(limits=c(0,20))

# building a model with weakly informative priors
model.priors <- brm(data = df, family = gaussian,
                    ciphers ~ 0 + experimental + control,
                    prior = c(prior(normal(10, 5), lb = 0, ub = 20, class = b),
                              prior(cauchy(0, 2.5), class = sigma)),
                    iter = 10000, warmup = 2000, chains = 4,
                    seed = 19)

# summarizing the model
summary(model.priors)

# visualizing the model parameters and checking the chains
plot(model.priors)

# building a new data.frame that includes samples from the posterior
# and adding a difference variable
posterior <- as.data.frame(model.priors) %>%
  mutate(diff = b_experimental - b_control)

# calculating the effect sizes and adding to the new data.frame
posterior <- posterior %>% 
  mutate(cohen.d = diff / sigma,
         cles = pnorm(0, diff, sigma, lower.tail=FALSE))

# numerical summary of parameter estimates 
numerical.summary <- posterior %>%
  select(-lp__) %>% # dropping the lp__ variable
  select_all(funs(gsub("_", ".", .))) %>%
  summarise_all(list(mean = mean, 
                     sd = sd,
                     CI.95.lower = ~ rethinking::HPDI(samples=., prob=.95)[1],
                     CI.95.upper = ~ rethinking::HPDI(samples=., prob=.95)[2])) %>%
  t %>%
  as.data.frame %>%
  tibble::rownames_to_column() %>%
  separate(rowname, into = c("parameter", "fun"), sep="_") %>%
  pivot_wider(names_from = fun, values_from = V1)

# probability that the experimental effect is greater than 0
round(mean(posterior$diff > 0), 3)

# visual summaries of posterior distribution

# visualizing prior and posterior distributions for each condition
df.plot.prior2posterior <- posterior %>%
  select(b_experimental, b_control) %>%
  mutate(prior_experimental = df.priors.only$b_experimental,
         prior_control = df.priors.only$b_control)

prior2posterior.1 <- ggplot(data=df.plot.prior2posterior) +
  geom_density(aes(x=b_experimental), fill = "gray60") +
  geom_density(aes(x=prior_experimental), linetype="dashed") +
  scale_y_continuous(NULL, breaks = NULL) +
  ggpubr::theme_pubr() +
  theme(axis.title.x=element_blank())

prior2posterior.2 <- ggplot(data=df.plot.prior2posterior) +
  geom_density(aes(x=b_control)) +
  geom_density(aes(x=prior_control), linetype="dashed") +
  scale_y_continuous(NULL, breaks = NULL) +
  ggpubr::theme_pubr() +
  theme(axis.title.x=element_blank())

figure <- ggpubr::ggarrange(prior2posterior.1, prior2posterior.2, 
                            ncol = 1, nrow = 2)

annotate_figure(figure,
                bottom = text_grob("Ciphers completed", color = "black",
                                   hjust = .5, x = .5, size = 12),
                left = text_grob("Probability",
                                 color = "black", rot = 90))

# visualizing prior and posterior distributions for the experimental effect
prior2posterior.3 <- ggplot(data=posterior) +
  geom_density(aes(x=diff), fill = "gray60") +
  geom_density(data=df.priors.only, aes(x=diff), linetype="dashed") +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous(limits=c(-10,10), breaks=c(seq(-10,10,2))) +
  theme_pubr() +
  theme(axis.title.x=element_blank())

annotate_figure(prior2posterior.3,
                top = text_grob("Experimental Effect",
                                color = "black", hjust = .5, size = 12),
                bottom = text_grob(expression(alpha["experimental"] - alpha["control"]), color = "black",
                                   hjust = .5, x = .5, size = 12),
                left = text_grob("Probability",
                                 color = "black", rot = 90))

#
## Comparing alternative approaches
#

# building a model with uninformative priors
model.uninformative.priors <- brm(ciphers ~ 0 + experimental + control, data=df, 
                                  iter = 10000, warmup = 2000, 
                                  chains = 4, seed = 19)
summary(model.uninformative.priors)

# building a model that does not assume equal variances
uneq_var.model <- brm(ciphers ~ condition, sigma ~ condition, 
                      data=df,
                      iter = 46000, 
                      warmup = 45000, 
                      chains = 4,
                      seed = 4)
summary(uneq_var.model)

# independent samples t-test and Cohen's d effect size
t.test(ciphers ~ condition, data=df, var.equal=TRUE, alternative="less")
effsizecohen.d(ciphers ~ condition, data=df)