library(boot)
library(ggplot2)
library(dplyr)
library(tidyr)

#logistic model y = a*x + b  ~Food Web Model
logistic <- function(x, a, b){
  yy <- inv.logit(a*x+b)
  return(yy)
}

x <- seq(-10, 10, length=100)
ggplot(mapping=aes(x=x, y=logistic(x, 1, 1))) +
  geom_line()


## parameters
mass_mean <- 1
mass_sd <- 1
num_species <- 10000
real_a <- 1 # ~Real Parameters
real_b <- -5 # ~Real Parameters



x <- rlnorm(num_species, mass_mean, mass_sd) # ~Real Body Sizes
# ~Real food web
p_interaction <- logistic(x, real_a,real_b) 

# ~Stochasticity in the system's biology
interaction <- numeric(length=num_species)
for(i in 1:num_species)
  interaction[i] <- rbinom(1, 1, p_interaction[i])

## deterministic biology version
#interaction <- ifelse(p_interaction<0.3, 0, 1)


## add multiple observations per interaction (perhaps body size dependent)


ggplot(mapping=aes(x=x, y=p_interaction)) +
  geom_point() + geom_line() +
  geom_point(mapping=aes(y=interaction), col="red")


## add observation error
p_1_1 <- 1  ## larger values give lower observation error
p_0_0 <- 1  ## larger values give lower observation error
rans <- runif(num_species)

obs_interaction <- ifelse(interaction==1 & rans>p_1_1, interaction-1,
                          ifelse(interaction==0 & rans>p_0_0, interaction+1, interaction))

ggplot(mapping=aes(x=x, y=p_interaction)) +
  geom_point() + geom_line() +
  geom_point(mapping=aes(y=interaction), col="red", size=3) +
  geom_point(mapping=aes(y=obs_interaction), col="green")


#Parameter Estimation (ABC)
numm <- 10000
guesses <- data.frame(a=runif(numm,-5,5),
                      b=runif(numm,-10,10),
                      Q = rep(NA,numm))

for(i in 1:numm)
{
  guesses$Q[i] <-  -sum(dbinom(obs_interaction, 1, logistic(x, a[i], b[i]), log = TRUE))
}

ggplot(gather(guesses), aes(x=value)) +
  facet_wrap(~key, scales="free_x") +
  geom_histogram()

guesses %>%
  filter(Q<150) %>%
  gather() %>%
  ggplot(aes(x=value)) +
     facet_wrap(~key, scales="free_x") +
     geom_histogram()


means <- guesses %>%
  filter(Q<50) %>%
  summarise_all(mean)
  


ggplot(mapping=aes(x=x, y=p_interaction)) +
  geom_point() + geom_line() +
  geom_point(mapping=aes(y=interaction), col="red", size=3) +
  geom_point(mapping=aes(y=obs_interaction), col="green") +
  geom_line(aes(y=logistic(x, means$a, means$b)))


m1 <- glm(obs_interaction ~ x, family=binomial)
summary(m1)
