library(boot)
library(ggplot2)

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
num_species <- 10
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

ggplot(mapping=aes(x=x, y=p_interaction)) +
  geom_point() + geom_line() +
  geom_point(mapping=aes(y=interaction), col="red")


## add observation error
p_1_1 <- 0.9
p_0_0 <- 0.9
rans <- runif(num_species)

obs_interaction <- ifelse(interaction==1 & rans>p_1_1, interaction-1,
                          ifelse(interaction==0 & rans>p_0_0, interaction+1, interaction))


noise <- rbinom(1,1,y)
y <- y + noise

x_y <- data.frame(x=x,y=y)
# ~Interactions observed by a human
prob <- 0.3
n_sample <- as.integer(prob*num)
obs <- x_y[sample(nrow(x_y), n_sample), ]

#Parameter Estimation (ABC)
numm <- 10000
a <- runif(numm,0,5)
b <- runif(numm,-10,10)
Q = rep(NA,numm)


for(i in 1:numm)
{
  Q[i] <-  sum(abs(obs[,2]-a[i]*obs[,1]-b[i])) 
}
datta <- data.frame(a=a,b=b)
p_datta <- datta[Q<100,]

meann = colMeans(p_datta)

# ~Estimated Parameters
obs_a <- meann[1]
obs_b <- meann[2]
plot(x,y)
lines(x,obs_a*x+obs_b,type='l')
