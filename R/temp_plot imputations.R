# library(gamlss.cens)
# gen.cens()
#
# x = with(lung, Surv(time, status))
# ggplot(data = data.frame(x = with(lung, Surv(time, status))), aes(x =x))+
#   stat_function(fun = dNOrc, args = list(mu = mean(x), sigma = sd(x)))
# dNOrc(x, mu = mean(x), sigma = sd(x))
# mean(x)

# (Visualize draws) ------------------------------------------------------------------------------
require(ggplot2)
require(reshape2)

df = r$imputex$imputations
df$observation = seq(1, nrow(df)) # longformat index
df <- melt(df ,  id.vars = 'observation', variable.name = 'proposalVec')

ggplot(df, aes(observation,value)) +
  geom_boxplot(aes(group=observation)) +
  ylab('Proposals for observation [i]')

ggplot(df, aes(observation,value)) +
  geom_point() +
  ylab('Proposals for observation [i]')
