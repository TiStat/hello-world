# (Visualize draws) ------------------------------------------------------------------------------
require(ggplot2)
require(reshape2)

# make it a class method!
plotimputations <- function(df, boxes = TRUE) {

  # Convert to Longformat
  df$observation = seq(1, nrow(df))
  df <- melt(df ,  id.vars = 'observation', variable.name = 'proposalVec')

  if (boxes) {
    return(ggplot(data = subset(df, df$proposalVec != 'imputedx'), aes(observation, value)) +
             geom_boxplot(aes(group = observation)) +
             geom_point(data = subset(df, df$proposalVec == 'imputedx'), aes(observation, value, color = 'red'))+
             ylab('Proposals for observation [i]')) # NOTE that red dots are based on mean, boxes display median.
  }else {
    return(ggplot() +
             geom_point(data = subset(df, df$proposalVec != 'imputedx'), aes(observation, value)) +
             geom_point(data = subset(df, df$proposalVec == 'imputedx'), aes(observation, value, color = proposalVec))+
             ylab('Proposals for observation [i]'))
  }
}

plotimputations(d2$imputations[, -ncol(r$imputex$imputations)],
                boxes = FALSE)  # the imputed value is ommited
