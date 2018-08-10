# (Visualize draws) ------------------------------------------------------------------------------
require(ggplot2)
require(reshape2)

# make it a class method!
plotimputations <- function(object, boxes = TRUE, quantiles = FALSE) {
  df = object$imputations
  if(quantiles == TRUE){
    print(ggplot(object$impquantiles,
                 aes(x = seq(1:nrow(d2$impquantiles)),
                     ymin=q5,
                     lower=q25,
                     middle=q50,
                     upper=q75,
                     ymax=q95)) +
            geom_boxplot(stat="identity")+
            xlab('draw'))
  }


  # Convert to Longformat
  df$observation = seq(1, nrow(df))
  df <- melt(df ,  id.vars = 'observation', variable.name = 'proposalVec')

  if (boxes) {
    print( ggplot(data = subset(df, df$proposalVec != 'imputedx'), aes(observation, value)) +
             geom_boxplot(aes(group = observation)) +
             geom_point(data = subset(df, df$proposalVec == 'imputedx'), aes(observation, value, color = 'red'))+
             ylab('Proposals for observation [i]')) # NOTE that red dots are based on mean, boxes display median.
  }else {
    print( ggplot() +
             geom_point(data = subset(df, df$proposalVec != 'imputedx'), aes(observation, value)) +
             geom_point(data = subset(df, df$proposalVec == 'imputedx'), aes(observation, value, color = proposalVec))+
             ylab('Proposals for observation [i]'))
  }
}

# quantiles should be highly skewed
plotimputations(d2, boxes = FALSE, quantiles = TRUE)  # the imputed value is ommited


