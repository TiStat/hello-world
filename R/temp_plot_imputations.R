#'@description Given an imputex object, this funciton plots the imputations and
#'  optionally plots the approximate averaged quantiles of the valid part of the
#'  censored conditonal distributions, from which the proposed vectors were
#'  'drawn', before the proposals were aggregated to become the imputed vector.
#'@param object imputex object.
#'@example plotimputations(d, boxes = FALSE, quantiles = TRUE)
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




