# NOT working example yet
# data = simulateData(method ='right' ,
#                     n = 1000,
#                     ymu.formula = ~ 0.5*x1,
#                     ysigma.formula ~ 0.3*x2,
#                     family = 'NO')
#
# d2 <- imputex(data = data$defected,
#               xmu_formula= x1~y+x2,
#               xsigma_formula = ~1,
#               xnu_formula = ~1,
#               xtau_formula = ~1,
#               xfamily = NO(mu.link = 'identity'),
#               indicator = "indicator",
#               censtype )

# (Visualize draws) ------------------------------------------------------------------------------
require(ggplot2)
require(reshape2)

# make it a class method!
plotimputations <- function(df, boxes = TRUE) {

  # Convert to Longformat
  df$observation = seq(1, nrow(df))
  df <- melt(df ,  id.vars = 'observation', variable.name = 'proposalVec')

  if (boxes) {
    return(ggplot(df, aes(observation, value)) +
             geom_boxplot(aes(group = observation)) +
             ylab('Proposals for observation [i]'))
  }else {
    return(ggplot(df, aes(observation, value)) +
             geom_point() +
             ylab('Proposals for observation [i]'))
  }
}

plotimputations(r$imputex$imputations[, -ncol(r$imputex$imputations)],
                boxes = FALSE)  # the imputed value is ommited

