BHMMSimulation.R is used to run simulation study for Bayesian hierarchical mediation model (BHMM). One can specify the assumption of independent mediators or correlated mediators, and use different Gibbs sampling algorithm to generate posterior distributions of parameters.  If Gibbs_Ind.R is used, one assumes that mediators are independent. If Gibbs_Cor.R is used, one assumes there is correlation among mediators. Rcpp.cpp produces the same results as Gibbs_Cor.R, but the computation cost is significantly lower than Gibbs_Cor.R. UMAsimulation.R is used to run simulation study for univariate mediation analysis(UMA) as comparison for BHMM.
