#  Evaluating and Testing for Actionable Treatment Effect Heterogeneity

Developing tools for estimating heterogeneous treatment effects and individualized treatment effects (HTE) has been an area  of active research over the past decade. While these tools have proven to be useful in many contexts, a concern when deploying such methods is the degree to which incorporating HTE into a prediction model provides an advantage over predictive methods which do not allow for variation in treatment effect across individuals. To address this concern, we propose  a procedure which evaluates the extent to which an HTE model provides actionable benefit. Specifically, our procedure targets the gain in predictive performance  from using a flexible predictive model incorporating HTE versus an alternative model which is similar to the full, HTE-utilizing model except that it is constrained to not allow variation in treatment effect. By drawing upon recent work in using nested cross-validation techniques for prediction error inference, we generate confidence intervals for this gain in predictive performance  which allows one to directly calculate the level at which one is confident  of a substantial HTE-modeling gain in prediction -- a quantity which we refer to as the h-value. Our procedure is generic and can be directly used to assess the HTE-modeling benefit of any method that incorporates treatment effect variation.

# Dataste
The raw dataset can be downloaded from [here](http://www.icpsr.umich.edu/icpsrweb/HMCA/studies/9795?paging.startRow=51). The dataset is titled: 'DS141: Transport Format SAS
Library Containing the 59 Evaluation'.The cleaned dataset containing the variables used in this study is available at [here](https://github.com/mahsaashouri/HTE-Model-Comparison).

# Reproducing results
* `Sim1Run.R`: results for linear conditional mean function simulation study (*Tables 2*, and *B1*)
* `SimTwoRunParallel.R`: results for non-linear conditional mean functions simulation study (*Tables 3*, *4*, and *B2*)
* `plot_example_compare.R`: simulated data illustration (*Figure B1*)
* `IHDP-Dataset.R`: histogram of h-values for IHDP dataset (*Figures B2*, and *B4*)
* `Par_Dep_Plot.R`: partial dependence plots using glmboost (*Figure B3*)
  
