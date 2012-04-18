#GTDistrib
This unit contains definitions of objects that represent probability distributions.
 The object hierarchy is as follows:

         tDistribution
             tContinuous_distribution
                 tNormal_distribution
                 tBeta_distribution
                 tGamma_distribution
                     tChi_square_distribution
                     tExponential_distribution
                 tF_distribution
                 tT_distribution
                 tCauchy_distribution
                 tDouble_exponential_distribution
                 tGumbel_distribution
                 tLogistic_distribution
                 tLognormal_distribution
                 tPareto_distribution
                 tUniform_distribution
                 tWeibull_distribution
                 tNon_central_Chi_square_distribution
                 tNon_central_Beta_distribution
                 tNon_central_F_distribution
                 tNon_central_t_distribution
             tDiscrete_distribution
                 tBinomial_distribution
                     tBernoulli_distribution
                 tHypergeometric_distribution
                 tPoisson_distribution
                 tDiscrete_Uniform_distribution
                 tNegative_Binomial_distribution
                     tGeometric_distribution

## For the user the most interesting methods will be:
         . cdf(x) (cumulative distribution function)
         . quantile(q) (inverse cdf)
         . random_value (gives a random draw from the distribution)
         ### Further methods include:
         . density(x) (for continuous distributions)
         . probability(x) (for discrete distributions)
         . in_support(x) (true if x is in the support of the distribution)
         . expectation_exists and variance exists (true or false)
         . mean and variance

## In Resource file:
	id_error, "Caution"
	err_val_le, "Value should be greater than "
	err_val_l, "Value should be at least "
	err_val_ge, "Value should be less than "
	err_val_int, "Value should be integer."
	err_val2_le, "Second parameter should be greater than "
	err_val2_l, "Second parameter should be at least "
	err_val2_ge, "Second parameter should be smaller than "
	err_val2_int, "Second parameter should be integer."
	err_val3_le, "Third parameter should be greater than "
	err_val3_l, "Third parameter should be at least "
	err_val3_ge, "Third parameter should be less than "
	err_val3_int, "Third parameter should be integer."
	err_val1_l_val2, "The second parameter should not be greater than the first."
	err_val1_l_val3, "The third parameter should not be greater than the first."
	err_val_g, "Value should be not greater than "
	err_overflow, "Numbers too large."
	err_expectation, "Expectation does not exist."
	err_variance, "Variance and standard deviation do not exist."