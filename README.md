# Rcode
This DP model aims at handling simple linear mixture regression model(i.e. design matrix is composed of an intercept and a predictor). 
An example is included in the Debugging section at the bottom of myDP_ori.R file. 
De-annotate Debugging section to get a quick start of myDP code.<br/>
dpreassign.R is a corrected version of the original dirichletpackage, which has revised the updating function for Gaussian kernel from conjugate to a semi-conjugate way. It can detect latent gaussian distributions in the dataset. Some other functions were also refined to boost training. See "A First Course in Bayesian Statistical Methods" Chapter 6.1 for detailed process. Author: Peter D. Hoff

