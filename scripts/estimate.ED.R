#
# estimate.ED
#
# Marco Del Giudice (2020). Version 1. Contact: marcodg@unm.edu
#
# A simple R function to compute estimates of effective dimensionality (ED), either from raw data or from a correlation/covariance matrix. The function returns four estimators, both uncorrected and corrected for small-sample bias if enough information is available. Disattenuated estimates can also be computed if desired, by supplying a vector of reliability coefficients. For more information see Del Giudice (2020).
#
# ED estimators: The n1 index is suitable as a general-purpose estimator of ED (Cangelosi & Goriely 2007, Roy & Vetterli 2007, Gnedenko & Yelnik 2016); n2 is more conservative and can be used as a lower bound estimate of dimensionality (Fraedrich et al. 1995, Bretherton et al. 1999, Pirkl et al. 2012). The estimators nC (Cheverud 2001, Wagner et al. 2008) and nInf (Kirkpatrick et al. 2009) are not recommended, and are only reported for comparison with the published literature. Whereas nC typically overestimates ED, nInf is extremely conservative and insensitive to changes in the correlation structure.  
#
# Small-sample correction: if raw data are supplied, the function uses the nlshrink package (Ramprasad 2016) to perform nonlinear shrinkage of the eigenvalues (Ledoit & Wolf 2012, 2015). If only the correlation or covariance matrix is supplied, the function applies the correction method by Mestre (2008). If the corrected ED is lower than the uncorrected ED (which may happen when the uncorrected ED is very close to the number of variables), the uncorrected value is returned in place of the corrected one. IMPORTANT: implementing Mestre's method is tricky; sometimes, the uniroot.all function may fail to find the correct mu values. When applying the small-sample correction to a correlation/covariance matrix, please check the results for unexpected or implausible values.
#
# Disattenuation: to correct for measurement error, disattenuated estimates of ED can be obtained by supplying a vector of reliability coefficients (e.g., Cronbach’s α, McDonald's omega-t or omega-h). For details see Del Giudice (2020). 
#
# Note: If the matrix supplied by the user is indefinite (or becomes indefinite after disattenuation), Higham's (2002) method is automatically applied to replace it with the nearest positive definite matrix (Matrix package by Bates & Maechler, 2019). A warning message is displayed in the output summary.
#
#
# References:
# Bates, D., & Maechler, M. (2019). Matrix v. 1.2-17. URL: https://CRAN.R-project.org/package=Matrix 
# Bretherton, C. S., Widmann, M., Dymnikov, V. P., Wallace, J. M., & Bladé, I. (1999). The effective number of spatial degrees of freedom of a time-varying field. Journal of Climate, 12, 1990-2009.
# Cangelosi, R., & Goriely, A. (2007). Component retention in principal component analysis with application to cDNA microarray data. Biology Direct, 2, 2.
# Cheverud, J. M. (2001). A simple correction for multiple comparisons in interval mapping genome scans. Heredity, 87, 52-58.
# Del Giudice, M. (2020). Effective dimensionality: A tutorial. Multivariate Behavioral Research. https://doi.org/10.1080/00273171.2020.1743631
# Fraedrich, K., Ziehmann, C., & Sielmann, F. (1995). Estimates of spatial degrees of freedom. Journal of Climate, 8, 361-369.
# Gnedenko, B., & Yelnik, I. (2016). Minimum entropy as a measure of effective dimensionality. SSRN, doi: 10.2139/ssrn.2767549.
# Higham, N. J. (2002). Computing the nearest correlation matrix—a problem from finance. IMA journal of Numerical Analysis, 22, 329-343.
# Kirkpatrick, M. (2009). Patterns of quantitative genetic variation in multiple dimensions. Genetica, 136, 271-284.
# Ledoit, O., & Wolf, M. (2012). Nonlinear shrinkage estimation of large-dimensional covariance matrices. Annals of Statistics, 40, 1024-1060.
# Ledoit, O., & Wolf, M. (2015). Spectrum estimation: A unified framework for covariance matrix estimation and PCA in large dimensions. Journal of Multivariate Analysis, 139, 360-384.
# Pirkl, R. J., Remley, K. A., & Patané, C. S. L. (2012). Reverberation chamber measurement correlation. IEEE Transactions on Electromagnetic Compatibility, 54, 533-545.
# Ramprasad, P. (2016). nlshrink v. 1.0.1. URL: https://CRAN.R-project.org/package=nlshrink
# Roy, O., & Vetterli, M. (2007). The effective rank: A measure of effective dimensionality. In M. Domański, R. Stasiński, & M. Bartkowiak (Eds.), 15th European Signal Processing Conference (pp. 606-610). Poznań, Poland: PTETiS.



# function estimate.ED (x, sample.size = NULL, rel.values = NULL, cov.mat = FALSE, small.sample.c = FALSE, round.digits = 2, print.summary = TRUE)
#
# Arguments
#
# x					data frame or correlation/covariance matrix
# sample.size		sample size (optional: required for small-sample bias correction from a correlation/covariance matrix)
# rel.values		vector of reliability coefficients (optional: only for disattenuation)
# cov.mat			if TRUE, ED estimators are computed from the covariance matrix; if FALSE (default), from the correlation matrix
# small.sample.c	if TRUE (default), the function returns estimators corrected for small-sample bias
# round.digits		rounding digits for the output (default is 2)
# print.summary		if TRUE (default), the function prints a summary of the analysis before returning the results
#
# Value
# returns a list object containing some or all of the following:
#
# n1				ED estimator based on the Shannon entropy H1 
# n2				ED estimator based on the quadratic entropy H2 
# nInf				ED estimator based on the min-entropy H_Inf; extremely conservative, not recommended
# nC				ED estimator by Cheverud (2001); not recommended
# n1.c				n1 corrected for small-sample bias
# n2.c				n2 corrected for small-sample bias
# nInf.c			nInf corrected for small-sample bias
# nC.c				nC corrected for small-sample bias



library(rootSolve)
# library(nlshrink)
library(Matrix)

# cov_mat <- Rfast::cova(train_dat[, useMarkers] %>% as.matrix())

estimate.ED <- function (x, sample.size = NULL, rel.values = NULL, cov.mat = TRUE, small.sample.c = TRUE, round.digits = 2, print.summary = TRUE) {

indefinite.matrix = FALSE

########## check the input format, compute the correlation or covariance matrix

if (is.data.frame(x) == TRUE) {
	sample.size = nrow(x)
	
	if (cov.mat == FALSE) {
		matrix = cor(x, use="pairwise.complete.obs")
	}
	else matrix = cov(x, use="pairwise.complete.obs")
}
else {
	matrix = x

	if ( ((sum(diag(matrix)) != nrow(matrix) ) | (var(diag(matrix)) != 0)) & (cov.mat == FALSE) ) {  # checks if the matrix is a covariance matrix
		matrix = cov2cor(matrix)
	}	
}

########## disattenuate the correlation matrix if a vector of reliabilities has been supplied

if ( (cov.mat == FALSE) & (is.null(rel.values) == FALSE) ) {
	rel.matrix = sqrt(crossprod(t(rel.values), t(rel.values)))
	diag(rel.matrix) = 1
	matrix = matrix/rel.matrix
}
	
output = list()

########## check if the matrix is indefinite; if so, find the nearest positive definite matrix

if (sum(eigen(matrix)$values < 0) > 0) {
	indefinite.matrix = TRUE
	matrix = nearPD(matrix, corr=!cov.mat)$mat
}

########## find the eigenvalues, sort them in ASCENDING order for later computations

eigen_val = sort(eigen(matrix)$values)

eigen_val.c = NULL

########## correct the eigenvalues for small-sample bias 

if (small.sample.c == TRUE) {

	if (is.data.frame(x) == TRUE) {   #  from raw data: nonlinear shrinkage by Ledoit & Wolf 2012, 2015
			
		if(cov.mat == FALSE) {
			
			if(is.null(rel.values) == TRUE) {
				suppress = file()
				sink(file = suppress)  #  suppresses the print output from tau_estimate()
				tmp = tau_estimate(as.matrix(scale(x)))
				sink()
				close(suppress)
				eigen_val.c = tmp
			}
			else {  #  disattenuation
				sink("NUL")  #  suppresses the print output from tau_estimate()
				tmp = nlshrink_cov(as.matrix(scale(x)))
				sink()
				matrix.c = tmp/rel.matrix  #   disattenuate the corrected correlation matrix
				if (sum(eigen(matrix.c)$values < 0) > 0) {
					indefinite.matrix = TRUE
					matrix.c = nearPD(matrix.c, corr=TRUE)$mat  #  check if the matrix is indefinite, find the nearest positive definite matrix
				}
				eigen_val.c = sort(eigen(matrix.c)$values)
			}
		}	
		else {
			suppress = file()
			sink(file = suppress)  #  suppresses the print output from tau_estimate()
			tmp = tau_estimate(as.matrix(x))
			sink()
			close(suppress)
			eigen_val.c = tmp			
		}
		
	}
	else if (is.null(sample.size) == FALSE) {   #  from correlation/covariance matrix: correction by Mestre 2008
				
		eigen_val.rounded = round(eigen_val, 6)  #  to avoid errors when detecting duplicates
		eigen_val.c = eigen_val
			
		function.body = paste(eigen_val.rounded[1],"/(",eigen_val.rounded[1],"-x)")  #  build equation to find mu values
		for(index in 2:length(eigen_val.rounded)) {
			function.body = paste(function.body, "+",eigen_val.rounded[index],"/(",eigen_val.rounded[index],"-x)" )
		}
		function.body = paste(function.body,"-",sample.size)
		
		eval(parse(text = paste("f = function(x) {", function.body,"}")))
		
		unique.eigen = sum(!duplicated(eigen_val.rounded[eigen_val.rounded != 0]) )  #   number of unique nonzero eigenvalues
		
		mu_val = uniroot.all(f, interval=c(0, max(eigen_val.rounded)), n=10000000)  #  find mu values
		mu_val = mu_val[which( abs(f(mu_val)) %in% sort(abs(f(mu_val)))[1:unique.eigen]  ) ]   
		mu.index = 1
			for (index in 1:length(eigen_val.rounded)) {   
				
				if (eigen_val.rounded[index] != 0) {   # correct nonzero eigenvalues
					if (duplicated(eigen_val.rounded)[index] == FALSE) {
						eigen_val.c[index] = sample.size * (eigen_val.rounded[index] - mu_val[mu.index])
						mu.index = mu.index + 1
					}
					else {
						eigen_val.c[index] = eigen_val.c[index-1]
						}
				}
			}
	}

}

########## compute ED estimators

K = length(eigen_val)	
eigen_sum = sum(eigen_val)
norm_eigen_val = eigen_val/eigen_sum
eigen_var = var(eigen_val)*((K-1)/K)

 
output$n1 = prod(norm_eigen_val^(-norm_eigen_val))  
output$n2 = (eigen_sum^2)/sum(eigen_val^2) 
output$nInf = eigen_sum/max(eigen_val)    
output$nC = K - ((K^2)/(eigen_sum^2))*eigen_var  

if ( (small.sample.c == TRUE) & (is.null(eigen_val.c) == FALSE) ) {  # estimators corrected for small-sample bias
	eigen_sum.c = sum(eigen_val.c)
	norm_eigen_val.c = eigen_val.c/eigen_sum.c
	eigen_var.c = var(eigen_val.c)*((K-1)/K)

 	output$n1.c = max(prod(norm_eigen_val.c^(-norm_eigen_val.c)), output$n1)
	output$n2.c = max((eigen_sum.c^2)/sum(eigen_val.c^2), output$n2)
	output$nInf.c = max(eigen_sum.c/max(eigen_val.c), output$nInf) 
	output$nC.c = max(K - ((K^2)/(eigen_sum.c^2))*eigen_var.c, output$nC) 
}

########## print summary of the analysis if requested

if (print.summary == TRUE) {

	description = "ED estimated from the"
	if (cov.mat == TRUE) {description = paste(description, "covariance matrix;")}
	else {description = paste(description, "correlation matrix;")}
	if (is.null(rel.values) == FALSE) {description = paste(description, "disattenuated;")}
	else {description = paste(description, "no disattenuation;")}
	if ( (small.sample.c == TRUE) & (is.null(eigen_val.c) == FALSE) ) {
		if (is.data.frame(x) == TRUE) {description = paste(description, "corrected for small-sample bias (Ledoit & Wolf 2015).")}
		else {description = paste(description, "corrected for small-sample bias (Mestre 2008).")}
	}
	else {description = paste(description, "no small-sample correction.")}
	if (indefinite.matrix == TRUE) {description = paste(description, "Warning: an indefinite matrix was detected and replaced with the nearest positive definite matrix.")} 
	
	print(description, quote = FALSE)
	cat("\n")
}

########## return output list

return(lapply(output,round,round.digits))
}


# res <- estimate.ED(x = cov_mat, cov.mat = TRUE, round.digits = 0)
