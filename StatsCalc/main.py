# Author: Tonia Sanzo
#
# Program that calculates a upper and lower bound of a confidence interval using z-interval

import math
import scipy.stats


# Determine the critical t* value 
invT = scipy.stats.t.ppf


# Determine confidence interval using the p-value (p = x/n "success/total"),
#     z* score, and sample size(n)
# 
# req: x >= 5 && n - x >= 5 (# of successes is greater than or equal to 5
#     # of non-successes is greater than or equal to 5)
def waldInterval(p, z, n):
    lb = 0.0
    ub = 0.0
    lb = p - (z * (math.sqrt((p * (1 - p))/n)))
    ub = p + (z * (math.sqrt((p * (1 - p))/n)))
    print("LB: ", lb, "UB: ", ub) 


# Determine confidence interval using the p-value (p = x/n "success/total"),
#     z* score, and sample size(n)
#
# req: n > 2 (sample size is larger than 2)
def wilsonInterval(p, z, n):
    lb = 0.0
    ub = 0.0
    n1 = p + ((z * z)/ (2 * n))
    n2 = z * math.sqrt(((p * (1 - p)) / n) + ((z * z) / (4 * (n * n))))
    d  = 1 + ((z * z) / n)
    lb = (n1 - n2)/ d 
    ub = (n1 + n2)/ d
    print("LB:", lb, "UB:", ub)


# Returns the proportion of the interval between a - b for a normally 
#     distributed graph with a specified mu(mean) and sigma(standard deviation)
#
# req: normal distribution
def normalCDF(a, b, mu, sigma):
    cdf = scipy.stats.norm.cdf
    return cdf(b, mu, sigma) - cdf(a, mu, sigma)


# Given a proportion, the mean, and standard deviation of a normally
#     distributed graoh, determine the value associated with that
#     proportion
#
# req: normal distribution
def invNorm(p, mu, sigma):
    return scipy.stats.norm.ppf(p, mu, sigma)


# Determine the confidence interval for a normally distributed graph 
#     using the one sample t-interval techniques
#
# req: normal distribution
def predictionInterval(cl, sigma, n, xbar):
    lb = 0.0
    ub = 0.0
    t = abs(invT((1 - cl) / 2, n - 1))
    lb = xbar - ((t * sigma / math.sqrt(n)) * math.sqrt(n + 1))
    ub = xbar + ((t * sigma / math.sqrt(n)) * math.sqrt(n + 1))
    print("LB:", lb, "UB:", ub, "with", cl * 100, "% confidence")
    
    
# Determine the confidence interval for a normally distributed graph 
#     using the one sample t-interval techniques
#
# req: normal distribution or n > 30
def oneSampleT_Interval(cl, sigma, n, xbar):
    lb = 0.0
    ub = 0.0
    t = abs(invT((1 - cl) / 2, n - 1))
    lb = xbar - ((t * sigma / math.sqrt(n)))
    ub = xbar + ((t * sigma / math.sqrt(n)))
    print("LB:", lb, "UB:", ub, "with", cl * 100, "% confidence")


# Used in t-interval for mu1 - mu2 (unpaired data)
def determineAB(sigma, n):
    return (sigma * sigma)/n


# Determine degree's of freedom for the t-Interval of tInterval_Unpaired
def satterthwaiteFormula(a, b, n1, n2):
    return math.floor((a + b) * (a + b)) / (((a * a) / (n1 - 1)) + ((b * b) / (n2 - 1)))


# Return the interval of mu1 - mu2 for 2 normally distributed unpaired data
#
# req: x1 and x2 have to be normally distributed
def tInterval_Unpaired(cl, sigma1, n1, xbar1, sigma2, n2, xbar2):
    a = sigma1**2 / n1
    b = sigma2**2 / n2
    df = satterthwaiteFormula(a, b, n1, n2)
    print("df:",df)
    t = abs(invT((1 - cl) / 2, df))
    print("t:",t)
    comp1 = xbar1 - xbar2
    comp2 = t * math.sqrt(a + b)
    lb = comp1 - comp2
    ub = comp1 + comp2
    print(lb, "<= mu1 - mu2 <=", ub)


# Determine the test statistic and p-val for a chi squared test
#
# req: expected frequency of each category is > 5
def chiSquared_pVal(e, o, df):
    x0 = 0.0
    i = 0

    while i < len(e):
        x0 += ((e[i] - o[i])**2)/e[i]
        i += 1

    pval = 1 - scipy.stats.chi2.cdf(x0, df)
    print("Test Statistic x0:", x0, "P-Value:", pval)
    

# Determine test statistic and p-value for hypothesis testing for a one-tailed
#     test
#
# req: normal distribution or n > 30
def oneSample_pVal_1Tail(xbar, mu, sigma, n):
	t = ((xbar - mu)/ sigma) * math.sqrt(n)
	pval = scipy.stats.t.sf(abs(t), n-1)
	print("t-statistic =", t, "p-val=",pval)
	
	
# Determine test statistic and p-value for hypothesis testing for a two-tailed
#     test
#
# req: normal distribution or n > 30
def oneSample_pVal_2Tail(xbar, mu, sigma, n):
        t = ((xbar - mu)/ sigma) * math.sqrt(n)
        pval = 2 * scipy.stats.t.sf(abs(t), n-1)
        print("t-statistic =", t, "p-val=",pval)


# Determine test statistic and p-value for hypothesis testing for a 2-tailed 
# test
#
# req: both samples are normal distribution or both samples sizes are > 30
def pVal_Unpaired(xbar1, mu1, sigma1, n1, xbar2, mu2, sigma2, n2):
	a = sigma1**2/n1
	b = sigma2**2/n2
	t = ((xbar1 - xbar2) - (mu1 - mu2))/math.sqrt(a - b)
	df = satterthwaiteFormula(a, b, n1, n2)
	pval = 2 * scipy.stats.t.sf(abs(t), n-1)
	print("t-statistic =", t, "p-val=",pval)
	
	
# Test hypotheses that two or more groups have the same population mean.
#
# req: 1) Each population is normally distributed (If samples are  too small look
#         at normal probablity plot of combined within-sample deviations)
#      2) Variances are equal across populations such that: 
#                    (sigma^2 max) / (sigma^2 min) <= 4
#                     [if samples are small <= 6 is OK]
def anova_oneWay(data):
    n = 0          # Total number of observations in data
    sse = 0        # Sum of Squares (error)
    temp_mean = 0  # Used to temporarly store a mean
    
    # Determine the number of groups
    k = len(data)
    for group in data:
        # Determine the total number of observations
        length = len(group)
        n += length
        i = 0
        
        # Determine the mean for a single group
        while i < length:
            temp_mean += group[i]
            # If the loop is in it's final iteration determine the mean, then 
            # calculate SSE and prep data for SSF
            if (i + 1) == length:
                temp_mean = temp_mean / length
                j = 0
                while j < length:
                    sse += (group[j] - temp_mean) ** 2
                    group[j] = temp_mean
                    j += 1
                temp_mean = 0
            i += 1
    print('n', n)
    print('k', k)
    # Determine degree's of freedom
    df_factor = k - 1
    df_error = n - k
    df_total = df_factor + df_error
    
    # Determine grand mean
    for group in data:
        for i in group:
            temp_mean += i
    temp_mean = temp_mean / n
    print('grand mean', temp_mean)

    sstr = 0
    # determine SSF
    for group in data:
        for i in group:
            sstr += (i - temp_mean) ** 2
    print('df[factor]', df_factor)
    print('df[error]', df_error)
    print('df[total]', df_total)
    print('SS[factor]', sstr)
    print('SS[error]', sse)
    
    ms_factor = sstr/ df_factor
    ms_error  = sse/ df_error
    f_val     = ms_factor / ms_error
    p_val     = 1 - scipy.stats.f.cdf(f_val, df_factor, df_error)
    print('ms[factor]', ms_factor)
    print('ms[error]', ms_error)
    print('f-value', f_val)
    print('p-value', p_val)


if __name__ == "__main__":
