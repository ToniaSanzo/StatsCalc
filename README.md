# StatsCalc
A calculator to help with the statistics of the sciences course.

> Wald Interval:
>> p = "Successes / Total Trials"
>> z = critical z-score for the confidence level
>> n = sample size
>>>waldInterval(p, z, n)

>Wilson Interval:
>> p = "Successes / Total Trials"
>> z = critical z-score for the cofidence level
>> n = sample size
>>>wilsonInterval(p, z, n)

> Normal CDF
>> a = point a of the normally distributed graph
>> b = point b of the normally distributed graph
>> mu = mean of the normally distributed graph
>> sigma = standard deviation of the normally distributed graph
>>>normalCDF(a, b, mu, sigma)

> Inverse Normal
>> p = proportion of the graph
>> mu = mean of the normally distributed graph
>> sigma = standard deviation of th enormally distributed graph
>>>invNorm(p, mu, sigma)

>Prediction Interval
>> cl = Confidence level 0 - 1
>> sigma = standard deviation of th enormally distributed graph
>> n = sample size
>> xbar = standard deviation of the normally distributed graph
>>> predictionInterval(cl, sigma, n, xbar)

>One Sample T-Interval
>> cl = Confidence level 0 - 1
>> sigma = standard deviation of graph
>> n = sample size
>> xbar = mean of graph
>>>oneSampleT_Interval(cl, sigma, n, xbar)

>Unpaired T-Interval
>> cl = Confidence level 0 - 1
>> sigma1 = standard deviation of graph 1
>> n1 = sample size of graph 1
>> xbar1 = mean of graph 1
>> sigma2 = standard deviation of graph 2
>> n2 = sample size of graph 2
>> xbar2 = mean of graph 2

>Chi Squared Test
>> e = list of expected values
>> o = list of observed values
>> df = Degrees of freedom
>>>chiSquared_pVal(e, o, df)

>One-Tailed Hypothesis Test
>> xbar = sample mean
>> mu = mean
>> sigma = standard deviation
>> n = sample size
>>>oneSample_pVal_1Tail(xbar, mu, sigma, n)

>Two-Tailed Hypothesis Test
>> xbar = sample mean
>> mu = mean
>> sigma = standard deviation
>> n = sample size
>>>oneSample_pVal_2Tail(xbar, mu, sigma, n)

> Hypothesis Test Unpaired
>> xbar1 = sample mean 1
>> mu1 = mean 1
>> sigma1 = standard deviation 1
>> n1 = sample size 1
>> xbar2 = sample mean 2
>> mu2 = mean 2
>> sigma2 = standard deviation 2
>> n2 = sample size 2
>>>pVal_Unpaired(xbar1, mu1, sigma1, n1, xbar2, mu2, sigma2, n2)

> One-Way Anova
>> data = List of Lists for each populations data
>>> anova_oneWay(data)
