"""
Operative file


"""


import pandas as pd
# import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sb
# from lifelines.plotting import plot_lifetimes
from lifelines import CoxPHFitter
# import sklearn
# from sklearn import preprocessing, linear_model
from lifelines import KaplanMeierFitter
# from lifelines import WeibullAFTFitter
# from pandas import DataFrame, plotting
# from numpy.random import uniform, exponential
# from lifelines.utils import median_survival_times
from lifelines.statistics import multivariate_logrank_test


# importing data and calling some functions
df = pd.read_csv("aug_3_new_genes.csv", sep=",")

kmf = KaplanMeierFitter()
cph = CoxPHFitter()


# doesnt work for whatever reason
# print("kmf median surv t: ")
# print(kmf.median_survival_time_)


# kaplan meier fitter, should give a graph. works!
T = df["duration"]
E = df["event"]


# This works, gives tabular view of coefficients and related stats
cph.fit(df, duration_col="duration", event_col="event")
print("cph summary: ")
cph.print_summary()


# log rank test tells you siginificance between survival curves
# ix = df["KRAS"] == "1"
# T_exp, E_exp = df.loc[ix, T], df.loc[ix, E]
# T_con, E_con = df.loc[~ix, T], df.loc[~ix, E]

results = multivariate_logrank_test(T, df["tp53"], E)
print("log rank: ")
results.print_summary()


# kaplan meier on differences within individual columns
groups = df["tp53"]
ix0 = (groups == 0)
ix1 = (groups == 1)

kmf.fit(T[ix0], E[ix0], label="No")
ax0 = kmf.plot()

kmf.fit(T[ix1], E[ix1], label="Yes")
ax1 = kmf.plot()
plt.show()


# cumulitiv density
# kmf.fit(T, event_observed=E)
# kmf.plot_cumulative_density()
# plt.show()


# box whisker plot shows differences between cohorts
cph.fit(df, duration_col="duration", event_col="event")
cph.plot()
plt.show()


# kaplan meier survival fx
kmf.fit(T, event_observed=E)
kmf.plot_survival_function()
plt.show()