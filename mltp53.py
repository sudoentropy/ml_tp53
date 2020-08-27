"""
Operative file for tp53 paper


"""


import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test


# importing data and calling functions
df = pd.read_csv("aug_3_new_genes.csv", sep=",")

kmf = KaplanMeierFitter()
cph = CoxPHFitter()


# kaplan meier fitter, should give a graph. works!
T = df["duration"]
E = df["event"]


# This gives tabular view of coefficients and related stats
cph.fit(df, duration_col="duration", event_col="event")
print("cph summary: ")
cph.print_summary()


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


# box whisker plot shows differences between cohorts
cph.fit(df, duration_col="duration", event_col="event")
cph.plot()
plt.show()


# kaplan meier survival fx
kmf.fit(T, event_observed=E)
kmf.plot_survival_function()
plt.show()
