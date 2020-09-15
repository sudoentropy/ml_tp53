"""
Operative code for tp53 paper


Kristin Syed BS. University of Toledo College of Medicine (UTCOM).
3000 Arlington Ave.
Toledo, OH 43614
419.383.4000
Vyshnavi Reddy MS. UTCOM.
Judy Daboul BS. UTCOM.
Nicholas Thompson BS.
Hamza Syed BS. Independent Data Scientist.
Gang Ren [credentials and affiliations]
Joseph Sferra MD. Promedica Toledo Hospital.
Jeffrey Sutton MD. UTCOM.
"""


import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test


# importing data and calling functions
df = pd.read_csv("all_genes_60months.csv", sep=",")

kmf = KaplanMeierFitter()
cph = CoxPHFitter()


# kaplan meier fitter
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
