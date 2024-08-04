import os
import pandas as pd
import glike

results = []
for i in range(1, 51):
  for j in range(1, 21):
    if not os.path.isfile(f"output_{i}/{j}.bestlhoods"):
      continue
    results_ = pd.read_csv(f"output_{i}/{j}.bestlhoods", delimiter = "\t")
    results_.insert(0, "group", i)
    results.append(results_)

results = pd.concat(results)
columns = ['group', 
           'T1$', 'T2$', 'T3$', 'T4$', 
           'R1$', 'R2$',
           'Nafr$', 'Neur$', 'Nasia$', 'Nadmix$', 'Nooa$', 'Nanc$',
           'GReur$', 'GRasia$', 'GRadmix$',
           'MaxEstLhood', 'MaxObsLhood']
results = results.loc[:, columns]
results.to_csv("results.csv")

results = results.loc[results["MaxEstLhood"] < 0, :]
idxmax = results.groupby(['group'])["MaxEstLhood"].transform(max) == results["MaxEstLhood"]
data = results.loc[idxmax, columns[1:-2]]
data.to_csv("data.csv")
