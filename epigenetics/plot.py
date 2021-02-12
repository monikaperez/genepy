import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plotAverageOfSamples(samples, folder="", showAll=False, maxv=None, minv=None):
  res = [] 
  plt.figure()
  plt.ylim(minv,maxv)
  for sample in samples:
    data = pd.read_csv(sample, sep='\t', skiprows=1, header=None, names=['chr', 'start', 'end', 'name', "foldchange","."]+list(range(600)))
    r = data[list(range(600))].mean().tolist()
    res.append(r)
    if showAll:
      sns.lineplot(data=np.array(r), color="#BFBFFF")
  sns.lineplot(data=np.array(res).mean(0))
  if folder:
    plt.savefig(folder+"_averageofsamples.pdf", color="#1F1FFF")
  return res
