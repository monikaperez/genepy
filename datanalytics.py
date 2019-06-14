import numpy as np


def getDFinfo(df):
  val = df
  print("sums over cell lines! ------mean, var, totmin, meanmin, totmax, meanmax")
  print(val.sum(1).mean(), val.sum(1).var(), val.sum(1).min(), val.mean(1).min(), val.sum(1).idxmin(), val.sum(1).max(), val.mean(1).max(), val.sum(1).idxmax())
  print("sums over features! ------mean, var, totmin, meanmin, totmax, meanmax")
  print(val.sum(0).mean(), val.sum(0).var(), val.sum(0).min(), val.mean(0).min(), val.sum(0).idxmin(), val.sum(0).max(), val.mean(0).max(), val.sum(0).idxmax())
  print("nans!")
  print(np.count_nonzero(np.isnan(val)))
