
import numpy as np
import matplotlib.pyplot as plt

def plot_hist():
  data = np.loadtxt("ufdata1.csv")

  fig, ax = plt.subplots(1,1)

  ax.hist(data, bins=30, histtype="step", density=True)

  fig.savefig("hist.png")
  pass

if __name__=="__main__":
  plot_hist()
  pass
