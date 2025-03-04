
import numpy as np
import matplotlib.pyplot as plt

def plot_hist():
  data1 = np.loadtxt("ufdata1.csv")
  data2 = np.loadtxt("ufdata2.csv")

  fig, (ax1, ax2) = plt.subplots(1,2)

  ax1.hist(data1, bins=30, histtype="step", density=True)
  ax2.hist(data2, bins=30, histtype="step", density=True, log=True)

  fig.savefig("hist.png")
  pass

if __name__=="__main__":
  plot_hist()
  pass
