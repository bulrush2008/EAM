
import numpy as np
import matplotlib.pyplot as plt

def plot_hist():
  data1 = np.loadtxt("ufdata1.csv")

  # RR 分布参数
  DBar = 0.04
  expn = 3
  dmin = 1.e-4
  dmax = 0.1

  K = 1.0 - np.exp((dmin/DBar)**expn - (dmax/DBar)**expn)

  diams = DBar * ((dmin/DBar)**expn - np.log(1.0 - K*data1))**(1.0/expn)

  fig, axes = plt.subplots(2,2, figsize=(16,16))

  axes[0][0].hist(data1, bins=100, histtype="step", density=True, label="pdf-uniform")
  axes[0][1].hist(data1, bins=100, histtype="step", density=True, cumulative=True, label="cdf-uniform")
  axes[0][0].legend()
  axes[0][1].legend()

  axes[1][0].hist(diams, bins=100, histtype="step", density=True)
  axes[1][1].hist(diams, bins=100, histtype="step", density=True, cumulative=True)

  fig.savefig("hist.png")
  pass

if __name__=="__main__":
  plot_hist()
  pass
