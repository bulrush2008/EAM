"""
方式1：一次性生成所有的随机数，备用
"""

import numpy as np
from matplotlib import pyplot as plt

class GenRandSeries(object):
  def __init__(self, N:int):
    self.data = np.zeros(N)

    self.N = N
    pass

  def gen(self):
    N = self.N
    self.data = np.random.rand(N)
    pass

  def plotHist(self):
    fig, ax = plt.subplots(1,1)

    d = self.data
    ax.hist(self.data, bins=30, histtype="step", density=True)

    fig.savefig("./hist.png")
    pass
  pass

if __name__=="__main__":
  g = GenRandSeries(1000)

  g.gen()

  g.plotHist()

  pass