from .glike import *

def threeway_admixture_lmp(t1, t2, r1, r2, n = 1e-4, m = 1e-4):
  times = [0, t1, t1 + t2]
  mss = [np.array([[0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, m, m],
                   [0, 0, m, 0, m],
                   [0, 0, m, m, 0]])] * 3
  nss = [[n, n, n, n, n]] * 3
  Ps = [np.identity(5),
        np.matrix([[0, r1, 0, 0, 1-r1], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]),
        np.matrix([[1, 0, 0, 0, 0], [0, 0, r2, 1-r2, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])]
  lmp = LMP(times, mss, nss, Ps = Ps)
  return lmp
