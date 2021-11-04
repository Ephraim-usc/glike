from .glike import *

def single_lmp(N):
  times = [0]
  mss = [np.array([[0]])]
  nss = [[1/N]]
  lmp = glike.LMP(times, mss, nss)
  return lmp


def twoway_admixture_lmp(t, r, N_ab, N_a, N_b, m = 1e-4):
  times = [0, t]
  mss = [np.array([[0, 0, 0],
                   [0, 0, m],
                   [0, m, 0]])] * 2
  nss = [[1/N_ab, 1/N_a, 1/N_b]] * 2
  Ps = [np.identity(5),
        np.matrix([[0, r, 1-r], [0, 1, 0], [0, 0, 1]])]
  lmp = LMP(times, mss, nss, Ps = Ps)
  return lmp


def threeway_admixture_lmp(t1, t2, r1, r2, N_abc, N_ab, N_a, N_b, N_c, m = 1e-4):
  times = [0, t1, t1 + t2]
  mss = [np.array([[0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, m, m],
                   [0, 0, m, 0, m],
                   [0, 0, m, m, 0]])] * 3
  nss = [[1/N_abc, 1/N_ab, 1/N_a, 1/N_b, 1/N_c]] * 3
  Ps = [np.identity(5),
        np.matrix([[0, r1, 0, 0, 1-r1], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]),
        np.matrix([[1, 0, 0, 0, 0], [0, 0, r2, 1-r2, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])]
  lmp = LMP(times, mss, nss, Ps = Ps)
  return lmp
