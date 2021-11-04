import msprime
from .glike import *

def single_lmp(N):
  times = [0]
  mss = [np.array([[0]])]
  nss = [[1/N]]
  lmp = glike.LMP(times, mss, nss)
  return lmp

def single_demography(N):
  demography = msprime.Demography()
  demography.add_population(name = "A", initial_size = N)
  return demography


def twoway_admixture_lmp(t, r, N_ab, N_a, N_b, m = 1e-4):
  times = [0, t]
  mss = [np.array([[0, 0, 0],
                   [0, 0, m],
                   [0, m, 0]])] * 2
  nss = [[1/N_ab, 1/N_a, 1/N_b]] * 2
  Ps = [np.identity(3),
        np.matrix([[0, r, 1-r], [0, 1, 0], [0, 0, 1]])]
  lmp = LMP(times, mss, nss, Ps = Ps)
  return lmp

def twoway_admixture_demography(t, r, N_ab, N_a, N_b, m = 1e-4):
  demography = msprime.Demography()
  demography.add_population(name = "AB", initial_size = N_ab)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.set_symmetric_migration_rate(["A", "B"], m)
  demography.add_admixture(time=t, derived="AB", ancestral=["A", "B"], proportions = [r, 1-r])
  return demography


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

def threeway_admixture_demography(t1, t2, r1, r2, N_abc, N_ab, N_a, N_b, N_c, m = 1e-4):
  demography = msprime.Demography()
  demography.add_population(name = "ABC", initial_size = N_abc)
  demography.add_population(name = "AB", initial_size = N_ab)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "C", initial_size = N_c)
  demography.set_symmetric_migration_rate(["A", "B"], m)
  demography.set_symmetric_migration_rate(["A", "C"], m)
  demography.set_symmetric_migration_rate(["B", "C"], m)
  demography.add_admixture(time=t1, derived="ABC", ancestral=["AB", "C"], proportions = [r1, 1-r1])
  demography.add_admixture(time=t1+t2, derived="AB", ancestral=["A", "B"], proportions = [r2, 1-r2])
  return demography

