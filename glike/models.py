import msprime
from .glike import *

def single_lmp(N):
  times = [0]
  mss = [np.array([[0]])]
  nss = [[1/N]]
  lmp = LMP(times, mss, nss)
  return lmp

def single_demography(N):
  demography = msprime.Demography()
  demography.add_population(name = "A", initial_size = N)
  return demography


def popsize_variation_lmp(ts, Ns):
  times = ts
  mss = [np.array([[0]]) for N in Ns]
  nss = [[1/N] for N in Ns]
  lmp = LMP(times, mss, nss)
  return lmp

def popsize_variation_demography(ts, Ns):
  demography = msprime.Demography()
  demography.add_population(name = "A", initial_size = Ns[0])
  for t, N in zip(ts[1:], Ns[1:]):
    demography.add_population_parameters_change(time = t, initial_size = N, population = "A")
  return demography


def twoway_admixture_lmp(t, r, N, N_a, N_b, N_ab):
  times = [0, t, 1e5]
  mss = [np.zeros((4, 4))] * 3
  nss = [[1/N, 1/N_a, 1/N_b, 1/N_ab]] * 3
  Ps = [np.identity(4),
        np.matrix([[0, r, 1-r, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
        np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 1]])]
  lmp = LMP(times, mss, nss, Ps = Ps)
  return lmp

def twoway_admixture_demography(t, r, N, N_a, N_b, N_ab):
  demography = msprime.Demography()
  demography.add_population(name = "O", initial_size = N)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "AB", initial_size = N_ab)
  
  demography.add_admixture(time=t, derived="O", ancestral=["A", "B"], proportions = [r, 1-r])
  demography.add_population_split(time=1e5, derived=["A", "B"], ancestral="AB")
  return demography


'''
demography = twoway_admixture_demography(10, 0.7, 2000, 10000, 20000, 5000)
trees = msprime.sim_ancestry({"O": 100}, sequence_length = 1e6, recombination_rate = 1e-8, 
                         demography = demography, ploidy = 1)
tree = trees.first()
labels = ['0'] * 100

loglike_tree(tree, labels, twoway_admixture_lmp(10, 0.6, 2000, 10000, 20000, 5000))[2]
loglike_tree(tree, labels, twoway_admixture_lmp(10, 0.7, 2000, 10000, 20000, 5000))[2]
loglike_tree(tree, labels, twoway_admixture_lmp(10, 0.8, 2000, 10000, 20000, 5000))[2]

loglike_tree(tree, labels, twoway_admixture_lmp(10, 0.7, 1500, 10000, 20000, 5000))[2]
loglike_tree(tree, labels, twoway_admixture_lmp(10, 0.7, 2000, 10000, 20000, 5000))[2]
loglike_tree(tree, labels, twoway_admixture_lmp(10, 0.7, 2500, 10000, 20000, 5000))[2]


'''




def twoway_admixture_lmp(t, r, N_ab, N_a, N_b, m = 1e-4):
  times = [0, t]
  mss = [np.array([[0, m, m],
                   [m, 0, m],
                   [m, m, 0]]),
         np.array([[0, 0, 0],
                   [0, 0, m],
                   [0, m, 0]])]
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
  demography.set_symmetric_migration_rate(["AB", "A"], m)
  demography.set_symmetric_migration_rate(["AB", "B"], m)
  
  demography.add_admixture(time=t, derived="AB", ancestral=["A", "B"], proportions = [r, 1-r])
  demography.add_symmetric_migration_rate_change(time=t, populations=["AB","A"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t, populations=["AB","B"], rate=0)
  return demography


def threeway_admixture_lmp(t1, t2, r1, r2, N_abc, N_ab, N_a, N_b, N_c, m = 1e-4):
  times = [0, t1, t2]
  mss = [np.array([[0, 0, m, m, m],
                   [0, 0, 0, 0, 0],
                   [m, 0, 0, m, m],
                   [m, 0, m, 0, m],
                   [m, 0, m, m, 0]]),
         np.array([[0, 0, 0, 0, 0],
                   [0, 0, m, m, m],
                   [0, m, 0, m, m],
                   [0, m, m, 0, m],
                   [0, m, m, m, 0]]),
         np.array([[0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, m, m],
                   [0, 0, m, 0, m],
                   [0, 0, m, m, 0]])]
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
  demography.set_symmetric_migration_rate(["ABC", "A"], m)
  demography.set_symmetric_migration_rate(["ABC", "B"], m)
  demography.set_symmetric_migration_rate(["ABC", "C"], m)
  demography.set_symmetric_migration_rate(["AB", "A"], 0)
  demography.set_symmetric_migration_rate(["AB", "B"], 0)
  demography.set_symmetric_migration_rate(["AB", "C"], 0)
  
  demography.add_admixture(time=t1, derived="ABC", ancestral=["AB", "C"], proportions = [r1, 1-r1])
  demography.add_symmetric_migration_rate_change(time=t1, populations=["ABC","A"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["ABC","B"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["ABC","C"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["AB","A"], rate=m)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["AB","B"], rate=m)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["AB","C"], rate=m)
  
  demography.add_admixture(time=t2, derived="AB", ancestral=["A", "B"], proportions = [r2, 1-r2])
  demography.add_symmetric_migration_rate_change(time=t2, populations=["AB","A"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t2, populations=["AB","B"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t2, populations=["AB","C"], rate=0)
  return demography




def twoway_split_lmp(t, N_ab, N_a, N_b, m = 1e-4):
  times = [0, t]
  mss = [np.array([[0, 0, 0],
                   [0, 0, m],
                   [0, m, 0]]),
         np.array([[0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0]])]
  nss = [[1/N_ab, 1/N_a, 1/N_b]] * 2
  Ps = [np.identity(3),
        np.matrix([[1, 0, 0], [1, 0, 0], [1, 0, 0]])]
  lmp = LMP(times, mss, nss, Ps = Ps)
  return lmp

def twoway_split_demography(t, N_ab, N_a, N_b, m = 1e-4):
  demography = msprime.Demography()
  demography.add_population(name = "AB", initial_size = N_ab)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  
  demography.set_symmetric_migration_rate(["A", "B"], m)
  demography.add_symmetric_migration_rate_change(time=t, populations=["A","B"], rate=0)
  
  demography.add_population_split(time=t, derived=["A", "B"], ancestral="AB")
  return demography


def threeway_split_lmp(t1, t2, N_abc, N_ab, N_a, N_b, N_c, m = 1e-4):
  times = [0, t1, t2]
  mss = [np.array([[0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, m, m],
                   [0, 0, m, 0, m],
                   [0, 0, m, m, 0]]),
         np.array([[0, 0, 0, 0, 0],
                   [0, 0, 0, 0, m],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, m, 0, 0, 0]]),
         np.array([[0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0]])]
  nss = [[1/N_abc, 1/N_ab, 1/N_a, 1/N_b, 1/N_c]] * 3
  Ps = [np.identity(5),
        np.matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1]]),
        np.matrix([[1, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0]])]
  lmp = LMP(times, mss, nss, Ps = Ps)
  return lmp

def threeway_split_demography(t1, t2, N_abc, N_ab, N_a, N_b, N_c, m = 1e-4):
  demography = msprime.Demography()
  demography.add_population(name = "ABC", initial_size = N_abc)
  demography.add_population(name = "AB", initial_size = N_ab)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "C", initial_size = N_c)
  
  demography.set_symmetric_migration_rate(["A", "B"], m)
  demography.set_symmetric_migration_rate(["A", "C"], m)
  demography.set_symmetric_migration_rate(["B", "C"], m)
  demography.set_symmetric_migration_rate(["AB", "C"], 0)
  
  demography.add_population_split(time=t1, derived=["A", "B"], ancestral="AB")
  demography.add_symmetric_migration_rate_change(time=t1, populations=["A","B"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["A","C"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["B","C"], rate=0)
  demography.add_symmetric_migration_rate_change(time=t1, populations=["AB","C"], rate=m)
  
  demography.add_population_split(time=t2, derived=["AB", "C"], ancestral="ABC")
  demography.add_symmetric_migration_rate_change(time=t2, populations=["AB","C"], rate=0)
  return demography
