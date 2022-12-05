from glike import *

N = 1000
l = 30000000

demography = threeway_admixture_demography(20, 50, 0.4, 0.7, 2000, 10000, 3000, 20000, 15000, 5000, 5e-3, 1e-3)
arg = msprime.sim_ancestry({"O": N}, sequence_length = l, demography = demography, ploidy = 1)
arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)

trees = [arg.at(pos).copy() for pos in range(0, l, 1000000)]



names = ["t1", "t2", "r1", "r2", "N", "N_a", "N_b", "N_c", "N_d", "N_e", "m_ab", "m_cd"]
values = [20, 50, 0.4, 0.7, 2000, 10000, 3000, 20000, 15000, 5000, 0, 0]
limits = [(0,"t2"),("t1",100),(0,1),(0,1),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(0,0),(0,0)]
fixed = ["m_ab", "m_cd"]

searchspace = Searchspace(names, values, limits, fixed)
estimate(trees, threeway_admixture_demo, searchspace)
