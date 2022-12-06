from glike import *
import tsinfer
import tsdate


N = 1000
l = 5000000

demography = threeway_admixture_demography(20, 50, 0.4, 0.7, 2000, 10000, 3000, 20000, 15000, 5000, 5e-3, 1e-3)
arg = msprime.sim_ancestry({"O": N}, sequence_length = l, recombination_rate = 1e-10, demography = demography, ploidy = 1)
arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
arg.dump("arg.trees")




write_relate_input(arg, "twa")


'''
~/bin/Relate --mode All -m 1.25e-8 -N 30000 --memory 20 --haps twa.haps --sample twa.sample --map twa.map -o twa
~/bin/RelateFileFormats --mode ConvertToTreeSequence -i twa -o twa
'''

arg_relate = tskit.load("twa.trees")
trees_relate = [arg_relate.at(pos).copy() for pos in range(0, l, 200000)]



write_tsinfer_input(arg, 'twa')
'''
~/bin/tsinfer infer twa.samples -p -t 4 --recombination-rate 1e-10
'''

arg_tsinfer = tskit.load("twa.trees")
arg_tsinfer_simplified = arg_tsinfer.simplify(filter_populations = False, filter_individuals = False, filter_sites = False, keep_unary = False)
arg_tsdate = tsdate.date(arg_tsinfer_simplified, Ne=14587, mutation_rate = 1e-8)
arg_tsdate.dump("arg_tsdate.trees")

trees = [arg.at(pos).copy() for pos in range(0, l, 1000000)]
trees_tsdate = [arg_tsdate.at(pos).copy() for pos in range(100000, l, 200000)]




names = ["t1", "t2", "r1", "r2", "N", "N_a", "N_b", "N_c", "N_d", "N_e", "m_ab", "m_cd"]
values = [40, 70, 0.4, 0.7, 2000, 10000, 3000, 20000, 15000, 5000, 0, 0]
limits = [(0,"t2"),("t1",100),(0,1),(0,1),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(0,0),(0,0)]
fixed = ["m_ab", "m_cd"]

searchspace = Searchspace(names, values, limits, fixed)
estimate(trees_tsdate, threeway_admixture_demo, searchspace)
estimate(trees_tsdate, threeway_admixture_soft_demo, searchspace)






