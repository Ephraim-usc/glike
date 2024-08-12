import math
import numpy as np
import pandas as pd
import gzip
from tqdm import tqdm

import tsinfer


def write_relate_input(arg, name, recomb_rate = 1e-8):
  haps_file = open(name + ".haps", "wt")
  prev_position = 0
  for variant in tqdm(arg.variants(), total = arg.num_mutations):
    if variant.genotypes.max() > 1:
      continue
    position = math.ceil(variant.position)
    if position == prev_position:
      continue
    prev_position = position
    
    string = "1 snp" + str(variant.index + 1) + " " + str(position) + " A" + " T "
    string = string + " ".join(map(str, variant.genotypes)) + "\n"
    haps_file.write(string)
  haps_file.close()
  
  sample_file = open(name + ".sample", "wt")
  sample_file.write("ID_1 ID_2 missing\n0 0 0\n")
  for sample in range(arg.num_samples):
    string = "UNR" + str(sample + 1) + " NA" + " 0\n"
    bytes = sample_file.write(string)
  sample_file.close()

  map_file = open(name + ".map", "wt")
  map_file.write("pos COMBINED_rate Genetic_Map\n")
  prev_position = 0
  for variant in arg.variants():
    position = math.ceil(variant.position)
    if position == prev_position:
      continue
    prev_position = position
    
    string = str(position) + " " + str(recomb_rate * 1e8) + " "
    string = string + str(variant.position * recomb_rate * 1e2) + "\n"
    map_file.write(string)
  map_file.close()


def write_tsinfer_input(arg, name):
  sample_data = tsinfer.SampleData(path = name + ".samples", sequence_length = math.ceil(arg.last().interval[1]))
  first = arg.first()
  for sample in arg.samples():
    sample_data.add_individual(time = first.time(sample))
  
  prev_position = 0
  for variant in tqdm(arg.variants(), total = arg.num_mutations):
    if variant.genotypes.max() > 1:
      continue
    position = math.ceil(variant.position)
    if position == prev_position:
      continue
    prev_position = position
    sample_data.add_site(position, variant.genotypes)
  
  sample_data.finalise()
  return sample_data


def plot_tree(tree):
  x = {}
  x[tree.root] = 0
  
  nodes = [tree.root]
  while nodes:
    node = nodes.pop()
    shift = x[node] - len(list(tree.samples(node)))/2
    prev = 0
    for child in tree.children(node):
      tmp = len(list(tree.samples(child)))/2
      x[child] = shift + prev + tmp
      shift = x[child]
      prev = tmp
      nodes.append(child)
  
  lines = []
  for node in tree.nodes():
    if tree.parent(node) != -1:
      lines.append((x[node], tree.time(node), x[node], tree.time(tree.parent(node))))
    if tree.children(node):
      xs = [x[child] for child in tree.children(node)]
      lines.append((min(xs), tree.time(node), max(xs), tree.time(node)))
  
  lines = pd.DataFrame(lines, columns = ("x", "y", "xend", "yend"))
  return lines


'''
# this function converts an msprime demography object into a glike Demo object
# this works for msprime v1.2.0 but may not always work if msprime updates
# note that all continuous migrations will be ignored
def demography_to_demo(demography):
    import msprime
    import pandas as pd
    demo = Demo()
    events = demography.events
    
    t = 0
    t_end = events[0].time if len(events) > 0 else np.inf
    ns = [1/population.initial_size for population in demography.populations]
    grs = [population.growth_rate for population in demography.populations]
    populations = [population.name for population in demography.populations]
    demo.add_phase(Phase(0, t_end, ns, grs, populations = populations))
    
    while event = demography.events.pop(0):
      if type(event) is msprime.demography.Admixture:
        populations_ = populations; populations = populations.copy(); 
        populations.remove(event.derived)
        P = pd.DataFrame(0, index = populations_, columns = populations)
        for population in populations:
          P.loc[population, population] = 1
        for anc, prop in zip(event.ancestral, event.proportions):
          P.loc[event.derived, anc] = prop
        P = P.values
      elif type(event) is msprime.demography.PopulationSplit:
        populations_ = populations; populations = populations.copy()
        for derived in event.derived:
          populations.remove(derived)
        P = pd.DataFrame(0, index = populations_, columns = populations)
        for population in populations:
          P.loc[population, population] = 1
        for derived in event.derived:
          P.loc[derived, event.ancestral] = 1
        P = P.values
      elif type(event) is msprime.demography.PopulationParametersChange:
        populations_ = populations; populations = populations.copy()
        P = None
      else:
        continue
      
      t_ = t; t = event.time
      t_end = events[0].time if len(events) > 0 else np.inf
      ns_ = ns.copy()
      ns = [None for _ in populations]
      for population in populations:
        n = ns_[populations_.index(population)]
        if type(n) is float:
          ns[populations.index(population)] = n
        else:
          ns[populations.index(population)] = (math.log(n[0]) + n[1] * (t-t_), n[1])
      
      if type(event) is msprime.demography.PopulationParametersChange:
        if event.population == -1:
          tmp = populations.copy()
        else:
          tmp = [event.population]
        for population in tmp:
          ns[populations.index(population)] = (1/event.initial_size, event.growth_rate) if event.growth_rate and event.growth_rate>0 else 1/event.initial_size
'''


# this function converts a glike Demo object into an msprime demography object
# this works for msprime v1.2.0 but may not always work if msprime updates
def demo_to_demography(demo):
  import msprime
  import pandas as pd
  
  # find all relevant populations and their initial sizes
  populations = []
  ns = []
  grs = []
  for phase in demo.phases:
    for population, n, gr in zip(phase.populations, phase.ns, phase.grs):
      if population not in populations:
        populations.append(population)
        ns.append(n)
        grs.append(gr)
  
  demography = msprime.Demography()
  for population, n, gr in zip(populations, ns, grs):
    demography.add_population(name = population, initial_size = 1/n, growth_rate = gr)
  
  for phase_, phase in zip(demo.phases[:-1], demo.phases[1:]):
    P = phase.P
    if (P.shape[0] == P.shape[1]) and (P == np.eye(P.shape[0])).all():
      continue
    P = pd.DataFrame(P, index = phase_.populations, columns = phase.populations)
    for source in phase_.populations:
      if source in phase.populations and P.loc[source, source] == 1:
        continue
      dests, proportions = [], []
      for dest in phase.populations:
        if P.loc[source, dest] > 0:
          dests.append(dest); proportions.append(P.loc[source, dest])
      if len(dests) == 1:
        demography.add_mass_migration(time = phase.t, source = source, dest = dests[0], proportion = 1)
      else:
        demography.add_admixture(time = phase.t, derived = source, ancestral = dests, proportions = proportions)
    for population, n, gr in zip(phase.populations, phase.ns, phase.grs):
      demography.add_population_parameters_change(time = phase.t, initial_size = 1/n, growth_rate = gr, population = population)
  
  return demography


# the coalescent times of a tree, in ascending order
# returns an 1D array of length N-1
def get_coalescent_times(tree):
  times = [tree.time(node) for node in tree.nodes()]
  progenies =  [len(tree.children(node)) for node in tree.nodes()]
  times_coal = [time for progeny, time in zip(progenies, times) for i in range(progeny - 1) if progeny >= 2]
  times_coal = sorted(times_coal)
  return np.array(times_coal)

def get_coalescent_times_trees(trees):
  return np.array([get_coalescent_times(tree) for tree in trees])

def get_coalescent_times_demo(demo, samples_msprime, sims = 10000):
  import msprime
  demography = demo_to_demography(demo)
  trees = [msprime.sim_ancestry(samples_msprime, sequence_length = 1, demography = demography, ploidy = 1).first() for _ in range(sims)]
  return get_coalescent_times_trees(trees)
