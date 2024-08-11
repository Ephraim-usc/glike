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


# this function converts a glike Demo object into an msprime demography object
# this works for msprime v1.2.0 but may not always work if msprime updates
def demo_to_demography(demo):
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
      
      demo.add_phase(Phase(t, ns, P = P, populations = populations))
