import math
import numpy as np
import pandas as pd
import gzip
from tqdm import tqdm

import tsinfer


def write_relate_input(arg, name):
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
    
    string = str(position) + " " + str(1.0) + " "
    string = string + str(variant.position / 1000000) + "\n"
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
  
