import numpy as np
import gzip


def write_relate_input(arg, name):
  haps_file = gzip.open(name + ".haps.gz", "wt")
  for variant in arg.variants():
    string = "1 snp" + str(variant.index + 1) + " " + str(variant.position) + " A" + " T "
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
  for variant in arg.variants():
    string = str(variant.position) + " " + str(1.0) + " "
    string = string + str(variant.position / 1000000) + "\n"
    map_file.write(string)
  map_file.close()
