import numpy as np
import gzip

# genotypes could be an MxN matrix or generator of size N vectors
# loci is a size M vector
write_relate_input(genotypes, loci, names):
  haps_file = gzip.open(name + ".haps", "wt")
  for i, locus in enumerate(loci):
    string = "1 snp" + str(i + 1) + " " + str(locus) + " A" + " T "
    string = string + " ".join(map(str, next(genotypes))) + "\n"
