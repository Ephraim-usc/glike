//Number of population samples (demes)
6
//Population effective sizes (number of genes)
Nadmix$
Nafr$
Neur$
Nasia$
Nooa$
Nanc$
//Sample sizes
50
25
25
25
0
0
//Growth rates: negative growth implies population expansion
GBRadmix$
0
GBReur$
GBRasia$
0
0
//Number of migration matrices: 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
7 historical event 
T1$ 0 1 R1$ 1 0 0
T1$ 0 2 R4$ 1 0 0
T1$ 0 3 1 1 0 0
T2$ 2 4 1 1 0 0
T2$ 3 4 1 1 0 0
T3$ 4 1 1 1 0 0
T4$ 1 5 0 1 0 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA 30000000 1.0e-8 1.0e-8
