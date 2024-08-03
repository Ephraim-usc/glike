Tutorial
========

This is a toy example that replicates the experiment in Figure 2B of the paper.

Import necessary packages to go through this tutorial

    import numpy as np
    import msprime
    import glike


Simulating the ARG by msprime
------------

The three-way admixture model as defined in the msprime language is

    def threeway_admixture_demography(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e):
      demography = msprime.Demography()
      demography.add_population(name = "O", initial_size = N)
      demography.add_population(name = "A", initial_size = N_a)
      demography.add_population(name = "B", initial_size = N_b)
      demography.add_population(name = "C", initial_size = N_c)
      demography.add_population(name = "D", initial_size = N_d)
      demography.add_population(name = "E", initial_size = N_e)
      demography.add_admixture(time=t1, derived="O", ancestral=["A", "B"], proportions = [r1, 1-r1])
      demography.add_admixture(time=t2, derived="B", ancestral=["C", "D"], proportions = [r2, 1-r2])
      demography.add_population_split(time=t3, derived=["A", "C", "D"], ancestral="E")
      return demography

The true demography is created using the true parameters

    x_true = {"t1":30, "t2":60, "t3":1e4, "r1":0.4, "r2":0.7, "N":2000, "N_a":20000, "N_b":3000, "N_c":30000, "N_d":10000, "N_e":5000}
    demography = threeway_admixture_demography(**x_true)

We simulate 1000 haplotypes on a 30Mb chromosome and select 10 equally distant trees

    arg = msprime.sim_ancestry({"O": 1000}, sequence_length = 3e7, recombination_rate = 1e-8, demography = demography, ploidy = 1)
    arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
    trees = [arg.at(pos).copy() for pos in range(3000000, 30000000, 3000000)]


Defining the demographic model
------------

The three-way admixture model as defined in the glike language is

    def threeway_admixture_demo(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e):
      demo = glike.Demo()
      demo.add_phase(glike.Phase(0, t1, [1/N]))
      demo.add_phase(glike.Phase(t1, t2, [1/N_a, 1/N_b], P = np.array([[r1, 1-r1]])))
      demo.add_phase(glike.Phase(t2, t3, [1/N_a, 1/N_c, 1/N_d], P = np.array([[1, 0, 0], [0, r2, 1-r2]])))
      demo.add_phase(glike.Phase(t3, np.inf, [1/N_e], P = np.array([[1], [1], [1]])))
      return demo

Where the demography consists of four phases, each defined by the starting and ending time, the list of inverse population sizes, and the mass migration matrices.


Checking out
------------

It is recommended to always check if the demographic model is written correctly. Here we create the true demography as an example:

    demo = threeway_admixture_demo(**x_true)
    demo.print()

The output is

    [phase from 0.0 to 30.0]
              A
    ns   0.0005
    grs  0.0000
    
    [phase from 30.0 to 60.0]
         A    B
    A  0.4  0.6
               A         B
    ns   0.00005  0.000333
    grs  0.00000  0.000000
    
    [phase from 60.0 to 10000.0]
         A    B    C
    A  1.0  0.0  0.0
    B  0.0  0.7  0.3
               A         B       C
    ns   0.00005  0.000033  0.0001
    grs  0.00000  0.000000  0.0000
    
    [phase from 10000.0 to inf]
       A
    A  1
    B  1
    C  1
              A
    ns   0.0002
    grs  0.0000

Which lists all phases in this demography. 
In each phase, the migration matrix -- if applicable -- is printed first.
Then the coalescence rates (i.e. inverse of population sizes) and growth rates are printed.
Note that because we did not specify names of the populations, alphabetical names (A, B, C, ...) are automatically created.

We can now play with gLike by

    glike.glike_trees(trees, demo)

Which outputs something like

    -80165.92668862805

This is the probability of `demo` to generate the genealogical trees in `trees`.
Feel free to change the parameters a bit and see how the likelihood changes (drops in most cases).


Estimating parameters
------------

gLike provides a function `maximize` that mimics the popular `scipy.optimize.minimize` function, but has been made convenient for optimizing demographic parameters. The most important grammatical difference is that `glike.maximize` works with named parameters. For example,

1. It accepts a dict for `x0` (rather than a 1D array), so that the function is called in the `fun(**x)` way.
2. It accepts string expressions containing parameter names in `bounds`, which will be interpreted by `eval` during runtime.

To use `glike.maximize` for parameter estimation, we need to first define a wrapper function, the initial values, and boundaries

    def fun(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e):
      demo = threeway_admixture_demo(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e)
      return glike.glike_trees(trees, demo)
    
    x0 = {"t1":10, "t2":20, "t3":5000.0, "r1":0.5, "r2":0.5, "N":10000, "N_a":10000, "N_b":10000, "N_c":10000, "N_d":10000, "N_e":10000}
    bounds = [(1, "t2"),("t1", "t3"),("t2", 1e5),(0.0,1.0),(0.0,1.0),(100,1000000),(100,1000000),(100,1000000),(100,1000000),(100,1000000),(100,1000000)]

Then running

    glike.maximize(fun, x0, bounds = bounds)

Would generate

    xxx

This usually takes around 2 hours on a personal computer.
When the function finishes, `x` will be the estimated parameters, and `logp` will be the maximum likelihood ever reached.


A more complicated example
------------

More advanced features (samples from mulitple pupulations, ancient DNA samples, etc.) can be demonstrated by replicating the experiment in Figure 5B of the paper.




