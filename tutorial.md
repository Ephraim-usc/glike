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

gLike provides a function `maximize` that mimics the popular `scipy.optimize.minimize` function, but has been made convenient for optimizing demographic parameters in several ways. The most important grammatical difference is that `glike.maximize` works with named parameters. Specifically,

1. It accepts a dict for `x0` (rather than a 1D array), so that the function is called in the `fun(**x)` way.
2. It accepts string expressions containing parameter names in `bounds`, which will be interpreted by `eval` during runtime.






There are several pieces of information needed to describe the estimation task.
`names` lists the names of each parameter

    names = ["t1", "t2", "t3", "r1", "r2", "N", "N_a", "N_b", "N_c", "N_d", "N_e"]

`values` gives an initial guess of the parameters

    values = [20, 80, 5000.0, 0.5, 0.5, 1000, 10000, 1000, 10000, 10000, 1000]

`limits` gives the list of tuples containing the lower and upper bound of each parameter

    limits = [(10, "t2"),("t1", 100),(1e3, 2e4),(0.001,0.999),(0.001,0.999),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000)]

The search task is the defined as

    search = glike.Search(names, values, limits, precision = 0.02)

Where `precision = 0.02` specifies the minimum step size in the hill-climbing optimization to be 2%.

The estimation is launched by

    x, logp = glike.estimate(trees, threeway_admixture_demo, search, prune = 0.5)

Which prints, for example
    
    [20, 80, 5000.0, 0.5, 0.5, 1000, 10000, 1000, 10000, 10000, 1000] -38918.96248255331
    [21.0, 78.0, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -38228.46344287755
    [22.65, 74.7, 6060.0, 0.4501, 0.48254, 1238.5, 12623.5, 1238.5, 12623.5, 12623.5, 1238.5] -37374.656869676255
    [25.49625, 69.0075, 7198.5, 0.4501, 0.48254, 1494.6625, 15441.2875, 1494.6625, 15441.2875, 15441.2875, 1494.6625] -36433.658071892474
    [30.72623, 58.54753, 9290.49375, 0.43326, 0.53671, 1965.36109, 20618.97203, 1965.36109, 20618.97203, 15441.2875, 1965.36109] -35762.026886215186
    [30.72623, 58.54753, 9290.49375, 0.40895, 0.53671, 1965.36109, 20618.97203, 2898.04164, 20618.97203, 12852.44523, 2898.04164] -35733.769256650114
    [30.72623, 58.54753, 9290.49375, 0.37453, 0.57572, 1965.36109, 20618.97203, 2898.04164, 25748.71504, 12852.44523, 4297.06246] -35717.3123988045
    [30.72623, 62.02519, 9290.49375, 0.32726, 0.57572, 1965.36109, 18054.10053, 3597.55205, 25748.71504, 12852.44523, 6395.59369] -35709.98502238328
    [30.72623, 62.02519, 9808.64961, 0.32726, 0.57572, 2081.94616, 18054.10053, 3597.55205, 30557.84911, 13659.4359, 6395.59369] -35701.13391522842
    [30.07854, 59.03019, 9808.64961, 0.29629, 0.56233, 2081.94616, 16370.90361, 3597.55205, 30557.84911, 12372.34882, 4821.69527] -35697.27254368971
    [30.07854, 59.03019, 9808.64961, 0.29629, 0.54161, 1989.04243, 16370.90361, 3269.65655, 30557.84911, 14119.72036, 6592.331] -35696.90557883543
    [30.07854, 59.03019, 9808.64961, 0.29629, 0.57417, 1989.04243, 15226.8557, 3269.65655, 30557.84911, 17113.97118, 6592.331] -35695.66375478204
    [30.07854, 60.04802, 9808.64961, 0.3068, 0.61954, 1989.04243, 13631.44514, 3269.65655, 29487.06535, 17113.97118, 4969.24825] -35691.35596414065
    [30.07854, 60.04802, 9984.8226, 0.32313, 0.55876, 2026.82328, 15772.16205, 3269.65655, 27937.35683, 17113.97118, 4969.24825] -35689.053793104555
    [30.07854, 60.04802, 9984.8226, 0.29733, 0.55876, 2026.82328, 15772.16205, 3206.26342, 27937.35683, 18476.62578, 4969.24825] -35688.672897231605
    [30.07854, 60.04802, 9984.8226, 0.29733, 0.55876, 1988.28681, 15772.16205, 3206.26342, 27937.35683, 16268.94069, 4969.24825] -35690.18506807915
    [30.07854, 60.04802, 9984.8226, 0.29733, 0.5852, 1988.28681, 16701.93143, 3268.38869, 27380.60969, 19182.6363, 5197.49426] -35688.99700226109
    ...

The estimation usually converges after 10 ~ 30 rounds, which takes around 2 hours.
When the function finishes, `x` will be the estimated parameters, and `logp` will be the maximum likelihood ever reached.
