Tutorial
========

This is a toy example that replicates the experiment in Figure 2B of the paper.


Simulating the ARG by msprime
------------

The three-way admixture model as defined in the msprime language is

    import msprime
    
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

    demography = threeway_admixture_demography(30, 60, 1e4, 0.4, 0.7, 2000, 20000, 3000, 30000, 10000, 5000)

We simulate a 30Mb chromosome and select 10 equally distant trees

    arg = msprime.sim_ancestry({"O": 1000}, sequence_length = 3e7, recombination_rate = 1e-8, demography = demography, ploidy = 1)
    arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
    trees = [arg.at(pos).copy() for pos in range(3000000, 30000000, 3000000)]


Defining the demographic model
------------

The three-way admixture model as defined in the glike language is

    def threeway_admixture_demo(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e):
      demo = Demo()
      demo.add_phase(Phase(0, [1/N]))
      demo.add_phase(Phase(t1, [1/N_a, 1/N_b], P = np.array([[r1, 1-r1]])))
      demo.add_phase(Phase(t2, [1/N_a, 1/N_c, 1/N_d], P = np.array([[1, 0, 0], [0, r2, 1-r2]])))
      demo.add_phase(Phase(t3, [1/N_e], P = np.array([[1], [1], [1]])))
      return demo

Where the demography consists of four phases, each defined by the starting time, the list of inverse population sizes, and the mass migration matrices.


Estimating parameters
------------

There are several pieces of information needed to describe the estimation task.
`names` lists the names of each parameter

    names = ["t1", "t2", "t3", "r1", "r2", "N", "N_a", "N_b", "N_c", "N_d", "N_e"]

`values` gives an initial guess of the parameters

    values = [20, 80, 5000.0, 0.5, 0.5, 1000, 10000, 1000, 10000, 10000, 1000]

`limits` gives the list of tuples containing the lower and upper bound of each parameter

    limits = [(10, "t2"),("t1", 100),(1e3, 2e4),(0.001,0.999),(0.001,0.999),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000)]

The search task is the defined as

    search = Search(names, values, limits, precision = 0.02)

Where `precision = 0.02` specifies the minimum step size in the hill-climbing optimization to be 2%.

The estimation is launched by

    x, logp = estimate(trees, threeway_admixture_demo, search, prune = 0.5)

Which prints, for example
    
    [20, 80, 5000.0, 0.5, 0.5, 1000, 10000, 1000, 10000, 10000, 1000] -38919.783993402685
    [20, 80, 5400.0, 0.5, 0.5, 1000, 10000, 1000, 10000, 10000, 1000] -38693.018071448445
    [20, 80, 5400.0, 0.4501, 0.5, 1000, 10000, 1000, 10000, 10000, 1000] -38518.86148675874
    [20, 80, 5400.0, 0.4501, 0.5499, 1000, 10000, 1000, 10000, 10000, 1000] -38511.40311310784
    [20, 80, 5400.0, 0.4501, 0.5499, 1090.0, 10000, 1000, 10000, 10000, 1000] -38245.01574334988
    [20, 80, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1000, 10000, 10000, 1000] -38231.61911714237
    [20, 80, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10000, 10000, 1000] -38212.70205963768
    [20, 80, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10000, 1000] -38210.17570146851
    [20, 80, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1000] -38207.1929041596
    [20, 80, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -38078.97551290889
    [20.5, 80, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -38053.10669115948
    [20.5, 79.0, 5400.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -38050.854700433134
    [20.5, 79.0, 6060.0, 0.4501, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -37721.66517975466
    [20.5, 79.0, 6060.0, 0.38273, 0.5499, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -37713.270123500275
    [20.5, 79.0, 6060.0, 0.38273, 0.61472, 1090.0, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -37711.62702384707
    [20.5, 79.0, 6060.0, 0.38273, 0.61472, 1238.5, 10990.0, 1090.0, 10990.0, 10990.0, 1090.0] -37396.84843228326
    ...

The estimation usually converges after 10~30 rounds (which means 110~330 lines of printed output).
When the function finishes, `x` will be the estimated parameters, and `logp` will be the maximum likelihood ever reached.
