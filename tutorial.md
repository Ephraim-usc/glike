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


Full likelihood of genealogical trees
------------

Core functionality

    logp = glike.glike_trees(trees, demo, samples, prune)
    
Where `trees` is any enumerable that contains `tskit` genealogical trees.
`demo` is the hypothesized Demography created manually or from provided models in `models.py`.
`samples` is the dict that contains sample:population pairs, which specifies which sample is collected from which population.
`prune` is the proportion of discarding low likelihood trees, this often reduces noise when dealing with reconstructed trees.

This function returns the probability that such genealogical trees are generated under the hypothesized demography.


Demography customization
------------

A Demography object is initialized with

    demo = glike.Demo()
    
And a number of Phases are added into this Demography

    demo.add_phase(phase)

A Phase is created by

    phase = glike.Phase(t, ns, grs, P, Q, populations)

Where `t` is the starting timee, `ns` is the vector of coalescent rates, `grs` is the vector of growth rates, `P` is the mass migration matrix at the beginning of this phase, `Q` is the continuous migration rate matrix, and `populations` is the vector of population names. Only `t` and `ns` are required , other arguments are optional. The number of populations in this phase should be consistent among parameters, so it is required that

    len(ns) == len(grs) == P.shape[1] == Q.shape[1] == len(populations)

When adding new Phases into Demogrpahy, the times and dimensions should match. Specifically, if the last added phase is `phases[i]`, and we are trying to add another phase `phases[i+1]`, then it is required that

    phases[i].t < phases[i+1].t
    phases[i].P.shape[1] == phases[i+1].P.shape[0]

gLike cannot be applied directly to a demography that contains continuous migrations, and requires it to be discretized first:

    demo = demo.discretize(delta)

Where `delta` is a float number indicating the time interval (that is, any `phase` with `phase.Q != None` will be sliced into chunks of length `delta`).

The resulting demography can be visualized by

    demo.print()


Demography parameter estimation
------------

To make a demographic model containing variable parameters, the idiom is

    def model(...):
      demo = glike.Demography()
      # add Phases that depend on the parameters
      return demo

We provide an function for estimating parameters 

    glike.estimate(trees, model, search)

Which runs a smart maximum likelihood protocol specially designed for glike to find the estimated parameters.
We use a Search object to tell the information about initial values and restrictions of the parameters. It is created with

    search = Search(names, values, limits, fixed)

Where `names` is a list of the names of the parameters, `values` is a list of initial parameter values, `fixed` is a list of the names of fixed parameters so that their values will keep untouched. `limits` is a list of tuples `(low, high)` where `low` and `high` could be a number or the name of another parameter. For example, if our model has three events happening between 0 and 100 generations ago, in order to estimate the times of these three events we could do

    names = ["t1", "t2", "t3"]
    values = [25, 50, 75] # or other initial values you see fit
    limits = [(0, "t2"), ("t1", "t3"), ("t2", 100)]

Empirically, if the output parameters take values very close to the lower or upper limits, it's likely that the estimation is stuck in a local optimal, or the proposed model is not well compatible with the genealogical trees. If that's the case, it's suggested to try other initial values or demography models.


About this project
-------------

This is an ongoing project, please refer to our poster at ASHG 2022 for a brief introduction.

If you are interested in knowing more, please let us know. Email the author: caoqifan@usc.edu

![](images/poster_ashg.png)
