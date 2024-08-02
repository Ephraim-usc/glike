gLike
========

Genealogical Likelihood (gLike) is a maximum likelihood method to infer the demographic history
of given populations that best explains the observed genealogical relationships between sample haplotypes. 
This is, to our knowledge, the first attempt to infer all types of parameters
(split times, admixture times and proportions, migration rates) in a complex demography model under a unified framework.


Installation
------------

download the package and install from local:

    git clone https://github.com/Ephraim-usc/glike.git
    
    python3 -m pip install ./glike

There are three pre-requisite packages: `tskit`, `numpy` and `scipy`.
Because glike contains a C extension module, a proper C environment is required.


Tutorial
------------

In the following sections, we introduce the API of the gLike model.

For an intuitive understanding of how gLike works, see a simple [tutorial](./tutorial.md) that walks though a toy example.


Full likelihood of genealogical trees
------------

Core functionality

    logp = glike.glike_trees(trees, demo, samples, prune)
    
Where `trees` is any enumerable that contains `tskit` genealogical trees.
Note that it is the user's duty to manually pick out trees that are selective neutral and independent, and to wrap them in a list or other iterable objects.
A `tskit.TreeSequence` object is not directly iterable and thus not a legitimate input. 
Although calling `ts.trees()` to enumerate all trees in the `tskit.TreeSequence` would work grammatically, it is not recommended in most cases, since neighboring trees are generally not independent from each other. 

`demo` is the hypothesized Demography created manually (see the following section) or from provided models in `models.py`.

`samples` is the dict that contains `sample:pop` pairs, which specifies which sample is collected from which population.
The population can be identified either by an integer representing the population index or a string denoting the population name.
Samples not mentioned in the dictionary are considered potentially from any available population at the time of the sample.
The default is an empty dictionary.

`prune` is the proportion of discarding low likelihood trees, enabling this feature often reduces noise when dealing with reconstructed trees.
The default is zero, meaning that all trees are reserved. 

This function returns the log probability that such genealogical trees are generated under the hypothesized demography.


Demography customization
------------

A Demography object is initialized with

    demo = glike.Demo()
    
And a number of Phases are added into this Demography

    demo.add_phase(phase)

A Phase is created by

    phase = glike.Phase(t, t_end, ns, grs, P, Q, populations)

Where `t` is the starting (most recent) time in generations. `t_end` is the ending (most ancient) time in generations. It is required that `t_end` is greater than `t`. 

`ns` is the 1D array of coalescent rates at the start of the phase (at time `t`). 

`grs` is the 1D array of growth rates (a positive growth rate means the population size is growing forward in time, or equivalently, the coalescent rate is growing backward in time). The default is a vector of all zeros, meaning all populations have constant sizes during this phase.

`P` is the 2D array of migration matrix at the beginning of this phase, where `P[i, j]` means the probability of a lineage from the `i`-th population (in the previous phase) to migrate into the `j`-th population in the current phase. The default is an identity matrix, meaning that no migration happens between population. Note that `P` must be provided if this phase has a different number of popuations than the previous phase.

`populations` is the 1D array of population names, this is only useful for visualizing the demography or conveniently assigning population identities to samples. If not provided, alphabatical names (A, B, C, ...) would be automatically given to populations. 

Only `t`, `t_end` and `ns` are required , other arguments are optional. The number of populations in this phase should be consistent among parameters, so it is required that

    len(ns) == len(grs) == P.shape[1] == len(populations)

When adding new Phases into Demogrpahy, the times and dimensions should match. Specifically, if the last added phase is `phases[i]`, and we are trying to add another phase `phases[i+1]`, then it is required that 

    phases[i].t_end == phases[i+1].t
    phases[i].P.shape[1] == phases[i+1].P.shape[0]

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

Where `names` is a list of the names of the parameters, `values` is a list of initial parameter values, `fixed` is a list of the names of fixed parameters so that their values will keep untouched. `limits` is a list of tuples `(low, high)` where `low` and `high` could be a number or a string expression involving names of the parameters (which will be evaluated at runtime by `eval()`). For example, if our model has three events happening between 0 and 100 generations ago, in order to estimate the times of these three events we could do

    names = ["t1", "t2", "t3"]
    values = [25, 50, 75] # or other initial values you see fit
    limits = [(0, "t2"), ("t1", "t3"), ("t2", 100)]

Empirically, if the output parameters take values very close to the lower or upper limits, it's likely that the estimation is stuck in a local optimal, or the proposed model is not well compatible with the genealogical trees. If that's the case, it's suggested to try other initial values or demography models.


About this project
-------------

This project is introduced [here](https://www.biorxiv.org/content/10.1101/2023.10.10.561787v1).

If you are interested in knowing more, please let us know. Email the author: fcq1116@gmail.com
