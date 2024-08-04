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

Find in [paper scripts](./paper/) codes for reproducing the figures in the paper.

Full likelihood of genealogical trees
------------

Core functionality

    logp = glike.glike_trees(trees, demo, samples, kappa, prune)
    
Where `trees` is any enumerable that contains `tskit` genealogical trees.
Note that it is the user's duty to manually pick out trees that are selective neutral and independent, and to wrap them in a list or other iterable objects.
A `tskit.TreeSequence` object is not directly iterable and thus not a legitimate input. 
Although calling `ts.trees()` to enumerate all trees in the `tskit.TreeSequence` would work grammatically, it is not recommended in most cases, since neighboring trees are generally not independent from each other. 

`demo` is a `glike.Demo` object, representing the hypothesized demography. It may be created manually (see the following section) or from provided models in `models.py`.

`samples` is a `dict` object that contains `sample:pop` pairs, which specifies which sample is collected from which population.
The population can be identified either by an integer representing the population index or a string denoting the population name.
For example, `{4:2, 13:"A"}` means that lineage 4 is in the second population (at time `tree.time(4)`), and lineage 13 is in the population named A (at time `tree.time(13)`).
In essense, the parameter `samples` restricts the graph of states (GOS) to contain only states that are compatible with these known population identities.
Although rarely useful, inner nodes in the genealogical trees may also be included in `samples`.
Lineages not mentioned in the dictionary are considered potentially from any available population at the time of the sample.
The default is an empty dictionary.

`kappa` is an integer that controls the maximum number of connections between layers (see the paper for a more technical explanation). 
It effectively controls the trade-off between accuracy and computational cost.
The default value is 10000. Experiments in the paper were conducted with this default value.

`prune` is a float number between 0 and 1, that specifies the proportion of discarding low likelihood trees.
Enabling this feature often reduces noise when dealing with reconstructed trees.
The default is 0.0, meaning that all trees are preserved. 

This function returns the log probability that such genealogical trees are generated under the hypothesized demography.


Demography customization
------------

Different from population-based languages to describe a demography (such as in `msprime`), gLike represents a demography as a series of phases (or intervals) going backward in time.
These phases are gapless in between, ensuring they collectively represent a complete history.
Each phase contains a certain number of populations, which can vary between phases.
At the transition from one phase to the next, a migration matrix defines the probabilities of lineages moving from one population to another.

A demography object is initialized with

    demo = glike.Demo()

And a number of phases are added into this Demography

    demo.add_phase(phase)

A phase is created by

    phase = glike.Phase(t, t_end, ns, grs, P, populations)

Where `t` is the starting (most recent) time in generations. `t_end` is the ending (most ancient) time in generations. It is required that `t_end` is greater than `t`. 

`ns` is the 1D array of coalescent rates at the start of the phase (at time `t`). 

`grs` is the 1D array of growth rates (a positive growth rate means the population size is growing forward in time, or equivalently, the coalescent rate is growing backward in time). The default is a vector of all zeros, meaning all populations have constant sizes during this phase.

`P` is the 2D array of migration matrix at the beginning of this phase, where `P[i, j]` means the probability of a lineage from the `i`-th population (in the previous phase) to migrate into the `j`-th population in the current phase. The default is an identity matrix, meaning that no migration happens between population. Note that `P` must be provided if this phase has a different number of popuations than the previous phase.

`populations` is the 1D array of population names, this is only useful for visualizing the demography or conveniently assigning population identities to samples. If not provided, alphabatical names (A, B, C, ...) would be automatically given to populations. 

Only `t`, `t_end` and `ns` are required , other arguments are optional. The number of populations in this phase should be consistent among parameters, so it is required that

    len(ns) == len(grs) == P.shape[1] == len(populations)

When adding new phases into a demogrpahy, the times and dimensions should match. Specifically, if the last added phase is `phases[i]`, and we are trying to add another phase `phases[i+1]`, then it is required that 

    phases[i].t_end == phases[i+1].t
    phases[i].P.shape[1] == phases[i+1].P.shape[0]

We recommend adding phases until `np.inf` for safety, although for gLike's purpose it suffices if the demography covers the oldest lineage in the input trees. 

The resulting demography can be visualized by

    demo.print()


Demography parameter estimation
------------

gLike provides a function `maximize` that mimics the popular `scipy.optimize.minimize` function, but has been made convenient for optimizing demographic parameters.

    x, logp = glike.maximize(fun, x0, bounds = bounds)

`fun` is the objective funciton that is called in the `fun(**x)` way, note that `x` must contain all parameters that `fun` does not have a default value.

`x0` is a `dict` object that contains the intial values represented as `param:value` pairs, where `param` is the name of the parameter (of type `str`), and `value` is the value (of type `float`).

`bounds` is a list of 2-tuples, of the same length as `x0`, where each tuple contains the lower and upper bound of the parameter. Specifically, it is a list of tuples `(low, high)` where `low` and `high` could be a number or a string expression including names of the parameters (which will be evaluated at runtime by `eval()`). For example, 

    x0 = {"t1":25, "t2":50, "t3":75]
    bounds = [(0, "t2"), ("t1", "t3"), ("t2", 100)]

Means that three parameters are being estimated, with names `t1`, `t2` and `t3`. Their initial values are `25`, `50` and `75`, respectively. 
The boundary conditions require that `0 < t1 < t2 < t3 < 100`.
Empirically, if the output parameters take values very close to the lower or upper limits, it's likely that the estimation is stuck in a local optimal, or the proposed model is not well compatible with the genealogical trees. If that's the case, it's suggested to try other initial values or demography models.

The wrapper function `fun` for estimating demographic parameters can usually be constructed like:

    def func(...):
      demo = glike.Demography()
      # add Phases that depend on the parameters
      # define samples if needed
      return glike.glike_trees(trees, demo, samples, kappa, prune)


Runtime considerations
------------

The computational time and memory cost for a gLike evaluation of a tree depends on the number of states and the number of connections between states, or, in the language of graph theory, the number of vertices and edges of the GOS. We note that there is not an intuitive way to predict the runtime of gLike on different demographic models (see the paper for a detailed discussion). In our experience, the following guidelines may be helpful:

(1) always simulate data and test gLike on true parameters to make sure the runtime is acceptable

(2) moving demographic events to the more recent times generally increases the number of states, because more lineages will be involved

(3) number of ancestral populations increase the number of states in an exponential manner

(4) samples from ancestral populations help determine population identities of lineages from admixed samples, therefore effectively reducing the number of states.

(5) adjust the parameter `kappa` to achieve a balance between accuracy and computational cost.


About this project
-------------

This project is introduced [here](https://www.biorxiv.org/content/10.1101/2023.10.10.561787v1).

If you are interested in knowing more, please let us know. Email the author: fcq1116@gmail.com
