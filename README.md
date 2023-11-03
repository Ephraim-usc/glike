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

See a simple [tutorial](./tutorial.md) to go though a toy example showing how glike works.


Full likelihood of genealogical trees
------------

Core functionality

    logp = glike.glike_trees(trees, demo, samples, prune)
    
Where `trees` is any enumerable that contains `tskit` genealogical trees.
It should be noted that a `tskit.TreeSequence` object is not directly iterable, and it's the user's duty to manually pick out trees that are selective neutral and mutually independent, and to wrap them in a list or other iterable objects. 
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

    phase = glike.Phase(t, t_end, ns, grs, P, Q, populations)

Where `t` is the starting (most recent) time, `t_end` is the ending (most ancient) time, `ns` is the vector of coalescent rates, `grs` is the vector of growth rates (a positive growth rate means the population size is growing forward in time, or equivalently, the coalescent rate is growing backward in time), `P` is the mass migration matrix at the beginning of this phase, `Q` is the continuous migration rate matrix, and `populations` is the vector of population names. Only `t`, `t_end` and `ns` are required , other arguments are optional. The number of populations in this phase should be consistent among parameters, so it is required that

    len(ns) == len(grs) == P.shape[1] == Q.shape[1] == len(populations)

When adding new Phases into Demogrpahy, the times and dimensions should match. Specifically, if the last added phase is `phases[i]`, and we are trying to add another phase `phases[i+1]`, then it is required that

    phases[i].t < phases[i].t_end = phases[i+1].t < phases[i+1].t_end
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

This project is introduced [here](https://www.biorxiv.org/content/10.1101/2023.10.10.561787v1).

If you are interested in knowing more, please let us know. Email the author: caoqifan@usc.edu
