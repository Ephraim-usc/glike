gLike
========

genealogical Likelihood (gLike) is a maximum likelihood method to infer the demographic history
of given populations that best explains the observed genealogical relationships between sample haplotypes. 
This is, to our knowledge, the first attempt to infer all types of parameters
(split times, admixture times and proportions, migration rates) in a complex demography model under a unified framework.


Installation
------------

download the package and install from local:

    git clone https://github.com/Ephraim-usc/glike.git
    
    pip3 install ./glike


Full likelihood of genealogical trees
------------

Core functionality

    glike_trees(trees, demo)
    
Where `trees` is any enumerable that contains `tskit` genealogical trees.
And `demo` is the hypothesized demography created manually or from provided models in `models.py`.
It returns the probability that such genealogical trees are generated under hypothesized demography.


About this project
-------------

This is an ongoing project, please refer to our poster at ASHG 2022 for a brief introduction.

If you are interested in knowing more, please let us know. Email the author: caoqifan@usc.edu

![](images/poster_ashg.png)
