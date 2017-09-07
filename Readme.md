# gretaDiscrete: module for sampling from discrete random variables in greta

The broad purpose of this module is to allow discrete random variables in greta models. A particular focus of this module is the development of variable selection tools in greta, using either stochastic search variable selection (SSVS) or a transdimensional sampler (e.g. reversible jump MCMC). Continuous random variables will be updated using greta's Hamiltonian Monte Carlo (HMC) samplers and discrete variables will be updated using an adaptive MCMC step outside of the HMC sampler. 

