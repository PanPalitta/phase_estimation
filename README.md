# Evolutionary algorithm for adaptive phase estimation

We implement evolutionary algorithm to the problem of adaptive phase estimation, which is an example of quantum control problems. The aim of this project is to create a library containing modules that streamlines the construction of an optimization algorithm for quantum control problems. Access to modules of optimization algorithms provides the building blocks that users can use to tweak the algorithm to their needs.

Features:

 * Library in C++
 * Support MPI
 * Support VSL and GPU for random number generation
 * Include modules for particle swarm optimization (PSO) and differential evolution (DE)
 * Include uniform and clustered method of initializing solution candidates
 * Include access to user specified accept-reject criteria
 * Include preliminary support for multi-objective calculation
 * Support Autotools

### Links

 * Latest release: 1.0.1 [download](https://github.com/PanPalitta/phase_estimation/releases)
 * Feedback and issue: [link](https://github.com/PanPalitta/phase_estimation/issues)

## Copyright and license
This is a free software made available under the [GNU GENERAL PUBLIC LICENSE](http://www.gnu.org/licenses/gpl-3.0.html), which means you can share, modify, and redistribute this software. While we endeavor to make this software as useful and as error-free as possible, we cannot make any such guarantee, and the software is hence released **without any warranty**.


## Acknowledgement

This software has been developed by [Pantita Palittapongarnpim](https://github.com/PanPalitta) and [Peter Wittek](https://github.com/peterwittek) with the financial support of [NSERC](http://www.nserc-crsng.gc.ca/index_eng.asp) and [AITF](http://www.albertatechfutures.ca/).

The computational work was enabled by support from [WestGrid](https://www.westgrid.ca/) and [Calcul Quebec](http://www.calculquebec.ca/en/) through [Compute Canada](https://www.computecanada.ca/).


## References

1. Pantita Palittapongarnpim, Peter Wittek and Barry C. Sanders. Controlling adaptive quantum phase estimation with scalable reinforcement learning. In *Proc. 24th European Symposium on Artificial Neural Networks, Computational Intelligence and Machine Learning* (ESANN 2016): 327--332, Apr 2016.

2. Pantita Palittapongarnpim, Peter Wittek and Barry C. Sanders. Single-shot adaptive measurement for quantum-enhanced metrology. In *Proc. of SPIE Quantum Communications and Quantum Imaging XIV*, **9980** :99800H, Sep 2016. [DOI](https://doi.org/10.1117/12.2237355) [arXiv:1608.06238](https://arxiv.org/abs/1608.06238)
  
3. Pantita Palittapongarnpim, Peter Wittek, Ehsan Zahedinejad, Shakib Vedaie and Barry C. Sanders. Learning in quantum control: High-dimensional global optimization for noisy quantum dynamics, *Neurocomputing* **268**: 116--126, Apr 2017. [DOI](https://doi.org/10.1016/j.neucom.2016.12.087) [arXiv:1607.03428](https://arxiv.org/abs/1607.03428)

4. Pantita Palittapongarnpim, Peter Wittek and Barry C. Sanders. Robustness of learning-assisted adaptive quantum-enhanced metrology in the presence of noise. In *Proc. 2017 IEEE International Conference on Systems, Man and Cybernetics* (2017 SMC): 294--299, Dec 2017. [DOI](https://doi.org/10.1109/SMC.2017.8122618)
 
5. Pantita Palittapongranpim and Barry C. Sanders, Robustness of Adaptive Quantum-Enhanced Phase Estimation. 2018. [arXiv:1809.05525](https://arxiv.org/abs/1809.05525)