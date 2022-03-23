# saSOSKS
Design of signal-adapted system of spectral graph kernels (saSOSKS), via coarse estimation of the ensemble energy spectral density (EESD) of a given set of signals defined on a graph. The implementation is based on the theory presented in this [paper](https://bme.lth.se/fileadmin/biomedicalengineering/Personal_folders/Hamid_Behjat/HBehjat_chapter2018.pdf):

> Behjat, H. and Van De Ville, D., 2019. Spectral Design of Signal-Adapted Tight Frames on Graphs. Vertex-Frequency Analysis of Graph Signals. Springer, Cham. [DOI](https://doi.org/10.1007/978-3-030-03574-7_4)

The original theory for designing signal-adapted tight frames that does not use the coarse estimation of EESD as presented here was presented in this [paper](https://doi.org/10.1109/TSP.2016.2591513):

> Behjat et al., 2016. Signal-adapted tight frames on graphs. IEEE Trans. Signal Process. [DOI](https://doi.org/10.1109/TSP.2016.2591513)

and its implementation is given on a seperate repository, [spg](https://github.com/aitchbi/spg) [the two repositories are to be merged in future].

`demo1` shows how saSOSKS can be constructed based on a coarse estimate of the EESD of the given graph signal set. The approach is particularly suitable for large graph for which the graph Fourier transform (GFT) of graph signals cannot be directly computed due to the sheer size of the graphs that makes it infeasible to compute the full eige-decomposition of the graph Laplacian matrix. An example 9-kernel saSOSKS and the associated energy-equalizing warping function that is used for generating the saSOSKS via warping a uniform SOSKS is shown below: 

![Sample energy equalizing warping function used to build saSOSKS](figs/demo1/warping_SampleGraphSignalSet4.jpg?raw=true)

![sample saSOSKS built on noisy signals](figs/demo1/saSOSKS_N09_SampleGraphSignalSet4.jpg?raw=true)

More examples can be found in folder `figs/demo1`. 

`demo2` shows how potential noise in the given signal set can be first estimated and subtracted from the signals to enable designing saSOSKS that are better tailored to EESD of the desired signals, not that of noise. For example, if you are given a noisy signal set, the resulting saSOSKS can look like the following if you just use the signals as they are:

![sample saSOSKS built on noisy signals](figs/demo2/saSOSKS_N12_SNR1over4_based_on_noisy_signals.jpg?raw=true)

However, if you reduce the noise level and then build saSOSKS, the result will look like this: 

![sample saSOSKS built on noise reduced signals](figs/demo2/saSOSKS_N12_SNR1over4_based_on_noise_reduced_signals.jpg?raw=true)

More examples can be found in folder `figs/demo2`.

 
