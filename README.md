# saSOSKS
Design of signal-adapted system of spectral graph kernels (saSOSKS), via coarse estimation of the enemble energy spectral density (EESD) of a given set of graph signals defined on a given graph.

`demo1` shows how saSOSKS can be constructed based on a coarse estimate of the EESD of the given graph signal set. The approach is particularly suitable for large graph for which the graph Fourier transform (GFT) of graph signals cannot be directly computed due to the sheer size of the graphs that renders computing the eigenvalue decomposition of the graph Laplacian matrix infusible. 

`demo2` shows how potential noise in the given signal set can be first estimated and subtracted from the signals to enable designing saSOSKS that are better tailored to the energy spectral density of the desired signal, not that of noise. 
