# saSOSKS
Design of signal-adapted system of spectral graph kernels (saSOSKS), via coarse estimation of the ensemble energy spectral density (EESD) of a given set of signals defined on a graph. The implementation is based on the theory presented in this [paper](https://bme.lth.se/fileadmin/biomedicalengineering/Personal_folders/Hamid_Behjat/HBehjat_chapter2018.pdf):

> Behjat, H., Van De Ville, D., 2019. Spectral Design of Signal-Adapted Tight Frames on Graphs. Vertex-Frequency Analysis of Graph Signals. Springer, Cham. [DOI](https://doi.org/10.1007/978-3-030-03574-7_4)

`demo1` shows how saSOSKS can be constructed based on a coarse estimate of the EESD of the given graph signal set. The approach is particularly suitable for large graph for which the graph Fourier transform (GFT) of graph signals cannot be directly computed due to the sheer size of the graphs that makes it infeasible to compute the full eige-decomposition of the graph Laplacian matrix. 

`demo2` shows how potential noise in the given signal set can be first estimated and subtracted from the signals to enable designing saSOSKS that are better tailored to EESD of the desired signals, not that of noise. 
