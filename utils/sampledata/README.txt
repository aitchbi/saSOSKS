The Graph represent the cerebral cortex of a single hemisphere of a healthy subject from the Human Connectome Project database. Graph nodes represent voxels represeting the cerebral cortex at the resolution of 2 mm cubic, and graph edges represent the 26-neighborhood connectivity of adjacent voxels.    

The sample graph signals were defined as follows: 
- S1: addition of three eigenvectors of the graph's normalized Laplacian matrix; eigenvectors corresponding to eigenvalues around 0.25, 0.5 and 0.75. 
- S2: G.A*S1
- S3: (G.A)^2*S1
- S4: S1+n, where n = 0.1*mean(std(S1))*randn(size(S1));
- S5: S1+n, where n = mean(std(S1))*randn(size(S1));


