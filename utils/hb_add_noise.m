function [Y,N] = hb_add_noise(X,SNR)
% HB_ADD_NOIE adds additive Gaussian noise to signal to obtain a noisy
% signal at a given snr.

[N_g,N_s] = size(X);

N = zeros(size(X));

for iX=1:N_s
    x = X(:,iX);
    n_std = (1/SNR)*std(x);
    N(:,iX) = n_std*randn(N_g,1);
end
Y = X+N;
end