function [ESD,g_unif,cents,cOrds,opts_cOrds] = hb_get_denoised_esd(Y,G,kernelSet,N_unif,cOrds,opts_cOrds)
% HB_REDUCE_NOISE estimates additive Gaussian noise present in the signals,
% substracts the associate energy spectral density (ESD) of the noise from
% the that of the signal, and return the resulting noise reduced ESD. This
% is done by first estimaging the energy spectral density of the given
% signal as well as a reference Gaussian signal and then matching a portion
% of the spectrum of the reference noise signal to that of the input
% signal. The portion of the spectrum selected should be spectral regions
% in which no notable signal energy content is expected, thus, assumed to
% be just related to noise. This method is valid if one can make the
% assumption that at the specifed spectral range only/mainly Gaussian noise
% is present, and that there is no signal component or other notable noise
% type present.
%
% Hamid Behjat

if ~exist('N_unif','var') || isempty(N_unif)
    N_unif = 50;
end
if ~exist('cOrds','var') || isempty(cOrds)
    cOrds = [];
end
if ~exist('opts_cOrds','var') || isempty(opts_cOrds)
    opts_cOrds = [];
end

[N_g,N_s] = size(Y);

ESD = struct;
ESD.signal_noisy    = zeros(N_unif,N_s); % ESD of input signals
ESD.signal_denoised = zeros(N_unif,N_s); % ESD of denoised signals
ESD.noise_estimated = zeros(N_unif,N_s); % ESD of estimated noise
ESD.noise_reference = zeros(N_unif,1);   % ESD of reference noise

% reference noise
rn = randn(N_g,1); 
[rn_E,~,g_unif,cents,cOrds,opts_cOrds] = hb_get_coarse_eesd(G,rn,N_unif,[],cOrds,opts_cOrds); % ESD of reference noise
        
% denoise each signal seperately
for iY=1:N_s
    
    y = Y(:,iY);
    
    y_E = hb_get_coarse_eesd(G,y,N_unif,[],cOrds,opts_cOrds);  % energy spectral density 
    
    d = max(y_E)/max(rn_E); % approximate match
    
    kk = 0.75*d:0.01:1.25*d; % scaling search range

    e = zeros(length(kk),1); % matching errors

    % search for best match 
    for iK = 1:length(kk)
        k = kk(iK);
        d = (k*rn_E)-y_E; 
        d = d(kernelSet);
        e(iK) = norm(d);
    end
    [~,I] = min(e);

    % best noise ESD match
    d = kk(I)*rn_E;
    ESD.noise_estimated(:,iY) = d;
    
    % denoised signal ESD
    ESD.signal_denoised(:,iY) = y_E-d;
end

% input signal ESD
ESD.signal_noisy = ESD.noise_estimated + ESD.signal_denoised;

ESD.noise_reference = rn_E;
end