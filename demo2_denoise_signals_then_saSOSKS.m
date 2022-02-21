%-saSOSKS: signal-adapted system of spectral kernels. 

% Given a noisy signal set, assuming that you know the noise type, or there
% is a reason to believe that a certain noise type is dominantly present in
% the data aside form other potentail sources of noise, the signals can be
% initially denoised, and consequently, saSOSKS be buit. In doing so, the
% resulting saSOSKS will be better tailored to the signal component in the
% data, and thus, will result in decompositions of the signals that are
% more useful.

close all
clear 
clc

%-User settings.
%--------------------------------------------------------------------------
N_sa = 12;  % number of kernels in saSOSKS
SNR  = 1/8; % defining noise level based on SNR

%-Optional user settings.
%--------------------------------------------------------------------------
N_unif   = []; % number of uniform kernels for estimating graph signal energies.
kernels  = []; % kernels from the uniform SOSKS to be used for estimating energy spectral density of Gaussian noise; a set of values from within 1:N_unif.
L_smooth = []; % length of moving average window for smoothing the warping; rule of thumb: 0.05-0.1*(graph size)

%-Stuff.
%--------------------------------------------------------------------------
addpath(genpath('utils'));
d = load('SampleGraph.mat');     
G = d.G;
d = load('SampleGraphSignals_1.mat'); % S1: clean data
X = d.S; 
if isempty(N_unif)
    N_unif = 50;
end
if isempty(kernels)
   kernels= ceil(N_unif/2):N_unif;
end
if isempty(L_smooth)
    L_smooth = 0.05*G.N;
end
assert(all(kernels<=N_unif),'Kernels should be less than N_unif');

%-Add noise.
%--------------------------------------------------------------------------
[Y,N] = hb_add_noise(X,SNR); % Y=X+N

%-Get noise reduced Energy Spectral Density (ESD).
%--------------------------------------------------------------------------
[ESD,g_unif,cents,cOrds,opts] = hb_get_denoised_esd(Y,G,kernels,N_unif);

%-Ensemble ESDs.
%--------------------------------------------------------------------------
EESD = struct;

d = hb_get_coarse_eesd(G,X,N_unif,[],cOrds,opts);
EESD.signal_groundtruth = d/sum(d);

d = hb_get_coarse_eesd(G,N,N_unif,[],cOrds,opts);
EESD.noise_added = d/sum(d);

d = mean(ESD.signal_noisy,2);
EESD.signal_noisy = d/sum(d);

d = mean(ESD.signal_denoised,2);
EESD.signal_denoised = d/sum(d);

d = mean(ESD.noise_estimated,2);
EESD.noise_estimated = d/sum(d);

%-Energy equalizing (EE) warpings based on EESDs. 
%--------------------------------------------------------------------------
% for noisy siganls
d = EESD.signal_noisy;
E_noisy = d./sum(d);
E_noisy = E_noisy(:)';
[w_noisy,~,w_ns_noisy] = hb_get_coarse_ee_warping(G,E_noisy,cents,L_smooth);

% for noise-reduced signals
d = abs(EESD.signal_denoised);
E = d./sum(d);
E = E(:)';
[w,e,w_ns,G] = hb_get_coarse_ee_warping(G,E,cents,L_smooth);

%-Signal-adapted SOSKS (saSOSKS).
%--------------------------------------------------------------------------
% based on on noisy data
g_noisy = spgg_filter_design(G.lmax,N_sa,...
    'designtype','signal_adapted','warping',w_noisy,'E',e);
% based on noise-reduced data
g = spgg_filter_design(G.lmax,N_sa,...
    'designtype','signal_adapted','warping',w,'E',e);


%-Plots: resulting saSOSKS etc.
%--------------------------------------------------------------------------

% uniform SOSKS (used for ESD estimations)
hf1 = figure(1);
clf(1);
set(hf1,'position',[620 50 600 200]);
hold on;
for k=1:N_unif
    plot(e,g_unif{k}(e));
end
title('Uniform SOSKS used for estimating ESD of graph signals')
xlabel('\lambda');
ylabel('amplitude');
xlim([0 G.lmax]);
ylim([0 1.05]);
grid on
box off

% estimated EESDs
hf2 = figure(2);
clf(2);
set(hf2,'position',[620 330 1000 460]);
hold on;
bar(cents,[EESD.signal_groundtruth(:) E_noisy(:) E(:)]);
title('Coarse Ensemble Energy Spectral Density (EESD) of graph signals')
legend(...
    'EESD of clean signals', ...
    'EESD of noisy signals', ...
    'EESD of noise reduced signals');
xlabel('\lambda');
ylabel('% of total signal power');
xlim([0 G.lmax]);
grid on
box off

% EE warpings
hf3 = figure(3);
clf(3);
set(hf3,'position',[1 50 600 200]);
hold on;
plot(e,e/G.lmax,'r:','displayname','no warping');
plot(e,w_ns,'m',...
    'displayname','non-smoothed warping (based on noise-reduced signals)');
plot(e,w,'k',...
    'displayname','smoothed warping (based on noise-reduced signals)');
plot(e,w_noisy,'b',...
    'displayname','smoothed warping (based on noisy signals)');
title('Warping functions');
xlabel('\lambda');
ylabel('warping');
xlim([0 G.lmax])
ylim([0 1.05])
legend('location','se');
grid on
box off

% saSOSKS (based on noisy data)
hf4 = figure(4);
clf(4);
set(hf4,'position',[1 330 600 200]);
axis;
ha = gca;
N_interp = 1000;
d = 0:1/N_interp:G.lmax;
spgg_view_design(g_noisy,[0,G.lmax],...
    'Graph',G,'eigsToInterp',d,'plotLineWidth',2,'guiHandle',ha);
title('saSOSKS based on energy spectral density of noisy signals')
xlabel('\lambda');
ylabel('amplitude');
xlim([0 G.lmax])
ylim([0 1.05])
grid on
box off

% saSOSKS (baed on noise-reduced data)
hf5 = figure(5);
clf(5);
set(hf5,'position',[1 600 600 200]);
axis;
ha = gca;
N_interp = 1000;
d = 0:1/N_interp:G.lmax;
spgg_view_design(g,[0,G.lmax],...
    'Graph',G,'eigsToInterp',d,'plotLineWidth',2,'guiHandle',ha);
title('saSOSKS based on energy spectral density of noise-reduced signals');
xlabel('\lambda');
ylabel('amplitude');
xlim([0 G.lmax])
ylim([0 1.05])
grid on
box off

% additional plots related to noise reduction 
showFigs = false; % true
if showFigs
    hf6 = figure(6); %#ok<UNRCH>
    clf(6);
    set(hf6, 'position',[1 100 1000 460])
    xx = cents;
    subplot(611)
    bar(xx,EESD.signal_groundtruth);
    title('signal (S)');
    
    subplot(612)
    bar(xx,EESD.noise_added);
    title('noise (N) [randn, std = std(S)/SNR]');
    
    subplot(613)
    bar(xx,EESD.signal_noisy);
    title(sprintf('noisy signal (NS): N+S [SNR = %d]',SNR));
    
    subplot(614)
    bar(xx,ESD.noise_reference);
    title('reference noise [randn, std = 1]');
    
    subplot(615)
    bar(xx,EESD.noise_estimated);
    title('estimated noise (EN)');
    
    subplot(616)
    d1 = EESD.signal_denoised;
    d2 = abs(EESD.signal_denoised);
    d3 = EESD.signal_groundtruth*sum(d2);
    bar(xx,[d1(:) d2(:) d3(:)]);
    title('cleaned signal: NS-EN');
    legend('cleaned signal (CS)','|CS|','S');
end
