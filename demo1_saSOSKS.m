%---saSOSKS: signal-adapted system of spectral kernels

% Given a signal set, S, a saSOSKS is constructred in which the kernels are
% distributed along the spectrum in such way that they each capture, on
% average, an equal amount of spectral energy based on the ensemble energy
% spectral density of S. The estimation of EESD is done in acoare way using
% estimations of energy content at different parts of the spectrum using a
% uniform SOSKS. The approach is particularly suitable for large graph for
% which the graph Fourier transform (GFT) of graph signals cannot be
% directly computed due to the sheer size of the graphs that renders
% computing the eigenvalue decomposition of the graph Laplacian matrix
% infusible.

close all
clear 
clc

%-User settings.
%--------------------------------------------------------------------------
N_sa     = 7; % number of signal-adapted kernels

%-Optional user settings.
%--------------------------------------------------------------------------
N_unif   = []; % number of uniform kernels for estimating graph signal energies
L_smooth = []; % length of moving average window for smoothing the warping; rule of thumb: 0.05-0.1*(graph size)

%-Stuff.
%--------------------------------------------------------------------------
addpath(genpath('utils'));
d = load('SampleGraph.mat');     
G = d.G;
d = load('SampleGraphSignals_4.mat'); % there are 4 sample signal sets to check
S = d.S; 
if isempty(N_unif)
    N_unif = 50;
end
if isempty(L_smooth)
    L_smooth = 0.05*G.N;
end

%-Coarse estimate of the ensemble ESD of S.
%--------------------------------------------------------------------------
[E,G,g_unif,cents] = hb_get_coarse_eesd(G,S,N_unif);
E = E./sum(E);
E = E(:)';

%-Energy equalizing (EE) warping based on EESD of S. 
%--------------------------------------------------------------------------
[w,e,w_ns,G] = hb_get_coarse_ee_warping(G,E,cents,L_smooth);

%-saSOSKS.
%--------------------------------------------------------------------------
g = spgg_filter_design(G.lmax,N_sa,...
    'designtype','signal_adapted','warping',w,'E',e);

%-Plots.
%--------------------------------------------------------------------------
% uniform SOSKS (used for signal energy estimation)
hf1 = figure(1);
clf(1);
set(hf1,'position',[650 50 600 200]);
hold on;
for k=1:N_unif
    plot(e,g_unif{k}(e));
end
xlabel('\lambda');
ylabel('amplitude');
xlim([0 G.lmax]);
ylim([0 1.05]);
grid on
box off

% estimated ensemble energy spectral density
hf2 = figure(2);
clf(2);
set(hf2,'position',[650 330 600 200]);
hold on;
bar(cents,E,'k');
xlabel('\lambda');
ylabel('% of total signal power');
xlim([0 G.lmax]);
grid on
box off

% energy equalizing warping function
hf3 = figure(3);
clf(3);
set(hf3,'position',[1 330 600 200]);
hold on;
plot(e,e/G.lmax,'r:','displayname','no warping');
plot(e,w_ns,'m','displayname','non-smoothed warping');
plot(e,w,'k','displayname','smoothed warping');
xlabel('\lambda');
ylabel('warping');
xlim([0 G.lmax])
ylim([0 1.05])
legend('location','se');
grid on
box off

% saSOSKS
hf4 = figure(4);
clf(4);
set(hf4,'position',[1 50 600 200]);
axis;
ha = gca;
N_interp = 1000;
d = 0:1/N_interp:G.lmax;
spgg_view_design(g,[0,G.lmax],...
    'Graph',G,'eigsToInterp',d,'plotLineWidth',2,'guiHandle',ha);
xlabel('\lambda');
ylabel('amplitude');
xlim([0 G.lmax])
ylim([0 1.05])
grid on
box off