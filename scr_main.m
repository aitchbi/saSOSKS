%---saSOSKS: signal-adapted system of spectral kernels

close all
clear 
clc

%-User settings.
%--------------------------------------------------------------------------
N_sa     = 9;    % # of signal-adapted kernels

%-Optional user settings.
%--------------------------------------------------------------------------
N_unif   = 50;   % # of uniform kernels for estimating graph signal energies
L_smooth = 2500; % length of moving average window for smoothing the warping
UnifKernelType = 'MeyerLike'; % type of uniform kernels for estimating energy

%-Stuff.
%--------------------------------------------------------------------------
addpath(genpath('utils'));
d = load('SampleGraph.mat');     
G = d.G;
d = load('SampleGraphSignals_4.mat');
S = d.S; 

%-Load precomputed Chebyshev orders if available. 
%--------------------------------------------------------------------------
[cOrds,opts_cOrds] = hb_get_precomp_cheby_ords(UnifKernelType,N_unif);

%-Build energy equalizing spectral warping.
%--------------------------------------------------------------------------
[w,e,w_ns,G,g_unif,E,cents] = hb_ee_warping(...
    G,S,L_smooth,N_unif,UnifKernelType,cOrds,opts_cOrds);

%-Build signal-adapted SOSKS.
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
