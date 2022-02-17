% This function is an extended version of spg_filter_design.m from: 
% https://github.com/aitchbi/spg, 
% which included a number of additional features as specified in the below
% change-log:
% 
% [Sep 2019]
% Included option for 'signal_adapted_spline_type'
%
% [Sep 2019]
% Included option for designing 'uniform_multires' sosks. 
% Example: 
% cd('/Users/hamid/Dropbox/gits/braingraph_working_folder/estimate_energy_spectra_based_on_chapter_method_June2019')
% addpath('./disi')
% addpath('./hb')
% addpath(genpath('/Users/hamid/Dropbox/my/my_toolboxes/spg'));
% addpath(genpath('/Users/hamid/Dropbox/my/my_toolboxes/spg_v1.01'));
% addpath('/Users/hamid/Dropbox/my/my_toolboxes/handy');
% opts.sz = 1e-4;
% opts.lmax = 2;
% opts.sOrder = 3;
% opts.arange = [0 lmax];
% opts.Nsbs_l = 400; % res: 0.005
% opts.Nsbs_u = 40; % res: 0.05
% opts.sbs_l = 1:20;
% opts.sbs_u = 4:Nsbs_u;
% [gg,g_l,g_u,cents] = spgg_filter_design([],[],...
%     'designtype','uniform_multires',...
%     'opts_uniMultRes',opts);
%
% [25 June 2017]
% added: 'vers' option. Default is 1, which gives 'uniform_meyer_type' and
% 'signal_adapted' in the standard way such that the sum of squared abs 
% values of the kernels form a partition of unity. If 'vers' 2 is used, 
% instead the sum of abs values of the kernels form a partion of unity. 
%
% [25 April 2017]
% In 'spectrum-adapted' option, modified gsp_spectrum_cdf_approx.m to 
% enable handling MatFile G input. Renamed the function to: 
% hb_gsp_spectrum_cdf_approx.m. 
% Thus, REQUIRED function: hb_gsp_spectrum_cdf_approx.m
%
% [23 April 2017]
% added option for binbased warping input; i.e. instead of specifying a
% warping that is G.N long, define one that is M long, where 
% M=numel(warpingBins), and warpingBins are the bin centers at which the
% warping is defined as. 
%
% [17 Feb 2017] 
% Note that 'subSampleWarping' is good to use when the eigenvalues are 
% somehow densly located, so if you skip every second or third or ..., 
% the resulting kernels will be smoother. 
%
% ON THE OTHER HAND, the new added option: 'fixForSparseSpectra' is good to 
% handle small graphs which have sparsly located eigenvalues. For instance, 
% the EMG correlation graph was like this: had 16 eigenvalues, lmax= 1, and 
% the distance between 1st and 2nd was 0.65!!! This leads to a 
% signal-adapted warping function that lies on the diagonal (y=x) for the 
% range [0 0.65], and this means, in this range, the Meyer-like uniform 
% translate kernels will be maintained, i.e. no warp. And this is really 
% bad, since not only there is no energy in this range, there is even no 
% eigenvalues in this part of the spectrum.      



% SPG_FILTER_DESIGN Construct system of spectral kernels for a given graph
%                               (and graph signal set).
%
% Baed on the provided inputs, SPG_FILTER_DESIGN constructs:
%
% (1) Uniform Meyer-type system of spectral kernels    
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','uniform_meyer_type','warping',warping,'E',E)
% 
% &
%
% (2) Signal-adapted system of spectral kernels
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','signal_adapted')
%
% as proposed in:
%
%                H. Behjat, et al., ''Signal-adapted tight frames on graphs'', 
%                IEEE Trans. Signal Process., 2016, 
%                doi: 10.1109/TSP.2016.2591513.
%
% (3) Spline-based system of spectral kernels (SGWT frame)
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','abspline3')
%
% as proposed in:
%
%                D.K. Hammond, et al., ''Wavelets on graphs via spectral graph
%                theory'', Appl. Comput. Harmon. Anal., vol. 30, 
%                pp. 129-150, 2011.   
%
% (4) Meyer-like system of spectral kernels 
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','meyer')
%
% as proposed in:
%
%                N. Leonardi, et al., ''Tight wavelet frames on multislice graphs'', 
%                IEEE Trans. Signal Process., vol. 61(13), pp. 3357-3367, 2013.  
%
% (5) Half-cosine unform translates system of spectral kernels
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','half_cosine_uniform_translates')
%
% &
%
% (6) Sectrum-adapted system of spectral kernels
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','spectrum_adapted','G',G)
% 
% as proposed in:
%
%                D. I. Shuman, et al., ''Spectrum-adapted tight graph wavelet
%                and vertex-frequency frames'', IEEE Trans. Signal Process., 
%                vol. 63(16), pp. 4223-4235, 2015.
%
% ------------------------------------------------------------
%
% Inputs:
%
% lmax             : upper bound on spectrum
%
% N_subbands   : number of subbands [default: 7]
%
%
% Optional input parameters :
%
% 'designtype'      : Type of system of spectral kernels to construct 
%                          [default: 'signal_adapted']
%
% 'warping'          : Energy-equalizing transformation used for warping 
%                          the uniform Meyer-type system of spectral kernels
%                          to construct a signal-adapted syetem of spectral 
%                          kernels.
%                        
% 'G'                   : graph structure containg at least the following 
%                         required fields:
%                         - N        : Number of graph vertices
%                         - L        : Graph Laplacian matrix
%                      
% 'E'                   : eigenvalues of the graph Laplacian matrix 
%
% 'subSampleWarping' : [default: 0]
%                                 If nonzero, constructs energy-equalizing transformation  
%                                 based on every N-th eigenvalue, 
%                                 where subSampleWarping = N. 
%
%
% Outputs:
%
% g{1}                 : function handle for first subband spectral kernel
% g{2}                 : function handle for second subband spectral kernel
% ...
% g{N_subbands} : function handle for last subband spectral kernel 
%
%
% ------------------------------------------------------------
% Copyright (C) 2016, Hamid Behjat.
% This file is part of SPG (signal processing on graphs) package - v1.00
%
% Download: miplab.epfl.ch/software/
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [g,varargout] = spgg_filter_design(lmax,N_subbands,varargin)
control_params = {...
    'designtype',[],...
    'warping',[],...
    'G',[],...
    'E',[],...
    'subSampleWarping',0,...
    'fixForSparseSpectra',0,... % added in v1.01 [17 Feb 2017]
    'minSparseFactor',[],... % added in v1.01 [17 Feb 2017]
    'extrapolStep',[],... % added in v1.01 [17 Feb 2017]
    'binbasedWarping',0,... % added in v1.01 [23 Feb 2017]
    'warpingBins',[],... % added in v1.01 [23 Feb 2017]
    'vers',1,... % [25 June 2017]
    'pou','over1stPower',... % [21 Jan 2018]
    'sOrder',3,... % [21 Jan 2018]
    'sz',0.001,... % [21 Jan 2018] 
    'opts_uniMultRes', [],... % Sep 2019
    'opts_lphp', []}; % Nov 2020
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
if ~isempty(N_subbands) % empty for designtype 'uniform_multires'
    N_subbands = N_subbands-1;
end
g = cell(N_subbands+1,1);

if subSampleWarping==1 && fixForSparseSpectra==1
    error('This scenario has not been implemented.')
    % I think this even does not make sense to implement. [17 Feb 2017]
end

switch designtype
    case 'meyer_lowpass_highpass'
        % Meyer-type spectral kernel pair with quadrature mirror property
        % at transition band
        d1 = opts_lphp.t; % transition point
        d2 = opts_lphp.w; % width of transition band [d2 should be < 2*d1]
        if isfield(opts_lphp,'vers')
            vers = opts_lphp.vers;
        end
        g{1} = @(x) spgg_kernel_meyer_lowpass(x,d1,d2,vers);
        g{2} = @(x) spgg_kernel_meyer_highpass(x,d1,d2,vers);
    case {'uniform_multires'}
        opts = opts_uniMultRes;
        [g,g_l,g_u] = spgg_get_multiresuniform(opts);
        
        if nargout>2
            varargout{2} = g_l;
            varargout{3} = g_u;
        end
        
    case {'uniform_spline_type','signal_adapted_spline_type'}
        
        switch designtype
            case 'uniform_spline_type'
                % pou = 'over2ndPower' or 'over1stPower'
                % can be either based on application.
                
            case 'signal_adapted_spline_type'
                pou = 'over2ndPower';
                if isempty(warping)
                    error('Warping function missing.');
                end
                if isempty(warpingBins) && isempty(E)
                    d = 'Either ''E'' or ''warpingBins'' needed as input.';
                    error(d);
                end
        end
        
        so = sOrder;
        N = N_subbands+1;
        rmax = N-1;
        sc = lmax/rmax;
        g = cell(1,N);
        
        if N<4 || lmax>(N+1)
            error('inappropriate N or lmax');
        end
        
        [xl,xh,pl,ph] = spgg_aux1_splinetype(so,sz,pou);
        for iN = 1:N
            [d3,d4] = spgg_aux2_splinetype(so,sz,iN,N,sc,lmax,pou,xl,xh,pl,ph);
            if ~isempty(warping)
                d3 = spgg_apply_warping(d3,warping,lmax,E,warpingBins,...
                    subSampleWarping,fixForSparseSpectra,...
                    minSparseFactor,extrapolStep);
            end
            g{iN} = @(k) interp1(d3,d4,k);
        end
        
    case 'uniform_meyer_type'
        % Uniform Meyer-type system of spectral kernels
        % sum of squared abs values forms a partition of unity
        
        g{1} = @(x) spgg_kernel_1(x,lmax,N_subbands,'vers',vers);
        for j = 1:N_subbands
            g{j+1} = @(x) spgg_kernel_2(x,j,lmax,N_subbands,'vers',vers);
        end
        
    case 'signal_adapted'
        % Signal-adapted system of spectral kernels
        
        if subSampleWarping
            % construct warping based on every N-th lambda;
            % where N = subSampleWarping .
            g{1} = @(x) spg_kernel_1(x,lmax,N_subbands,...
                'warping',warping,'E',E,...
                'subSampleWarping',subSampleWarping);
            
            for j = 1:N_subbands
                g{j+1} = @(x) spg_kernel_2(x,j,lmax,N_subbands,...
                    'warping',warping,'E',E,...
                    'subSampleWarping',subSampleWarping);
            end
        elseif fixForSparseSpectra
            % Account for the fact that eigenvalues might be too sparse, 
            % and appropriately adjust the design of the warping function.
            
            g{1} = @(x) spgg_kernel_1(x,lmax,N_subbands,...
                'warping',warping,'E',E,...
                'fixForSparseSpectra',fixForSparseSpectra,...
                'vers',vers);
            
            for j = 1:N_subbands
                g{j+1} = @(x) spgg_kernel_2(x,j,lmax,N_subbands,...
                    'warping',warping,'E',E,...
                    'fixForSparseSpectra',fixForSparseSpectra,...
                    'vers',vers);
            end
            
        elseif binbasedWarping
            % Enables inputing a warping function which is not specified at
            % every G.E; instead, it is defined at a set of M bins, with 
            % centers given by 'warpingBins', and M=numel(warpingBins). 
            
            if ~isempty(warpingBins)
                g{1} = @(x) spgg_kernel_1(x,lmax,N_subbands,...
                    'warping',warping,'warpingBins',warpingBins,...
                    'vers',vers);
                for j = 1:N_subbands
                    g{j+1} = @(x) spgg_kernel_2(x,j,lmax,N_subbands,...
                        'warping',warping,'warpingBins',warpingBins,...
                        'vers',vers);
                end
            else
                error('Bin centers should be specified.')
            end
        else
            g{1} = @(x) spgg_kernel_1(x,lmax,N_subbands,...
                'warping',warping,'E',E,'vers',vers);
            for j = 1:N_subbands
                g{j+1} = @(x) spgg_kernel_2(x,j,lmax,N_subbands,...
                    'warping',warping,'E',E,'vers',vers);
            end
        end
   
    case 'abspline3'
        % Spline-based system of spectral kernels (SGWT frame)
        
        g = sgwt_filter_design(lmax,N_subbands,'designtype','abspline3');
        
    case 'meyer'
        % Meyer-like system of spectral kernels 
        
        gb = @(x) sgwt_meyer(x);
        gb2 = @(x) sgwt_mey_h(x);
        gb3 = @(x) sgwt_meyer_end(x);
        
        t = sgwt_setscales_mey(lmax,N_subbands);
        
        % wavelet kernels
        for j = 1:N_subbands-1
            g{j+1} = @(x) gb(t(j)*x);
        end
        g{N_subbands+1} = @(x) gb3(t(end)*x);
        
        % sacling function kernel
        g{1} = @(x) gb2(t(1)*x);
        
        
    case 'half_cosine_uniform_translates'
        % Half-cosine unform translates system of spectral kernels
        
        g = gsp_filter_design('uniform_translates',N_subbands+1,lmax);
        
    case 'spectrum_adapted'
        % Spectrum-adapted system of spectral kernels
        
        G = hb_gsp_spectrum_cdf_approx(G);
        wParam.warp_function_type = 'custom';
        wParam.warp_function = G.spectral_warp_fn;
        wParam.upper_bound_translates = 1;
        
        g = gsp_filter_design('warped_translates',N_subbands+1,lmax,wParam);
        
    otherwise
        error('Unrecognized sosks design type.');
end

if nargout==1
    return;
end

% Kernel centers; varargout{1} --------------------------------------------
switch designtype
    case {'uniform_multires'}
        d1 = hb_get_kernel_cents(g,opts.sz,opts.lmax);
        d2 = hb_get_kernel_cents(g_l,opts.sz,opts.lmax);
        d3 = hb_get_kernel_cents(g_u,opts.sz,opts.lmax);
        if any([d1(:);d2(:);d3(:)]<=0)
            error('fishy..');
        end
        d = struct;
        d.g = d1;
        d.g_l = d2;
        d.g_u = d3;
    otherwise
        d = hb_get_kernel_cents(g,sz,lmax);
        if any(d<=0)
            error('fishy..');
        end
end
varargout{1} = d;

