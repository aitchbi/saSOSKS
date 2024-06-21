% An extended version of of spg_filter_design.m from: 
% https://github.com/aitchbi/spg, 
% which includes a number of additional features. 

function [g,varargout] = spgg_filter_design(lmax,N_subbands,varargin)
control_params = {...
    'designtype',[],...
    'warping',[],...
    'G',[],...
    'E',[],...
    'subSampleWarping',0,...
    'fixForSparseSpectra',0,... 
    'minSparseFactor',[],... 
    'extrapolStep',[],... 
    'binbasedWarping',0,... 
    'warpingBins',[],... 
    'vers',1,... 
    'pou','over1stPower',... 
    'sOrder',3,... 
    'sz',0.001,... 
    'opts_uniMultRes', [],... 
    'opts_lphp', []}; 
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
if ~isempty(N_subbands) % empty for designtype 'uniform_multires'
    N_subbands = N_subbands-1;
end
g = cell(N_subbands+1,1);

if subSampleWarping==1 && fixForSparseSpectra==1
    error('This scenario has not been implemented.')
    % HB: this may even not make sense to implement.
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
        % Uniform Meyer-type system of spectral kernels wherein sum of
        % squared abs values forms a partition of unity
        
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