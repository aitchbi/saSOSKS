% Also computes and returns the std of the energy. This is good to see how
% much the ensemble energy is consistent across signals, and if this can be
% used as a feature to distinguish differnet signal sets. 
% [18 Feb 2017]
%
%


%SPG_CONSTRUCT_WARPING Construct energy-equalizing transformation for 
% constructing signal-adapted system of spectral kernels as proposed in:
%
% Behjat, et al., 'Signal-adapted tight frames on graphs', 
% IEEE Trans. Signal Process., 2016,  doi: 10.1109/TSP.2016.2591513.
%
% To see example usage and the required inputs, please refer to 
% either of the following demo functions:
%
% spg_demo_ (minnesota | alameda |cerebellum)
% spg_demo_ (your_data_frame | your_data_decompose)
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

function [signalSets, G] = spgg_construct_warping(signalSets,G,varargin)
control_params = {'sFactor',11,'saveSpectrums','no','stdNorm',0};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

uniqueEigs =unique(G.E);
eigCounts = histc(G.E,uniqueEigs);
repetitiveEigs = uniqueEigs(eigCounts>1);

% Repetitive eigenvalues and their counts
for iRep =1:numel(repetitiveEigs)
    repCount(iRep) = numel(find(G.E==repetitiveEigs(iRep)));
end

if numel(repetitiveEigs)>0
    G.repetitiveEigs = repetitiveEigs;
    G.repetitiveEigsCounts = vec(repCount);
else
    G.repetitiveEigs = 'none';
end

for kSignalSet=1:numel(signalSets)
     
    nSignals = size(signalSets{kSignalSet}.signals,2);
    
    signalSets{kSignalSet}.En  = zeros(G.N,1); 
    signalSets{kSignalSet}.C_energy = [];      
    signalSets{kSignalSet}.warpE = [];         
    signalSets{kSignalSet}.warpE_smooth = []; 
    signalSets{kSignalSet}.energies  = zeros(G.N,1);
    
    for iSignal = 1:nSignals
        % STORE signal realization and sigma of added noise
        dummy2 = signalSets{kSignalSet}.signals(:,iSignal);
        
        % step 4. normalize
        % -- Note: this normalization ensures that all signals 
        % contribute equally to the construction of the warping function.
        dummy2 = dummy2/norm(dummy2,2);
        
        % step 5. compute spectrum
        if ~numel(find(isnan(dummy2))) % skip the signal in case its corupted
            dummy3 = G.U'*dummy2; 
            dummy4 = vec(dummy3).^2;
            signalSets{kSignalSet}.energies(:,iSignal) = dummy4;
            signalSets{kSignalSet}.En = signalSets{kSignalSet}.En + dummy4;
        end
        
        % save normalized spectrums
        if strcmp(saveSpectrums,'yes')
            signalSets{kSignalSet}.spectrums(:,iSignal) = vec(dummy3);
        end
    end
    
    % Ensemble energy
    signalSets{kSignalSet}.En  = signalSets{kSignalSet}.En/(nSignals);
    
    % Cumulative energy
    signalSets{kSignalSet}.C_energy = cumsum(signalSets{kSignalSet}.En);
    
    % std of energies
    signalSets{kSignalSet}.stdE = std(signalSets{kSignalSet}.energies,stdNorm,2);
    
    % warping
    signalSets{kSignalSet}.warpE = (signalSets{kSignalSet}.C_energy - signalSets{kSignalSet}.C_energy(1)) / ...
        (signalSets{kSignalSet}.C_energy(end)-signalSets{kSignalSet}.C_energy(1));
    
    % smoothed warping
    signalSets{kSignalSet}.warpE_smooth = smooth(signalSets{kSignalSet}.warpE,sFactor);
    
    % adjust for repetitive eigenvalues - such that the same warp value reulsts for repetitive equal eigenvalues 
    for iRep =1:numel(repetitiveEigs)
        indDummy = find(G.E==repetitiveEigs(iRep));
        signalSets{kSignalSet}.warpE(indDummy) = mean(signalSets{kSignalSet}.warpE(indDummy));
        signalSets{kSignalSet}.warpE_smooth(indDummy) = mean(signalSets{kSignalSet}.warpE_smooth(indDummy));
    end
    
end

end