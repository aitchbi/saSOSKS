function x = spgg_apply_warping(x,warping,lmax,E,warpingBins,...
    subSampleWarping,fixForSparseSpectra,minSparseFactor,extrapolStep)

if isempty(minSparseFactor) 
    minSparseFactor = lmax/10;
    if isempty(extrapolStep)
        extrapolStep = minSparseFactor/2;
    end
end

if ~isempty(E)
    E = E(:);
end

if ~ischar(warping)
    warping = warping(:);
end

if strcmp(warping,'none')
    
elseif ~isempty(warpingBins)
    interp_x = warpingBins;
    interp_y = warping*lmax;
    x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,x));
else
    theEigs = E;
    diffOfEigs = theEigs(2:end)-theEigs(1:end-1);
    indNonEqual = [1; find(logical(diffOfEigs(:))) + 1]; % June 2017: added (:)
    
    if subSampleWarping
        dummy1 = E(indNonEqual);
        interp_x = dummy1(1:subSampleWarping:end-1);
        interp_x(end+1) = dummy1(end);
        
        dummy2 = warping(indNonEqual)*lmax;
        interp_y = dummy2(1:subSampleWarping:end-1);
        interp_y(end+1) = dummy2(end);
        
        x = x(x <= interp_x(end));
        x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,x));
        
    elseif fixForSparseSpectra
        
        interp_x = E(indNonEqual(1)); % should be 0
        interp_y = lmax*warping(1);   % should be 0
        
        for i = 2:numel(indNonEqual)
            
            check = E(indNonEqual(i)) - E(indNonEqual(i-1));
            
            if  check < minSparseFactor
                interp_x = [interp_x, E(indNonEqual(i))]; %#ok<*AGROW> 
                interp_y = [interp_y, lmax*warping(indNonEqual(i))];
            else
                dummy = E(indNonEqual(i-1))+extrapolStep:extrapolStep:E(indNonEqual(i));
                interp_x = [interp_x, dummy];
                interp_y = [interp_y, lmax*repmat(warping(indNonEqual(i-1)),1,numel(dummy))];
            end
        end
        x = abs(gsp_mono_cubic_warp_fn(vec(interp_x),vec(interp_y),x));
    else
        interp_x = E(indNonEqual);
        interp_y = warping(indNonEqual)*lmax;
        x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,x));
    end
end