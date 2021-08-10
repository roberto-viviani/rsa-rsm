function str = rsa_defaults()
% RSA_DEFAULTS Parameters used in the rsa package.
% Change the code below to change these parameters.
%
% 2020-21 Roberto Viviani - December 2020
% Institute of Psychology, University of Innsbruck
% Last modified: March 23th 2021

    %the minimum number of voxels included in the searchlight. Searchlights
    %that are smaller than this are excluded from the estimated volume.
    %Ignored when the searchlight size is less than 2 voxels (box) or 1
    %voxels (sphere); in this case, the min voxels threshold is set to 7.
    str.minsearchlightsize = 27;
    
    %the off-diagonal minimum value for a representational map RSM to be
    %considered as carrying information. If max(abs(triu(RSM))) is less
    %than this value, the RSM is rejected an an error is thrown.
    str.offdiagtol = 0.001;
    
    %the tolerance value for the largest scale ratio among the
    %representational maps. A large scale ratio may occur if sums of
    %squares and cross products from data are used as confound maps, and
    %the magnitude of these data leads to excessively large values. Set
    %this value to a negative value to switch off this diagnostic. Then,
    %maps will be scaled to ensure stable estimation of correlations.
    str.predscaletol = 10000;
    
    %the tolerance value to declare the set of predictors from the
    %representational maps to be rank-deficient. This value was set here to
    %a much larger value than the MATLAB proposed default, reflecting
    %statistical rather than numerical criteria of stability. To
    %repristinate the default MATLAB tolerance threshold, set this value to
    %zero.
    str.ranktol = 1 / 10000;

