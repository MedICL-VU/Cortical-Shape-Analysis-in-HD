% Random field theory for FWER p-value correction for chi statistical fields
%
% [pfwer, peak, clus] = rft_fwer(chi, df, alpha, fwhm, mask, slm)
%
% Input:
% chi   : array corresponding to the chi field
% df    : degree of freedom
% alpha : p-value threshold
% fwhm  : FWHM of the smoothing parameter
% mask  : mask of the ROI
% slm   : SurfStat linear model
%
% Output:
% pfwer.C = cluster extends
% pfwer.P = local maxima
%
%
% Author: Kilian Hett (kilian.hett@vanderbilt.edu)


function [pfwer, peak, clus] = rft_fwer(chi, df, alpha, fwhm, mask, slm)

[~,reselspvert,~] = SurfStatResels(slm,mask);

% Cluster extend (Frisoni et al. 1996)
pfwer.C = ones(size(chi))';
W = fwhm/sqrt(4*log(2));
p = chi2cdf(chi, df, 'upper');
[~,~,~,q] = fdr_bh(p,alpha);

bw = findcompcon((q<alpha) .* (mask'>0), slm.tri);

m = bw.NumObjects;
S = sum(mask>0);
u = chi2inv(1-alpha, df);
EN = S * exp(-u^2/2) / (u*sqrt(2*pi));
Em = S * (2*pi)^((-3)/2) * W^(-2) * u * exp(-(u^2)/2);
B = gamma(2)*(Em/EN);
for i=1:m
    s = sum(reselspvert(bw.VoxelIdxList{i}));
    pfwer.C(bw.VoxelIdxList{i}) = exp(-B*s); 
end


% Peak (Worsley et al. 1999)
[peaks, locs] = findpeaks(chi,'MinPeakDistance',fwhm^2);
pfwer.P = ones(size(chi))';

idx = find(peaks>u);
for i=1:length(idx)
    pfwer.P(locs(idx(i))) = sum(mask(:)>0)*reselspvert(locs(idx(i))) .*...
                            expected_ec_noniso(peaks(idx(i)), df);
end
pfwer.P = pfwer.P .* (pfwer.P<alpha) + (pfwer.P>=alpha);

clus = [];
peak = [];

end


% Expected euler characteristic
function ec = expected_ec_noniso(t, df)
    p1 = (4*log(2)) ./ (2*pi);
    p2 = (t.^((df-2)/2) .* exp(-t./(2))) ./ (2.^((df-2)/2) .* gamma(df./2));
    p3 = t - (df-1);
    ec =  p1 .* p2 .* p3;
end

