function [allthresholdsmet] = div_test_thresholds(cnts, p, coveragethresholds, controlsample)


%Tami Lieberman April 2012, could be written more matlaby, without going
%one sample at a time --- this hooks up nicely into a single sample version
%used during initial data curation

if nargin < 4
    controlsample=1;
end


NSamples=size(cnts,3);
NPositions=size(cnts,2);

allthresholdsmet=zeros(NPositions,NSamples);

%do control first
[controlFreq, ~, ~] = div_major_allele_freq(cnts(:,:,controlsample));


%assess others, using controlMAF
for sample=1:NSamples;
    v=div_single_sample_test_thresholds(cnts(:,:,sample), p, controlFreq, coveragethresholds(:,sample));
    allthresholdsmet(:,sample)=v;
end

 
