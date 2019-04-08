%% 1. Variable parameters Parameters
% Much of this will need to vary according to your reference genome,
% coverage, particular samples, and how well the experiment was done.
% clear all

% 1.1. Parameters for finding de novo mutations between samples
min_average_coverage_to_include_sample = 10;      % Filter out samples have lower than that number of coverage
max_fraction_ambigious_samples = .33;            % If more than x% of the samples have ambiguous NT, discard the candidate location
min_median_coverage_position = 10;               % Remove candidate locations have lower than this coverage
min_qual_for_call = 60;                          % Remove sample*candidate that has lower than this quality
min_maf_for_call = .9;                           % Remove sample*candidate
if strcmp(Donor,'S10');
    min_maf_for_call = .95;                           % Raise the threshold to exclude false positive, only for S10
end
min_cov_for_call_per_strand = 7;                 % Remove sample*candidate
% max_ambigious_locations = 0.6;                   % Remove samples if more than this fraction of locations are ambigious

FQ_cutoff=60;      %min -FQ that samples supporting both alleles must have
% 1.2. Promoter parameter: how far upstream of the nearest gene to annotate
promotersize=300;

