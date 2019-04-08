% Meant to run from inside the case folder
% Pulls from sample_names.csv
% Pulls genome length form mutation_table_



run_postfix='14_06';
load(['mutation_table_' run_postfix])

SampleInfo = read_sample_names ;

% get all coverage per bp and the mode coverage for each sample
[all_coverage_per_bp, coverage_modes] = get_all_coverage(SampleInfo, GenomeLength);
save('coveragematrix.mat', 'all_coverage_per_bp', 'coverage_modes', '-v7.3'); 
