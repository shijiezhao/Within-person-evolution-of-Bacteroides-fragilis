%% Timepoint for figure S5G-S5J
tmQ=[0;6;9;12;15;22];
load('competition_exp/FigureS5GJ_Sequencing.mat')
SCRIPTSDIRECTORY = ['./scripts_generating_figures/FigureS5_prep'];
path(SCRIPTSDIRECTORY,path);  
%% About Quantify:
% First and Third columns: reads aligned to the ancestral allele
% Second and Fourth columns: reads aligned to the mutant allele
% Fifth column: time point
% Sixth column: combinations
% Seventh column: different types of amplicon
% Eighth column: technical replicates
% Combinations:
% #1-#3: Figure S5G
% #4-#15: Figure S5H, S5I
% #16-#21: Figure S5J

Comb1_3_S1a_vs_S1b()
Comb4_6_S1a_vs_S2A()
Comb7_9_S1b_vs_S2A()
Comb10_12_S1a_vs_S2B()
Comb13_15_S1b_vs_S2B()
Comb16()
Comb17()
Comb18()
Comb19()
Comb20()
Comb21()
Trio_ratio_adjustment()
close all force;

Ratios{1} = Ratio_1_mean_ts;
Ratios{2} = Ratio_2_mean_ts;
Ratios{3} = Ratio_3_mean_ts;
Ratios{4} = Ratio_4_mean_ts;
Ratios{5} = Ratio_5_mean_ts;
Ratios{6} = Ratio_6_mean_ts;
Ratios{7} = Ratio_7_mean_ts;
Ratios{8} = Ratio_8_mean_ts;
Ratios{9} = Ratio_9_mean_ts;
Ratios{10} = Ratio_10_mean_ts;
Ratios{11} = Ratio_11_mean_ts;
Ratios{12} = Ratio_12_mean_ts;
Ratios{13} = Ratio_13_mean_ts;
Ratios{14} = Ratio_14_mean_ts;
Ratios{15} = Ratio_15_mean_ts;
%%
FigureS5K()
%%
FigureS5GHIJ()
%%