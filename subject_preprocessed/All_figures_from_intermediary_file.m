SCRIPTSDIRECTORY = ['./scripts_generating_figures'];
path(SCRIPTSDIRECTORY,path);       % Add the srcipt folder to the search path
majordir = char(pwd);
stop
%% Step 0: prepare intermediary files
RefGenomes={'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'};
Donors= {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'};
ClosestGenomes={'S07','S03','S02','S11','S08','S08','S08','S05','S10','S09','S04','S04'};
% 1.1. Reference genome & Donor name
for i  = 1:12;
RefGen = RefGenomes{i};Donor = Donors{i};closest = ClosestGenomes{i};
PrePare();cd(majordir);close all force;
end

%: Phylogenetic ree files for Figure 2A, 3A, 4A, 5A and S1-S3 can be found 
%from the individual donor folder

%% Step 1: Plot Figure 1a and Figure 1b
Figure1();cd(majordir);
phytreewrite('Figures/Figure_1a.tree', PhyloTree,'BranchNames',0)
figure(1);saveas(gcf,'Figures/Figure_1b.png'); 
xlim([16000 19000]);ylim([0 50]);saveas(gcf,'Figures/Figure_1b_zoom.png'); close all force;

%% Step 2: Plot Figure 2b + Figure 2e
Figure2bd();cd(majordir);
figure(1);saveas(gcf,'Figures/Figure_2b.png'); 
figure(2);saveas(gcf,'Figures/Figure_2e.png'); close all force;
%% Step 3: Plot Figure 2c
Figure2c();cd(majordir);
saveas(gcf,'Figures/Figure_2c.png'); close all force;
%% Step 4: Plot Figure 2f
Figure2f();cd(majordir);
saveas(gcf,'Figures/Figure_2f.png'); close all force;
%% Step 5: Plot additional figure about molecular clock and dMRCA
Additional_figure1();cd(majordir);
saveas(gcf,'Figures/Additional_figure_1.png'); close all force;
%% Step 6: Plot Figure S7
FigureS7();cd(majordir);
saveas(gcf,'Figures/Figure_S7.png'); close all force;
%% Step 7:Prepare simulation for plotting figure 4 and figure S4
Simulation_fig4_figs4();cd(majordir);
%% Step 8: Plot Figure 4c
Figure4c();cd(majordir);
saveas(gcf,'Figures/Figure_4c.png'); close all force;
%% Step 9: Plot Figure 3b
Figure4d();cd(majordir);
saveas(gcf,'Figures/Figure_4d.png'); close all force;
%% Step 10: Plot Figure S4
FigureS4();cd(majordir);
saveas(gcf,'Figures/Figure_S4.png'); close all force;
%% Step 13: Plot Figure 5a,5b,5c
Figure5();cd(majordir);
figure(1);saveas(gcf,'Figures/Figure_5b.png'); 
figure(2);saveas(gcf,'Figures/Figure_5c_Muller.png'); 
figure(3);saveas(gcf,'Figures/Figure_5c_PreSampling.png'); close all force;
%% Step 14: Plot Figure S5A, S5B, S5C
FigureS5abc();cd(majordir);
figure(1);saveas(gcf,'Figures/Figure_S5abc.png'); close all force;
%% Step 15: Plot Figure S5E, S5F
FigureS5ef();cd(majordir);
figure(1);saveas(gcf,'Figures/Figure_S5E.png'); 
figure(2);saveas(gcf,'Figures/Figure_S5F.png'); close all force;
%% Step 16: Plot Figure 6C
Figure6c();cd(majordir)
figure(1);saveas(gcf,'Figures/Figure_6C.png'); close all force;
%% Step 17: Plot Figure 5E-5H 
Figure5EFGH();cd(majordir)
figure(3);saveas(gcf,'Figures/Figure_5FGH.png'); close all force;

%% Step 18: Plot Figure 3B
Figure3b();cd(majordir)
figure(1);saveas(gcf,'Figures/Figure_3B.png'); close all force;

%% Step 19: Plot Figure 3C
Figure3c();cd(majordir)
figure(1);saveas(gcf,'Figures/Figure_3C.png'); close all force;

%% Step 20: Plot Figure S5G-S5K
FigureS5S();cd(majordir)
figure(22);saveas(gcf,'Figures/Figure_S5G.png');
figure(233445);saveas(gcf,'Figures/Figure_S5H.png'); 
figure(23);saveas(gcf,'Figures/Figure_S5I.png'); 
figure(234);saveas(gcf,'Figures/Figure_S5J.png'); 
figure(2);saveas(gcf,'Figures/Figure_S5K.png'); close all force;
%% Step 21: Plot
FigureS5L();cd(majordir)
figure(32);saveas(gcf,'Figures/Figure_S5L.png');close all force;