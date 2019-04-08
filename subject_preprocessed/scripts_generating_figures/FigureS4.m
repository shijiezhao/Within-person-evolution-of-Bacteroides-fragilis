%% Plot the panels for supplementary figure 8
%% Subplot 2: Plot M>=3 mut across all subjects
sim_tot_across = [];    % Simulation results
tot_across=0;           % observed
figure(1);   % Initiate the figure
subplot(2,3,2); hold on; set(gca,'TickDir','out')
% Parameters
M=3;
M3_ltg=[];
% Count for simulation
for sim = 1:numtrials;
    sltgs = [];
    sta = 0;
    for i = 1:12;
        RAND=randi(numtrials)
        sltgs = [sltgs; TotSim{i}{1}{RAND,1}];
    end
    [sultgs, srp] = unique(sltgs(:,1)); % srp: the position of the first appearance
    for i = 1:length(sultgs);
        if sultgs(i)>0 & sum(sultgs(i)==sltgs(:,1))>(M-1) & sum(sultgs(i)==sltgs(:,1))/sltgs(srp(i),2) > 1/genelengththreshold;
            sta = sta + 1;
        end
    end
    sim_tot_across = [sim_tot_across;sta];
end

tltgs=[]
for i = 1:12;
    tltgs = [tltgs; TotReal{i}{1}{1}];
end
[ultgs, rp] = unique(tltgs(:,1));
for i = 1:length(ultgs);
    if ultgs(i)>0 & sum(ultgs(i)==tltgs(:,1))>(M-1) & sum(ultgs(i)==tltgs(:,1))/tltgs(rp(i),2) > 1/genelengththreshold;
        tot_across = tot_across + 1;
        M3_ltg=[M3_ltg;ultgs(i)];
    end
end

% Plot the simulation bars
for sta = 0:max(sim_tot_across);
    h = bar(sta,sum(sim_tot_across==sta)/length(sim_tot_across),0.5);set(h,'FaceColor',Grey2);
    set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end
% Plot the observed bar
h=bar(tot_across,1,0.25,'r');set(h,'FaceColor',Blue1);set(h,'EdgeAlpha',0);
% Format
xlim([-1 15]);ylim([0 0.6])
set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1]);set(gca,'Xtick',[ 0 5 10 15])
set(gca,'Xticklabel',{ '0'  '5'  '10'  '15'});set(gca,'FontSize',10)


%% Subplot 2: Plot M>=2 mut across all subjects
sim_tot_across = [];tot_across=0;
subplot(2,3,1);hold on;set(gca,'TickDir','out')
% Initiate parameters
M=2;
M2_ltg=[];
% Create the simulation records
for sim = 1:numtrials;
    sltgs = [];
    sta = 0;
    for i = 1:12;
        RAND=randi(numtrials)
        sltgs = [sltgs; TotSim{i}{1}{RAND,1}];
    end
    [sultgs, srp] = unique(sltgs(:,1)); % srp: the position of the first appearance
    for i = 1:length(sultgs);
        if sultgs(i)>0 & sum(sultgs(i)==sltgs(:,1))>(M-1) & sum(sultgs(i)==sltgs(:,1))/sltgs(srp(i),2) > 1/genelengththreshold;
            sta = sta + 1;
        end
    end
    sim_tot_across = [sim_tot_across;sta];
end
tltgs=[]
for i = 1:12;
    tltgs = [tltgs; TotReal{i}{1}{1}];
end
[ultgs, rp] = unique(tltgs(:,1));
for i = 1:length(ultgs);
    if ultgs(i)>0 & sum(ultgs(i)==tltgs(:,1))>(M-1) & sum(ultgs(i)==tltgs(:,1))/tltgs(rp(i),2) > 1/genelengththreshold;
        tot_across = tot_across + 1;
        M2_ltg=[M2_ltg; ultgs(i)];
    end
end

% Plot the simulation
for sta = 0:max(sim_tot_across);
    h = bar(sta,sum(sim_tot_across==sta)/length(sim_tot_across),0.5);
    set(h,'FaceColor',Grey2);
    set(h,'LineWidth',0.05);
    set(h,'EdgeAlpha',0);
end
% Plot observed
h=bar(tot_across,1,0.25,'r');set(h,'FaceColor',Blue1);set(h,'EdgeAlpha',0);
% Format
xlim([-1.25 40]);ylim([0 0.2]);set(gca,'Ytick',[0 0.05 0.1 0.15 0.2])
set(gca,'Xtick',[ 0 5 10 15 20 25 30 35 40]);set(gca,'Xticklabel',{ '0' '5' '10' '15' '20' '25' '30' '35' '40'})
set(gca,'TickDir','out');set(gca,'FontSize',10)

%% Subplot 3: Plot m>=3
sim_tot_across = [];tot_across=0;
subplot(2,3,4);hold on;set(gca,'TickDir','out')
% Initialize paramter
m=3;
m3_ltg=[];
% Get simulation results
for sim = 1:numtrials;
    sltgs = [];
    distinct_ltgs=[]; genelength=[];
    sta = 0;
    for i = 1:12;
        RAND=randi(numtrials);
        tmp=TotSim{i}{1}{RAND,1};
        [sultgs, srp] = unique(tmp(:,1));       % Let each locustag appears only once from each donor
        distinct_ltgs = [distinct_ltgs;sultgs];
        genelength = [genelength;tmp(srp,2)];
    end
    [tultgs,srp]=unique(distinct_ltgs); % find the unique locustags
    for i = 1:length(tultgs);
        if tultgs(i)>0 & sum(tultgs(i)==distinct_ltgs)>(m-1) & sum(tultgs(i)==distinct_ltgs)/genelength(srp(i)) > 1/genelengththreshold;
            sta = sta + 1;
        end
    end
    sim_tot_across = [sim_tot_across;sta];
end
distinct_ltgs=[];genelength=[];
for i = 1:12;
    tmp = TotReal{i}{1}{1};
    [sultgs, srp] = unique(tmp(:,1));       % Let each locustag appears only once from each donor
    distinct_ltgs = [distinct_ltgs;sultgs];
    genelength = [genelength;tmp(srp,2)];
end
[tultgs,srp]=unique(distinct_ltgs); % find the unique locustags
for i = 1:length(tultgs);
    if tultgs(i)>0 & sum(tultgs(i)==distinct_ltgs)>(m-1) & sum(tultgs(i)==distinct_ltgs)/genelength(srp(i)) > 1/genelengththreshold;
        tot_across = tot_across + 1;
        m3_ltg=[m3_ltg;tultgs(i)];
    end
end

% Plot simulation;
for sta = 0:max(sim_tot_across);
    h = bar(sta,sum(sim_tot_across==sta)/length(sim_tot_across),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end
% Plot observed
h=bar(tot_across,1,0.25,'r');set(h,'FaceColor',Blue1);;set(h,'EdgeAlpha',0);
xlim([-1 5]);ylim([0 0.8]);set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1]);set(gca,'Xtick',[ 0 1 2 3 4 5])
set(gca,'TickDir','out');set(gca,'FontSize',10);

%% Subplot 4: Plot m>=2
sim_tot_across = [];tot_across=0;
subplot(2,3,3);hold on;set(gca,'TickDir','out')
% Initialize parameter
m=2;
m2_ltg=[];
for sim = 1:numtrials;
    sltgs = [];
    distinct_ltgs=[]; genelength=[];
    sta = 0;
    for i = 1:12;
        RAND=randi(numtrials)
        tmp=TotSim{i}{1}{RAND,1};
        [sultgs, srp] = unique(tmp(:,1));       % Let each locustag appears only once from each donor
        distinct_ltgs = [distinct_ltgs;sultgs];
        genelength = [genelength;tmp(srp,2)];
    end
    [tultgs,srp]=unique(distinct_ltgs); % find the unique locustags
    for i = 1:length(tultgs);
        if tultgs(i)>0 & sum(tultgs(i)==distinct_ltgs)>(m-1) & sum(tultgs(i)==distinct_ltgs)/genelength(srp(i)) > 1/genelengththreshold;
            sta = sta + 1;
        end
    end
    sim_tot_across = [sim_tot_across;sta];
end
distinct_ltgs=[];genelength=[];
for i = 1:12;
    tmp = TotReal{i}{1}{1};
    [sultgs, srp] = unique(tmp(:,1));       % Let each locustag appears only once from each donor
    distinct_ltgs = [distinct_ltgs;sultgs];
    genelength = [genelength;tmp(srp,2)];
end
[tultgs,srp]=unique(distinct_ltgs); % find the unique locustags
for i = 1:length(tultgs);
    if tultgs(i)>0 & sum(tultgs(i)==distinct_ltgs)>(m-1) & sum(tultgs(i)==distinct_ltgs)/genelength(srp(i)) > 1/genelengththreshold;
        tot_across = tot_across + 1;
        m2_ltg=[m2_ltg;tultgs(i)];
    end
end
% Plot simulations
for sta = 0:max(sim_tot_across);
    h = bar(sta,sum(sim_tot_across==sta)/length(sim_tot_across),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end
% Plot observed
h=bar(tot_across,1,0.25,'r');set(h,'FaceColor',Blue1);;set(h,'EdgeAlpha',0);
xlim([-1 40]);ylim([0 0.2]);set(gca,'Ytick',[0 0.05 0.1 0.15 0.2])
set(gca,'Xtick',[0 5 10 15 20 25 30 35 40]);set(gca,'Xticklabel',{'0' '5' '10' '15' '20' '25' '30' '35' '40'});set(gca,'TickDir','out');set(gca,'FontSize',10)

%% Supbplot 5: Inter-genic-region
subplot(2,3,5);hold on;
sim_int = [];
for sim = 1:numtrials;
    intg = 0;
    for i = 1:12;
        RAND=randi(numtrials)
        tmp=TotSim{i}{1}{RAND,1};
        intg=intg+sum(tmp(:,1)==0);
    end
    sim_int=[sim_int;intg];
end
real_int=0;
for i = 1:12;
    tmp = TotReal{i}{1}{1};
    intg=sum(tmp(:,1)==0);
    real_int=real_int+intg;
end
% Plot simulation
for sta = 0:max(sim_int);
    h = bar(sta,sum(sim_int==sta)/length(sim_int),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end
% Plot observed
h=bar(real_int,1,0.25,'r');set(h,'FaceColor',Blue1);set(h,'EdgeAlpha',0);
xlim([-1 120]);ylim([0 0.1]);set(gca,'Ytick',[0 0.05 0.1 0.15 0.2])
set(gca,'Xtick',[0 20 40 60 80 100 120]);set(gca,'TickDir','out');set(gca,'FontSize',10)

%%
% expectedNS(annotation_full_all) = 2.1544

Red1=[244 113 126]/256; Red2=[212 160 156]/256;
Blue1=[224 236 234]/256; Blue2=[193 219 216]/256; Blue3=[100 142 180]/256;
Grey1=[64 64 64]/256; Grey2=[200 200 200]/256;

%%
subplot(2,3,6);hold on; set(gca,'TickDir','out');ylim([log(0.1) log(128)])
set(gca,'Ytick',[log(0.125) log(0.25) log(.5) log(1) log(2) log(4) log(8) log(16) log(32) log(64) log(128)])
set(gca,'Yticklabel',{'0.125' '0.25', '.5' '1' '2' '4' '8' '16' '32' '64' '128'})
% Cross:
N=0;S=0;
for i = 1:length(annotation_full_all);
    if annotation_full_all(i).type=='N';
        N=N+1;
    else if annotation_full_all(i).type=='S';
        S=S+1;
        end
    end
end
[x,y,z]=binomialCIdNdS(N,S,expected_dnds_across_lineages)
% Plot the bar with CI
h=bar(1,log(z),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1); plot([1 1],[log(x) log(y)],'Color','black','LineWidth',1)
% M2
M2_ind=[];
for i = 1:length(annotation_full);
    if length(annotation_full(i).locustag)>5;
        tmp = strsplit(annotation_full(i).locustag,'_')
        tmpn = str2num(tmp{2});
        if sum(tmpn==M2_ltg)>0;
            M2_ind=[M2_ind;i];
        end
    end
end
[dNdS_M2,sim_M2,x,y] = Prepare_dNdS(annotation_full,M2_ind)
% Plot the bar with CI
h=bar(3,log(dNdS_M2),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);plot([3 3],[log(x) log(y)],'Color','black','LineWidth',1)

% M3
M3_ind=[];
for i = 1:length(annotation_full);
    if length(annotation_full(i).locustag)>5;
        tmp = strsplit(annotation_full(i).locustag,'_')
        tmpn = str2num(tmp{2});
        if sum(tmpn==M3_ltg)>0;
            M3_ind=[M3_ind;i];
        end
    end
end

[dNdS_M3,sim_M3,x,y] = Prepare_dNdS(annotation_full,M3_ind)
% Plot the bar with CI
h=bar(4,log(dNdS_M3),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);plot([4 4],[log(x) log(y)],'Color','black','LineWidth',1)

% m2
m2_ind=[];
for i = 1:length(annotation_full);
    if length(annotation_full(i).locustag)>5;
        tmp = strsplit(annotation_full(i).locustag,'_')
        tmpn = str2num(tmp{2});
        if sum(tmpn==m2_ltg)>0;
            m2_ind=[m2_ind;i];
        end
    end
end
[dNdS_m2,sim_m2,x,y] = Prepare_dNdS(annotation_full,m2_ind)
% Plot the bar with CI
h=bar(5,log(dNdS_m2),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);plot([5 5],[log(x) log(y)],'Color','black','LineWidth',1)

% m3
m3_ind=[];
for i = 1:length(annotation_full);
    if length(annotation_full(i).locustag)>5;
        tmp = strsplit(annotation_full(i).locustag,'_')
        tmpn = str2num(tmp{2});
        if sum(tmpn==m3_ltg)>0;
            m3_ind=[m3_ind;i];
        end
    end
end
[dNdS_m3,sim_m3,x,y] = Prepare_dNdS(annotation_full,m3_ind)
% Plot the bar with CI
h=bar(6,log(dNdS_m3),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);plot([6 6],[log(x) log(y)],'Color','black','LineWidth',1)

% within_ltg: WITHIN LINEAGE PARALLEL
wi_ind=[];
for i = 1:length(annotation_full);
    if length(annotation_full(i).locustag)>5;
        tmp = strsplit(annotation_full(i).locustag,'_')
        tmpn = str2num(tmp{2});
        if sum(tmpn==within_ltg)>0;
            wi_ind=[wi_ind;i];
        end
    end
end
[dNdS_wi,sim_wi,x,y] = Prepare_dNdS(annotation_full,wi_ind)
% Plot the bar with CI
h=bar(7,log(dNdS_wi),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);plot([7 7],[log(x) log(y)],'Color','black','LineWidth',1)


% Singly mutated genes
si_ind=[];
for i = 1:length(annotation_full);
    if sum(i==wi_ind)==0 & sum(i==m2_ind)==0 & sum(i==m3_ind)==0 & sum(i==M2_ind)==0 & sum(i==M3_ind)==0;
    	si_ind=[si_ind;i];
	end
end
[dNdS_si,sim_si,x,y] = Prepare_dNdS(annotation_full,si_ind)
% Plot bar with CI
h=bar(2,log(dNdS_si),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);plot([2 2],[log(x) log(y)],'Color','black','LineWidth',1)

% M>=3 but mmax <2;
addi=[];
for i = 1:length(annotation_full);
    if sum(i==wi_ind)==0 & sum(i==M3_ind)==1;
    	addi=[addi;i];
	end
end
[dNdS_addi,sim_addi,x,y] = Prepare_dNdS(annotation_full,addi)
% Plot bar with CI
h=bar(8,log(dNdS_addi),0.5);set(h,'FaceColor',Grey2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);plot([8 8],[log(x) log(y)],'Color','black','LineWidth',1)
