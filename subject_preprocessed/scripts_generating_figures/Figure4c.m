%% 0. Set up colors
Red1=[244 113 126]/256; Red2=[212 160 156]/256;
Blue1=[64 134 153]/256; Blue2=[193 219 216]/256; Blue3=[100 142 180]/256;
Grey1=[64 64 64]/256; Grey2=[200 200 200]/256;

c1 = Grey1; c2=Grey2;c3=Blue2;c4=Red2;
mycolors=([c1;c2;c3;c4]);


%% 1. Major paper figure: Check total within parallels
sim_totwithin = []
figure(45); hold on;
numtrials = 1000;
for sim = 1:numtrials;
    stw = 0;
    for i = 1:12;
        RAN=randi(numtrials)
        stw = stw + TotSim{i}{1}{RAN,3}; % get num simulated number of genes under parallel evolution
    end
    sim_totwithin = [sim_totwithin; stw];
end

for stw = 0:max(sim_totwithin);
    h = bar(stw,sum(sim_totwithin==stw)/length(sim_totwithin),0.5); % Plot the simulations
    set(h,'FaceColor',Grey2);
    set(h,'LineWidth',0.05);
    set(h,'EdgeAlpha',0);
end

% Plot the observed values
within_ltg=[];
for i = 1:12;
    tmp=TotReal{i}{1}{1};
    [tmpltg, srp]=unique(tmp(:,1));
    for j=1:length(tmpltg);
        if tmpltg(j)>0 & sum(tmpltg(j)==tmp(:,1))>1 & sum(tmpltg(j)==tmp(:,1))/tmp(srp(j),2)>1/genelengththreshold;
            within_ltg=[within_ltg;tmpltg(j)];
        end
    end
end
within_ltg=unique(within_ltg);
    


%% Format 

h = bar(length(within_ltg),1,0.25)
set(h,'FaceColor',Red2);
set(h,'EdgeAlpha',0);


xlim([-1.25 20])
ylim([0 0.3])
set(gca,'Ytick',[0 0.1 0.2 0.3])
set(gca,'Xtick',[-1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20])
set(gca,'Xticklabel',{'' '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'})
set(gca,'TickDir','out')
set(gca,'FontSize',10)