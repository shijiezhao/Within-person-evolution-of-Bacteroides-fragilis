cd('S03_metagenomes');
%% Step 1: load data and scripts
load('coveragematrix.mat');
load('candidate_mutation_table.mat');
load('TimeSeries.mat');
load('Isolates.mat');
load('optc.mat');   % Include the information of ancestor allele
%% Step 3: Set up colors pallete and parameters
% Test color pallete
Blue1=[225 237 235]/256;Blue2=[209 229 227]/256;Blue3=[194 220 217]/256;Red1=[246 220 212]/256;Red2=[213 161 157]/256;Red3=[242 205 192]/256;Grey1=[183 175 173]/256;
Grey2=[153 145 143]/256;RedD=[172 66 55]/256;Blue4=[100 139 173]/256;Blue5=[139 175 198]/256;white1=[210 210 210]/256;
c1=Red2; c2=Red1; c3=Blue1; c4=Grey1*1.1; c5=Grey2*1.1; c6=Blue3; c7=Blue5;
c10=[181 175 198]/255;
a={c10,c4,c3,c4,c7,c6,c7,'black',white1};
cx=[0.3 0.3 0.3]*1.1;

%% Plot Extended Data Figure 9e
figure(1); hold on;
dg=0;sp=15; % First smoothing parameter
SNP_trajectories = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
SM_Traj= {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
bin=1;smt=0;N=0;mk=0;MaxD=0;

box on;
xlim([-5.3 150]);ylim([-0.0 1.0]);set(gca,'xtick',[0 25 50 75 100 125 144 150])
set(gca,'Xticklabel',{});set(gca,'ytick',[0 0.5 1]);set(gca,'LineWidth',1);set( gca, 'TickDir', 'out' );
% Formating
for i = -4.5:0.5:-0.6;
    h=bar(i,1,0.5);
    h.FaceColor=[1 1 1];
    h.EdgeColor=[1 1 1];
    h.LineWidth=1;
end
for i = 144.5:0.5:150;
    h=bar(i,1,0.5);
    h.FaceColor=[1 1 1];
    h.EdgeColor=[1 1 1];
    h.LineWidth=1;
end

KKK=[]; 
for i = 1:21;
    [MIC,Days] = GetSNPtrajectory(i,counts,bin,N,Date,optc,smt);
    label = optc(i,4) % Where the mutation is
    if label <8;
%         mk=mk+1;
        DS=Days;CS=MIC;
%         plot(Days(1:length(Days)),MIC(1:length(Days)),'-','LineWidth',5,'Color',a{optc(i,4)},'LineWidth',2);hold on;
        CS = smooth(DS,CS,sp,'sgolay',dg);
        CS(CS>1)=1;CS(CS<0)=0;
        if label==7;CS(DS<50)=0;end
        
        if label==1 | label==2 | label==5;
        plot(DS,CS,'Color',a{optc(i,4)},'LineWidth',2.5);
        end
        KKK=[KKK;mean(CS)];
    end
    label = optc(i,4) % Where the mutation is
    if label < 8;
        SNP_trajectories{label} = [SNP_trajectories{label} [DS;CS']];
        MaxD=max(MaxD,max(DS));
    end
    
end

% Plot isolates inference
plot([0 0],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
plot([144 144],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
%Plot the isolates sequencing results
jig=60;ms=10.5;mk='o';mke=[0.8 0.8 0.8]/1.5;

for i = 1:2;
    for j = 1:5;
        if j==1 | j==2 | j==5;
            h=plot(Isolate_date(i),Isolate_data(i,j),mk,'MarkerSize',ms,'Color',a{j});set(h,'MarkerEdgeColor',mke,'MarkerFaceColor',a{j});
        end
    end
end

set(gca,'Fontsize',15)
% box on;

%% Continuous
FOOS={}
for i = 1:5;
    XX=[SNP_trajectories{i}(1,:)'];YY=[SNP_trajectories{i}(2,:)'];
    foo = fit(XX,YY,'fourier7');FOOS{i}=foo;
end


%% Step 9: Fit and plot Muller
Y=[];
Days=0:0.1:144;
for day = 0:0.1:144;
    f1 = max(FOOS{1}(day),0);f2 = max(FOOS{2}(day),0);f3 = max(FOOS{3}(day),0);f4 = max(FOOS{4}(day),0);
    f5 = max(FOOS{5}(day),0)
    f1=f1-f5;
    total=f1+f2+f5;
    f1=f1/total;f2=f2/total;f5=f5/total;
    f3=0;f4=0;
    Y=[Y;[f2 f1/2 f5 f1/2]];
end

figure(2); hold on;
bb=bar(Y,'stacked','Barwidth',1,'EdgeAlpha',0);
c = [a{2};a{1};a{5};a{1}];
xlim([0 1440]);ylim([0 1])
bb(1).Parent.Parent.Colormap = c;
box off;set(gca,'xtick',[]);set(gca,'ytick',[])
