cd('S01_metagenomes')
%% Step 1: load data and scripts
load('coveragematrix.mat');
load('candidate_mutation_table.mat');
load('TimeSeries.mat');
load('optc.mat');   % Include the information of ancestor allele
load('Isolate_Info.mat');
load('change.mat');
load('IsolateMatrix.mat');
counts=counts(:,change,:);
Quals=Quals(change,:);
p=p(change);

%% Step 3: Set up colors pallete and parameters
% Test color pallete
Blue1=[225 237 235]/256;Blue2=[209 229 227]/256;Blue3=[194 220 217]/256;Red1=[246 220 212]/256;Red2=[213 161 157]/256;Red3=[242 205 192]/256;Grey1=[183 175 173]/256;
Grey2=[153 145 143]/256;RedD=[172 66 55]/256;Blue4=[100 139 173]/256;Blue5=[139 175 198]/256;white1=[210 210 210]/256;
c1=Red2; c2=Red1; c3=Blue1; c4=Grey1*1.1; c5=Grey2*1.1; c6=Blue3; c7=Blue5;
c1=[181 175 198]/255;
a={c1,c2,c3,c4,c5,c6,c7,'black',white1};cx=[0.3 0.3 0.3]*1.1;

%% plot raw trajectories
figure(1); hold on;
dg=0;sp=30; % First smoothing parameter
SNP_trajectories = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
SM_Traj= {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
bin=1;smt=0;N=0;mk=0;MaxD=0;    % Initiate parameters
box on;
xlim([-21.3 621.3]);ylim([-0.0 1.0]);set(gca,'xtick',[0 100 200 300 400 500 600])
set(gca,'Xticklabel',{});set(gca,'ytick',[0 0.5 1]);set(gca,'LineWidth',1);set( gca, 'TickDir', 'out' );
% Just for formating purpose...
for i = -19.5:1:-0.5;
    h=bar(i,1,1);
    h.FaceColor=[1 1 1];
    h.EdgeColor=[1 1 1];
    h.LineWidth=1;
end
for i = 600.5:1:619.5;
    h=bar(i,1,1);
    h.FaceColor=[1 1 1];
    h.EdgeColor=[1 1 1];
    h.LineWidth=1;
end
% Plot trajectories for each SNP
KKK=[];
for i = 1:82;
    [MIC,Days] = GetSNPtrajectory(i,counts,bin,N,Date,optc,smt);    % No smoothing introduced in this step
    label = optc(i,4) % Where the mutation is
    KKK=[KKK;mean(MIC)];
    if label <8;
        DS=Days;CS=MIC;
        CS = smooth(DS,CS,sp,'sgolay',dg);
        CS(CS>1)=1;CS(CS<0)=0;
        if label==7;CS(DS<50)=0;end
        plot(DS,CS,'Color',a{optc(i,4)},'LineWidth',2.5);
    end
    label = optc(i,4) % Where the mutation is
    if label < 8;
        SNP_trajectories{label} = [SNP_trajectories{label} [DS;CS']];
        MaxD=max(MaxD,max(DS));
    end
    
end

% Plot the frequencies inferred by counting isolates
plot([mean(Isolate_date(1:2)) mean(Isolate_date(1:2))],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
plot([mean(Isolate_date(3:4)) mean(Isolate_date(3:4))],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
plot([mean(Isolate_date(5:6)) mean(Isolate_date(5:6))],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
plot([mean(Isolate_date(7)) mean(Isolate_date(7))],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
plot([mean(Isolate_date(8:9)) mean(Isolate_date(8:9))],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
plot([600 600],[-0.2 2],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
%Plot the isolates sequencing results
jig=60;
ms=10.5;
mk='o';
mke=[0.8 0.8 0.8]/1.5;

Isolate_data=zeros(6,7);
for i = 1:6;
    tmp1 = IsoMatrix{i};
    for j = 1:7;
        selectedisolates = find(optc(:,4)==j);
        relatedIsoMatrix = tmp1(selectedisolates,:);
        tmp2=[];
        for k=1:length(selectedisolates);
            tmp2 = [tmp2;sum(relatedIsoMatrix(k,:)~=optc(selectedisolates(k),3))/size(tmp1,2)];
        end
        Isolate_data(i,j)=mean(tmp2);
    end
end
for i = 1:7;
    MM=mean(Isolate_data(1,i))+mean(Isolate_data(1,i))*randn/jig;MM(MM>1)=1;MM(MM<0)=-MM;
    h=plot(mean(Isolate_date(1:2)),MM,mk,'MarkerSize',ms,'Color',a{i});set(h,'MarkerEdgeColor',mke,'MarkerFaceColor',a{i});
    MM=mean(Isolate_data(2,i))+mean(Isolate_data(2,i))*randn/jig;MM(MM>1)=1;MM(MM<0)=-MM;
    h=plot(mean(Isolate_date(3:4)),MM,mk,'MarkerSize',ms,'Color',a{i});set(h,'MarkerEdgeColor',mke,'MarkerFaceColor',a{i});
    MM=mean(Isolate_data(3,i))+mean(Isolate_data(3,i))*randn/jig;MM(MM>1)=1;MM(MM<0)=-MM;
    h=plot(mean(Isolate_date(5:6)),MM,mk,'MarkerSize',ms,'Color',a{i});set(h,'MarkerEdgeColor',mke,'MarkerFaceColor',a{i});
    MM=mean(Isolate_data(4,i))+mean(Isolate_data(4,i))*randn/jig;MM(MM>1)=1;MM(MM<0)=-MM;
    h=plot(mean(Isolate_date(7)),MM,mk,'MarkerSize',ms,'Color',a{i});set(h,'MarkerEdgeColor',mke,'MarkerFaceColor',a{i});
    MM=mean(Isolate_data(5,i))+mean(Isolate_data(5,i))*randn/jig;MM(MM>1)=1;MM(MM<0)=-MM;
    h=plot(mean(Isolate_date(8:9)),MM,mk,'MarkerSize',ms,'Color',a{i});set(h,'MarkerEdgeColor',mke,'MarkerFaceColor',a{i});
    MM=mean(Isolate_data(6,i))+mean(Isolate_data(6,i))*randn/jig;MM(MM>1)=1;MM(MM<0)=-MM;
    h=plot(600,MM,mk,'MarkerSize',ms,'Color',a{i});set(h,'MarkerEdgeColor',mke,'MarkerFaceColor',a{i});
end
set(gca,'Fontsize',15)

%% Step 5: Get a continuous trajactory for each type of SNP, by Fourier curve fitting
FOOS={};
for i = 1:7;
    XX=[SNP_trajectories{i}(1,:)'];YY=[SNP_trajectories{i}(2,:)'];
    foo = fit(XX,YY,'fourier8');FOOS{i}=foo;
end

%% Step 9: Compute the frequencies of each sublineages
Y=[];Days=0:539;
for day = 0:539;
    f1 = max(FOOS{1}(day),0);f2 = max(FOOS{2}(day),0);f3 = max(FOOS{3}(day),0);f4 = max(FOOS{4}(day),0);
    f5 = max(FOOS{5}(day),0);f6 = max(FOOS{6}(day),0);f7 = max(FOOS{7}(day),0);
    % Compute the frequency 
    f1 = f1-f2;f3=f3-f4-f5-f6;f6=f6-f1; f1=f1-f7; 
    if f3<0; f3=0;end;if f1<0;f1=0;end;if f4<0;f4=0;end;if f6<0;f6=0.0;end;if f7<0;f7=0;end;if f5<0;f5=0;end;
    if f6<0 & day<490; f6=0.005;end
    if day>180; f3=0;end;   % We believe it's noise if it's positive
    if day<150; f7=0;end;   % Positive values are fitting noises
    if day>400; f4=0;f5=0;f6=0;f1=0;f3=0;end    % Positive values are believed to be noises
    Y=[Y;[f2 f4 f5 f7 f1 f6 f3]];
end

%% Eliminate the discontinuity
FOOS2={};Days2=1:539;
for i = 1:7;
    foo = fit(Days',Y(:,i),'fourier8');
    FOOS2{i}=foo;
end


%% Plot the muller plot by bar-stacking 
Yxz=[];Yx=[];Days=0:539;
LW=0.004;DLW=8;
for day = 0:1:539;
    f2 = FOOS2{1}(day);f4 = FOOS2{2}(day);f5 = FOOS2{3}(day);
    f7 = FOOS2{4}(day);f1 = FOOS2{5}(day);f6 = FOOS2{6}(day);f3 = FOOS2{7}(day);
    total = f1+f2+f3+f4+f5+f6+f7;
    f2=f2/total;f3=f3/total;f4=f4/total;f5=f5/total;f6=f6/total;f7=f7/total;f1=f1/total;
    % Separate for plotting
    f1o1 = (0.5*(450-day)/400)*f1; f1o2=f1-f1o1;

    Yx=[Yx;[f2 f4 f5 f1o1 f7 f1o2 f6 f3]];
    if mod(day,DLW*2)<DLW;
        if f4>LW;
            Yxz=[Yxz;[f2-LW LW*2 f4-LW f5 f1o1 f7 f1o2 f6 f3]];
        else if f4+f5>LW;
                Yxz=[Yxz;[f2-LW LW*2 0 f5+f4-LW f1o1 f7 f1o2 f6 f3]];
            else if f4+f5+f1o1>LW;
                    Yxz=[Yxz;[f2-LW LW*2 0 0 f1o1+f4+f5-LW f7 f1o2 f6 f3]];
                else;
                    Yxz=[Yxz;[f2-LW LW*2 0 0 0 f1o1+f4+f5+f7-LW f1o2 f6 f3]];
                end
            end
        end
    else
        Yxz=[Yxz;[f2 0 f4 f5 f1o1 f7 f1o2 f6 f3]];
    end
end

figure(2);
b=bar(Yxz,'stacked','Barwidth',1);
c = [a{2};cx;a{4};a{5};a{1};a{7};a{1};a{6};a{3}];
xlim([0 539])
ylim([0 1])
b(1).Parent.Parent.Colormap = c;
box off
set(gca,'xtick',[])
set(gca,'ytick',[])
% 


%% Step 10: Plot historical inference
h2=Yx(1,1);h4=Yx(1,2);h5=Yx(1,3);h1=Yx(1,4)+Yx(1,6);h6=Yx(1,7);h3=Yx(1,8); h3t=1-h2;
h61=h6+h1;
m61=0.6;s61=0.05;
m1=0.8;s1=0.02;
m5=0.7;s5=0.02;
m4=0.7;s4=0.02;
m3=0.5;s3=0.05;
m2=0.5;s2=0.05;
Ypre=[]
%%
for pt = 0:0.001:1;
    f2 = h2*exp((pt-m2)/s2)/(1+exp((pt-m2)/s2)); f021=(h2-f2)/2; f022=(h2-f2)/2;
    f3 = h3t*exp((pt-m3)/s3)/(1+exp((pt-m3)/s3)); f031=(h3t-f3)/2; f032=(h3t-f3)/2;
    Ypre=[Ypre;[0 0 0 0 0 0 0 0 0 0 0 0 0 0 f021 f2 f022 f031 f3 f032]];% 0 0 0 0 0 0 0 0 ]];
end

%%
DLW=0.04;
for pt = 0.3:0.0002:1;
    f2 = h2;
    f4 = h4*exp((pt-m4)/s4)/(1+exp((pt-m4)/s4)); f341=(h4-f4)/2; f342=(h4-f4)/2;
    f5 = h5*exp((pt-m5)/s5)/(1+exp((pt-m5)/s5)); f351=(h5-f5)/2; f352=(h5-f5)/2;
    f61 = h61*exp((pt-m61)/s61)/(1+exp((pt-m61)/s61)); f361=(h61-f61)/2; f362=(h61-f61)/2;
    f1 = h1*exp((pt-m1)/s1)/(1+exp((pt-m1)/s1));
    f611 = 1*(f61-f1)*(1-pt); f612 = f61-f1-f611;
    f611 = max((0.5*h1-f361-0.5*f1),0);f612=f61-f1-f611;
    f3 = h3;
    if mod(pt*1000,DLW*2000)<DLW*1000;
        if f341>LW;
            Ypre = [Ypre;[f2-LW 2*LW f341-LW f4 f342 f351 f5 f352 f361 f611 f1 f612 f362 f3 0 0 0 0 0 0]];
        else
            Ypre = [Ypre;[f2-LW 2*LW 0 f4+f341-LW f342 f351 f5 f352 f361 f611 f1 f612 f362 f3 0 0 0 0 0 0]];
        end
    else
        Ypre = [Ypre;[f2 0 f341 f4 f342 f351 f5 f352 f361 f611 f1 f612 f362 f3 0 0 0 0 0 0]];
    end

end
%%
for pt = 1:10;
    Ypre = [Ypre; zeros(1,20)];
end

figure(3); 
b=bar([Ypre],'stacked','Barwidth',1);
c = [a{2};cx;a{3};a{4};a{3};a{3};a{5};a{3};a{3};a{6};a{1};a{6};a{3};a{3};a{9};a{2};a{9};a{9};a{3};a{9}];
b(1).Parent.Parent.Colormap = c;
xlim([0 4500])
ylim([0 1])
box off
set(gca,'xtick',[])
set(gca,'ytick',[])
