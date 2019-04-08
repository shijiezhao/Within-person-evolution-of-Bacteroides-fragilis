% Read tree file, and format the tree
tr = phytreeread(treename);
phytreeviewer(treename)
nwk=getnewickstr(tr);
piece = strsplit(nwk,':');
% Create heatmap
MEpattern=[];
NamesonTree={};
j=0;

CoveragePattern={}; markers=1;
%%
if strcmp(Donor,'D44');

    for i = 1:length(piece);
        if findstr('S1_',piece{i});
            j=j+1
            strpos=findstr('S1_',piece{i});
            NamesonTree{j}=piece{i}(strpos:length(piece{i}));
            sampleid=strmatch(NamesonTree{j},AllSampNames);
            sampleid
            MEpattern=[MEpattern TotalPattern(:,sampleid(1))];
            CoveragePattern{markers} = {};%% Markers: which donor it is referring to
            for tty = 1:length(TotalTYPE); %% tty: which linkage group 
                CoveragePattern{markers}{tty} = []; 
                for cont = 1:size(TotalTYPE{tty},1); %% cont: which specific mobile element region
                    start = TotalTYPE{tty}(cont,1); eend = TotalTYPE{tty}(cont,2);
                    covv = double(all_coverage_per_bp(sampleid,start:eend))/meancoverage(sampleid);
                    CoveragePattern{markers}{tty} = [CoveragePattern{markers}{tty} covv];
                end
            end
            markers=markers+1;  


        end
        if findstr('D2_Bact',piece{i});
            j=j+1
            strpos=findstr('D2_Bact',piece{i});
            NamesonTree{j}=piece{i}(strpos:length(piece{i}));
            sampleid=strmatch(NamesonTree{j},AllSampNames);
            sampleid
            MEpattern=[MEpattern TotalPattern(:,sampleid(1))];
            CoveragePattern{markers} = {};%% Markers: which donor it is referring to
            for tty = 1:length(TotalTYPE); %% tty: which linkage group 
                CoveragePattern{markers}{tty} = []; 
                for cont = 1:size(TotalTYPE{tty},1); %% cont: which specific mobile element region
                    start = TotalTYPE{tty}(cont,1); eend = TotalTYPE{tty}(cont,2);
                    covv = double(all_coverage_per_bp(sampleid,start:eend))/meancoverage(sampleid);
                    CoveragePattern{markers}{tty} = [CoveragePattern{markers}{tty} covv];
                end
            end
            markers=markers+1;  
        end
        if findstr('Bfrag',piece{i});
            j=j+1
            strpos=findstr('Bfrag',piece{i});
            NamesonTree{j}=piece{i}(strpos:length(piece{i}));
            sampleid=strmatch(NamesonTree{j},AllSampNames);
            sampleid
            MEpattern=[MEpattern TotalPattern(:,sampleid(1))];
            CoveragePattern{markers} = {};%% Markers: which donor it is referring to
            for tty = 1:length(TotalTYPE); %% tty: which linkage group 
                CoveragePattern{markers}{tty} = []; 
                for cont = 1:size(TotalTYPE{tty},1); %% cont: which specific mobile element region
                    start = TotalTYPE{tty}(cont,1); eend = TotalTYPE{tty}(cont,2);
                    covv = double(all_coverage_per_bp(sampleid,start:eend))/meancoverage(sampleid);
                    CoveragePattern{markers}{tty} = [CoveragePattern{markers}{tty} covv];
                end
            end
            markers=markers+1;  
        end
    end

else
    for i = 1:length(piece);
        if findstr(Donor,piece{i});
            j=j+1
            strpos=findstr(Donor,piece{i});
            if findstr('D131N134',piece{i});
                NamesonTree{j}=[piece{i}(strpos:length(piece{i})-2) '_' piece{i}(length(piece{i})-1:length(piece{i}))];
            else
                NamesonTree{j}=piece{i}(strpos:length(piece{i}));
            end
            sampleid=strmatch(NamesonTree{j},AllSampNames);
            sampleid = sampleid(1);
%             MEpattern=[MEpattern TotalPattern(:,sampleid(1))];
            CoveragePattern{markers} = {};%% Markers: which donor it is referring to
            for tty = 1:length(TotalTYPE); %% tty: which linkage group 
                CoveragePattern{markers}{tty} = []; 
                for cont = 1:size(TotalTYPE{tty},1); %% cont: which specific mobile element region
                    start = TotalTYPE{tty}(cont,1); eend = TotalTYPE{tty}(cont,2);
                    covv = double(all_coverage_per_bp(sampleid,start:eend))/meancoverage(sampleid);
                    CoveragePattern{markers}{tty} = [CoveragePattern{markers}{tty} covv];
                end
            end
            markers=markers+1;  
        end

    end
end
%%
% hmap=HeatMap(flipud(MEpattern'),'Colormap',redbluecmap);
% hFig = plot(hmap);
% saveas(hFig,'MobileElementsInformation/Heatmap','pdf');
% csvwrite(['MobileElementsInformation/' Donor '.csv'],MEpattern)
%%
close all force;
figure(66); hold on;
yl = length(CoveragePattern); xl = length(CoveragePattern{1});
Covs={}; for j=1:xl; Covs{j}=[];end;
MaxCov={};
for i = 1:yl;  % sample
    for j = 1:xl; % mobile element
        if i==1; 
            MaxCov{j}=max(smooth(CoveragePattern{i}{j},length(CoveragePattern{i}{j})/100));
        else
            MaxCov{j}=max(MaxCov{j},max(smooth(CoveragePattern{i}{j},length(CoveragePattern{i}{j})/100)));
        end
    end
end

coveragess=[];
ylt=yl;
% yl=5;
for i = 1:yl;  % sample
    i
    ccc=[];
    for j = 1:xl; % mobile element
%         subplot(yl,xl,(i-1)*xl+j);
        subplot('Position',[(j-1+0.5)/(xl+1) 1-(i+0.5)/(yl+1) 1/(xl+1) 1/(yl+1)]);hold on;
        h=area(smooth(CoveragePattern{i}{j},length(CoveragePattern{i}{j})/30));
        h(1).FaceColor = (1-min(mean(CoveragePattern{i}{j})/5,1))*[1 1 1];
        h(1).EdgeColor = (1-min(mean(CoveragePattern{i}{j})/5,1))*[1 1 1];
        ylim([0 max(1,1.2*max(smooth(CoveragePattern{i}{j},length(CoveragePattern{i}{j})/100)))]);
        xlim([0 length(CoveragePattern{i}{j})]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(get(gca,'xlabel'),'color',[0 0 0]);
        set(get(gca,'ylabel'),'color',[0 0 0]);
%         box on;
        Covs{j} = [Covs{j} smooth(CoveragePattern{i}{j},100)];
        ccc = [ccc;mean(smooth(CoveragePattern{i}{j},length(CoveragePattern{i}{j})/30))];
    end
    coveragess=[coveragess ccc];
end

%%
for i=1:yl;
    for j=1:xl;
        subplot('Position',[(j-1+0.5)/(xl+1) 1-(i+0.5)/(yl+1) 1/(xl+1) 1/(yl+1)]);hold on;
        box on;
        set(gca,'LineWidth',1);
        set(gca,'layer','top')%% Important... set the axis on top
    end
end

saveas(gca,'Densitymap.pdf')