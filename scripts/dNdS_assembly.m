Names = {'D14','D26','D44','D52','D55','D57','D66','D77','D97','D128','D131','D440_normal','D44_hyper'};

A=[]

for i = 1:11;
name = i;
Donor = Names{name};
Donor
% workingdir=char(pwd);   % set up current directory as working dir.
% cd('../../..')                % move up two directories
% masterdir=char(pwd);    % set up the parent dir. as masterdir
% REFGENOMEFOLDER=[masterdir '/Reference_Genomes/' Donor '_SP201708'];  
% SCRIPTSDIRECTORY = [masterdir '/scripts'];
% path(SCRIPTSDIRECTORY,path);    % Add the srcipt folder to the search path
% cd(workingdir)
% load([REFGENOMEFOLDER '/cds_sorted.mat'])
% [ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);

%
% clear all;
K={'D14prokka.mat','D26prokka.mat','D44prokka.mat','D52prokka.mat','D55prokka.mat','D57prokka.mat','D66prokka.mat','D77prokka.mat','D97prokka.mat','D128prokka.mat','D131prokka.mat','D440prokka_normal.mat','D440prokka_hyper.mat'};
% K={'D14.mat'}
%K={'D440Hyper.mat'}
%K={'D440Normal.mat'}
%load('allgenes.across.donors.mat')
KK=[14 26 44 52 55 57 66 77 97 128 131 4401 4402];
for i = name;
    load(K{i});
    G=annotation_full;
    for k = 1:numel(G);
        G(k).Subject=KK(i);
        if length(G(k).locustag)>0;
        tmp=strsplit(G(k).locustag,'_');
        G(k).locusnum = str2double(tmp(2));
        else;
            G(k).locusnum = nan;
        end
    end
    if ~ismember('locustag1',fieldnames(G))
        G(1).locustag1=[];
        G(1).locustag2=[];
        G(1).protein1=[];
        G(1).protein2=[];
        G(1).gene1=[];
        G(1).gene2=[];
        G(1).distance1=nan;
        G(1).distance2=nan;
    end

    A=[A G];
end
end
annotation_full = A;

%% Test dnds
LCSN = []
for i = 1:length(annotation_full);
    LCSN = [LCSN;annotation_full(i).locusnum];
end
[a,b] = unique(LCSN);
triple_mut = [];
for i = 1:length(a);
    if a(i)>0;
        if sum(LCSN==a(i))>2;
            triple_mut = [triple_mut; find(LCSN==a(i))];
        end
    end
end

triple_annot = annotation_full(triple_mut);
expectedNS(triple_annot)

%% How dNdS change from variable genome to core genome...
MovWin = {[];[];[];[];[];[];[];[];[];[];[];[];[]};
k=0;
for i = 1:length(annotation_full);
    if length(annotation_full(i).locustag)>0;
        i
        tmp1 = strsplit(annotation_full(i).locustag,'_');
        tmp2 = strsplit(tmp1{1},'R');
        if tmp2{1}=='C';
            MovWin{13} = [MovWin{13}; i];
        end
        if length(tmp2)>1;
        if length(tmp2{2})>0;
            cat = str2double(tmp2{2});
            MovWin{cat} = [MovWin{cat};i];
        end
        end
    end
end
%%
figure(440);hold on;
for i = 1:12;
    subS=[];
    for k = 1:i;
        subS=[subS;MovWin{k}];
    end
    subSgene = annotation_full(subS);
    eNS = expectedNS(subSgene);
    N=0;S=0;
    for j = 1:length(subSgene);
        if subSgene(j).type=='S';
            S=S+1;
        end
        if subSgene(j).type=='N';
            N=N+1;
        end
    end
    dNdS=(N/S)/eNS;
    plot(i,dNdS,'O');
    plot(i,N/S,'O');
    r = eNS/(1+eNS);
    sims=[];
    for sim = 1:1000;
        tmp=0;n=0;s=0;
        for j = 1:(N+S);
            if rand(1)<r;
                n=n+1;
            else
                s=s+1;
            end
        end
        sims = [sims;(n/s)/eNS];
    end
    plot([i i],[quantile(sims,0.025),quantile(sims,0.975)],'-')

end
%%
figure(440);hold on;
for i = 13;
    subS=[];
    for k = i;
        subS=[subS;MovWin{k}];
    end
    subSgene = annotation_full(subS);
    eNS = expectedNS(subSgene);
    N=0;S=0;
    for j = 1:length(subSgene);
        if subSgene(j).type=='S';
            S=S+1;
        end
        if subSgene(j).type=='N';
            N=N+1;
        end
    end
    dNdS=(N/S)/eNS;
    plot(i,dNdS,'O');
    r = eNS/(1+eNS);
    sims=[];
    for sim = 1:1000;
        tmp=0;n=0;s=0;
        for j = 1:(N+S);
            if rand(1)<r;
                n=n+1;
            else
                s=s+1;
            end
        end
        sims = [sims;(n/s)/eNS];
    end
    plot([i i],[quantile(sims,0.025),quantile(sims,0.975)],'-')
    plot(i,N/S,'O');

end





