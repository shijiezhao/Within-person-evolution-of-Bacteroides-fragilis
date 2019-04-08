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
annotation_11L = A;
%%
A=[];
for i = 12;
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
annotation_Norm = A;
%%
A=[];
for i = 13;
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
annotation_Hyper = A;


%% Compare core genome percentage
Mat=[0 0; 0 0; 0 0];
for i = 1:length(annotation_Hyper);
    if length(annotation_Hyper(i).locustag)>0;
    if annotation_Hyper(i).locustag(1) == 'C';
        Mat(1,1) = Mat(1,1) + 1;
    end
    if annotation_Hyper(i).locustag(1) == 'M';
        Mat(1,2) = Mat(1,2) + 1;
    end
    end
end
for i = 1:length(annotation_Norm);
    if length(annotation_Norm(i).locustag)>0;
    if annotation_Norm(i).locustag(1) == 'C';
        Mat(2,1) = Mat(2,1) + 1;
    end
    if annotation_Norm(i).locustag(1) == 'M';
        Mat(2,2) = Mat(2,2) + 1;
    end
    end
end
for i = 1:length(annotation_11L);
    if length(annotation_11L(i).locustag)>0;
    if annotation_11L(i).locustag(1) == 'C';
        Mat(3,1) = Mat(3,1) + 1;
    end
    if annotation_11L(i).locustag(1) == 'M';
        Mat(3,2) = Mat(3,2) + 1;
    end
    end
end




