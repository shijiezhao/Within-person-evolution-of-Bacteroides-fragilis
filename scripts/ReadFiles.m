
clear all;
A=[]
K={'D14prokka.mat','D26prokka.mat','D44prokka.mat','D52prokka.mat','D55prokka.mat','D57prokka.mat','D66prokka.mat','D77prokka.mat','D128prokka.mat','D131prokka.mat','D97prokka.mat','D440prokka.mat'};
% K={'D14.mat'}
%K={'D440Hyper.mat'}
%K={'D440Normal.mat'}
%load('allgenes.across.donors.mat')
KK=[14 26 44 52 55 57 66 77 97 128 131 440]
for i = 1:numel(K);
    load(K{i});
    G=annotation_full;
    for k = 1:numel(G);
        G(k).Subject=KK(i);
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
annotation_full = A;

666666666666666