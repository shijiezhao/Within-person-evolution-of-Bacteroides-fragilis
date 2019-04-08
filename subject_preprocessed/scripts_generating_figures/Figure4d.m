%% 0. Load the mutations that are divergent across donors
load('summary_data/allgenes.across.donors.mat')
% load('all.mat');
%% DEFINE THE COLORS
Red1=[244 113 126]/256; Red2=[212 160 156]/256;Blue1=[224 236 234]/256; 
Blue2=[193 219 216]/256; Blue3=[100 142 180]/256;Grey1=[64 64 64]/256; Grey2=[200 200 200]/256;
c1 = Grey1; c2=Grey2;c3=Blue2;c4=Red2; mycolors=([c1;c2;c3;c4]);
%% Plot the cross-lineage mutations
figure(1);hold on;
set(gca,'TickDir','out')
N=0;S=0;    % Numbers of observed Non-sysnonymous mutation and Synonymous mutation
for i = 1:length(annotation_full_all);
    if annotation_full_all(i).type=='N';
        N=N+1;
    else if annotation_full_all(i).type=='S';
        S=S+1;
        end
    end
end

%% Plot the bar with CI for cross-lineage mutations
% d = expectedNS(annotation_full_all)  %% This step takes long time, the
% result from this annotation_full_all datastructure is 2.1544
expected_dnds_across_lineages = 2.1544;
[a1,b1,c1]=binomialCIdNdS(N, S, expected_dnds_across_lineages) % 2.1544 is the expected N/S ratio of the cross-species 
h=bar(3,log(c1),0.5);
set(h,'FaceColor',mycolors(1,:));set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);
plot([3 3],[log(a1) log(b1)],'Color','black','LineWidth',1)

%% Format the figure
ylim([log(0.02) log(64)])
set(gca,'Ytick',[log(0.06) log(0.25) log(1) log(4) log(16) log(64)])
set(gca,'Yticklabel',{'0.06' '0.25'  '1' '4' '16'  '64'})

%% Plot SNPs from the 16 that are under selection
wi_ind=[];
for i = 1:length(annotation_full);
    if length(annotation_full(i).locustag)>5;
        tmp = strsplit(annotation_full(i).locustag,'_')
        tmpn = str2num(tmp{2});
        if sum(tmpn==within_ltg)>0;  %% find mutations from genes with parallel mutation within individuals
            wi_ind=[wi_ind;i];
        end
    end
end
Subwi=annotation_full(wi_ind);      %% subset of SNPs from genes with parallel mutation within individuals
NS_wi=expectedNS(Subwi);
%
N6=0;S6=0;
for i = 1:length(Subwi);
    if Subwi(i).type=='N';
        N6=N6+1;
    else if Subwi(i).type=='S';
        S6=S6+1;
        end
    end
end
[a2,b2,c2]=binomialCIdNdS(N6,S6,NS_wi);
h=bar(1,log(c2),0.5);
set(h,'FaceColor',mycolors(4,:));
set(h,'LineWidth',0.05);
set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);
plot([1 1],[log(a2) log(b2)],'Color','black','LineWidth',1)


%% Plot SNPs from genes other than the 16 that are under selection
si_ind=[];
for i = 1:length(annotation_full);
    if sum(i==wi_ind)==0;% & sum(i==m2_ind)==0 & sum(i==m3_ind)==0 & sum(i==M2_ind)==0 & sum(i==M3_ind)==0;
    	si_ind=[si_ind;i];
	end
end
Subsi=annotation_full(si_ind);
NS_si=expectedNS(Subsi);
%
Nsi=0;Ssi=0;
for i = 1:length(Subsi);
    if Subsi(i).type=='N';
        Nsi=Nsi+1;
    else if Subsi(i).type=='S';
        Ssi=Ssi+1;
        end
    end
end
[a3,b3,c3]=binomialCIdNdS(Nsi,Ssi,NS_si);
h=bar(2,log(c3),0.5);
set(h,'FaceColor',mycolors(2,:));
set(h,'LineWidth',0.05);
set(h,'EdgeAlpha',0);
h(1).BaseValue = log(1);
plot([2 2],[log(a3) log(b3)],'Color','black','LineWidth',1)