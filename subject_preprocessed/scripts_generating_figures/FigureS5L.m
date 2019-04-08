%% Import phage plaque assay results (same with Table S6)
%% Phage_01, Phage_02, Phage_03: results from 3 different media
%% Donors: donor isolates with ME01-02 or not
%% Recipient: recipient isolate with ME01-02 or not
load('competition_exp/Phage_plaque_results.mat')
DATA={}
for d = 1:3;
    for r = 1:2;
        DATA{d,r}={[],[],[],[]};
    end
end
PHAGE={Phage_01,Phage_02,Phage_03};
%%
for d = 1:41;
    for r = 1:25;
        if Donors(d)>=0 & Recipients(r)>=0;
            tmp = []
            for m = 1:3;
                if PHAGE{m}(r,d)>=0;
                    xx = PHAGE{m}(r,d);
                    DATA{Donors(d)+1,Recipients(r)+1}{m}=[DATA{Donors(d)+1,Recipients(r)+1}{m};xx];
                    tmp = [tmp; PHAGE{m}(r,d)];
                end
            end
            if length(tmp)>0;
                DATA{Donors(d)+1,Recipients(r)+1}{4}=[DATA{Donors(d)+1,Recipients(r)+1}{4};sum(tmp)];
            end
        end
    end
end

%%
jit=0.1
figure(32);hold on
mk = 1;
h = zeros(3, 1);
YS={};
C1=[.1 .4 .9];
C2=[.9 .1 .22];
C3=[.42 .9 .1];
SC=20
for d = 1:3;
    for r = 1:2;
        Y1 = DATA{d,r}{1};
        Y2 = DATA{d,r}{2};
        Y3 = DATA{d,r}{3};
        if mk>5; mk=5;end
        h(1)=scatter(mk+randn(size(Y1))*jit,Y1+randn(size(Y1))*jit,SC,'MarkerEdgeAlpha',.5,'MarkerEdgeColor',C1,'MarkerFaceAlpha',.5,'MarkerFaceColor',C1)
        h(2)=scatter(mk+randn(size(Y2))*jit,Y2+randn(size(Y2))*jit,SC,'MarkerEdgeAlpha',.5,'MarkerEdgeColor',C2,'MarkerFaceAlpha',.5,'MarkerFaceColor',C2)
        h(3)=scatter(mk+randn(size(Y3))*jit,Y3+randn(size(Y3))*jit,SC,'MarkerEdgeAlpha',.5,'MarkerEdgeColor',C3,'MarkerFaceAlpha',.5,'MarkerFaceColor',C3)
%         plot(mk+randn(size(Y4))*jit,Y4+randn(size(Y4))*jit,'oc')
        YS{d,r}=[Y1;Y2;Y3]
        mk=mk+1;
    end
end

xticks([])
% xticklabels({'D-,R-','D-,R+','D+,R-','D+,R+','Negative Control'})
yticks([0 50 100 150 200])
yticklabels({'0','50','100','150','200'})
set(gca,'FontSize',15);
set(gca,'LineWidth',1);
xlim([0.2 5.8])
ylim([0 170])
box on;
legend(h,'BHI','BPRM','BPRM+bile');
% xlabel('Donor-Recipient combinations')
% ylabel('Number of plaques')
YS{3,1}=[YS{3,1};YS{3,2}];
XS={YS{1,1},YS{1,2},YS{2,1},YS{2,2},YS{3,1}};
%%
P1=zeros(5,5); % between D+R- and others
for i = 1:5;
    for j = 1:5;
        [a b]=ranksum(XS{i},XS{j});
        P1(i,j)=a;

    end
end


