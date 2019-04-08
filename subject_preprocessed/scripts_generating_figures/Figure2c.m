%% Read all intermediary files
figure(101);hold on

L01ts=[];L01mdSNPs=[];
ts=[];mdSNPs=[];
KKK={'o','^','<','h','s','d','p'};

for i = 1;
    load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    for ltp = 2:length(t);
        L01ts = [L01ts t(ltp)/365];
        L01mdSNPs=[L01mdSNPs mean(dSNP{ltp})];
    end
end
shapes=[];
for i = 2:12;
    if i < 10;
        load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    else
        load(['L' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    end
    for ltp = 2:length(t);
        ts = [ts t(ltp)/365];
        shapes=[shapes i];
        mdSNPs=[mdSNPs mean(dSNP{ltp})];
    end
end


% Fitting!
x=[L01ts ts];y=[L01mdSNPs mdSNPs];
[B,S] = polyfit(x,y, 1);

p = polyfit(x, y, 1);
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal
p

% Plot the fitting curve first
C3=[1 1 1]*0.6;
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',2);

% for L01:
C1=[1 1 1]*0.35;
h=plot(L01ts,L01mdSNPs,'o','Color',C1,'MarkerSize',10);
set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
% for non-L01
C2=[1 1 1]*0.8;
for i = 1:length(ts);
    h=plot(ts(i),mdSNPs(i),KKK{shapes(i)},'Color',C2,'MarkerSize',10);
    set(h,'MarkerEdgeColor',C2,'MarkerFaceColor',C2);

end
xlim([0,2])
ylim([0 2])
set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});

set(gca,'Fontsize',10)
set( gca, 'TickDir', 'out' );
set(gca,'linewidth',1)


%% Plot LEGENDS
% figure(242);hold on;
% h=plot(0,8-1,KKK{1},'Color',C1,'MarkerSize',10);
% set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
% for i = 2:7;
%     h=plot(0,8-i,KKK{i},'Color',C2,'MarkerSize',10);
%     set(h,'MarkerEdgeColor',C2,'MarkerFaceColor',C2);
% end
% xlim([-0.1 0.1])
% ylim([-1 9])
