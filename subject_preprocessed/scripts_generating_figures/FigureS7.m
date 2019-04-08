%% Plot the collector curves
load('summary_data/Colors.mat') % load color scheme
SNPs=[];dMRCATO=[];
figure(1);hold on;
jig=0.05;
TotIsolates=0;
for i = 1:12;
    if i < 10;
        load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    else
        load(['L' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    end
    figure(1);
    %Collector curve for dMRCA
    subplot(6,4,i);hold on;
    set(gca,'TickDir','out');
    ml=0;yml=0; % 
    for j = 1:length(dMRCA_collector);
        for k = 1:length(dMRCA_collector{j});
            scatter(k+jig*j,mean(dMRCA_collector{j}{k}),7,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',Colors{i}(j,:));
            plot([k+jig*j k+jig*j],[mean(dMRCA_collector{j}{k})-std(dMRCA_collector{j}{k}) mean(dMRCA_collector{j}{k})+std(dMRCA_collector{j}{k})],'-','Color',Colors{i}(j,:));
            if ml<k;ml=k;end; 
            if yml<mean(dMRCA_collector{j}{k})+std(dMRCA_collector{j}{k});
                yml=mean(dMRCA_collector{j}{k})+std(dMRCA_collector{j}{k});
            end
        end
    end
    
    xlim([0 (floor((ml-1)/5)+1)*5])
    ylim([0 yml*1.05])
    
    % Collector curve for SNPs
    subplot(6,4,i+12);hold on;
    set(gca,'TickDir','out');
    ml=0;yml=0;
    for j = 1:length(SNP_collector);
        for k = 1:length(SNP_collector{j});
            scatter(k+jig*j,mean(SNP_collector{j}{k}),7,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',Colors{i}(j,:));
            plot([k+jig*j k+jig*j],[mean(SNP_collector{j}{k})-std(SNP_collector{j}{k}) mean(SNP_collector{j}{k})+std(SNP_collector{j}{k})],'-','Color',Colors{i}(j,:));
            if ml<k;ml=k;end; 
            if yml<mean(SNP_collector{j}{k})+std(SNP_collector{j}{k});
                yml=mean(SNP_collector{j}{k})+std(SNP_collector{j}{k});
            end
        end
    end
    xlim([0 (floor((ml-1)/5)+1)*5])
    ylim([0 yml*1.05])

end
