function div_bar_charts(c, names, showall)
%showall is an optional argument to display more of counts matrix


LABEL_SIZE= 8; 


figure(660);clf; hold on;

if nargin < 3
    showall=0;
    num_call = subplot(1,1,1);
else
    num_call = subplot(4,1,1);
end


%Plot number of calls
a=squeeze(c(1:8,:));
bar(reshape([a; nan(4,size(a,2))],4,[])','stacked')
legend('A','T','C','G', 'Location', 'BestOutside')
ylabel('Number of reads')
set(num_call,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
set(num_call,'Xtick',2:3:(3*numel(names)-1))
xlim([0 (3*numel(names)+3)])


%Plot other options
if showall > 0
    %Plot call quality
    call_qual = subplot(4,1,2); hold on;
    title('Average call quality');
    a=squeeze(c(9:16,:));
    if nargin>3
        plot([0 3*size(a,2)], [params.min_bq params.min_bq], 'k:')
    end
    bar(reshape([a; nan(4,size(a,2))],4,[])','grouped', 'LineStyle', 'none')
    legend('Aq','Tq','Cq','Gq', 'Location', 'BestOutside')
    ylabel('Average Base Quality')
    set(call_qual,'Xtick',2:3:(3*numel(names)-1))
    set(call_qual,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
    xlim([0 (3*numel(names)+3)])
    set(call_qual, 'XTick', []);
    
    
    %Plot mapping quality
    map_qual = subplot(4,1,3); hold on;
    title('Average mapping quality');
    a=squeeze(c(17:24,:));
    if nargin>3
        plot([0 3*size(a,2)], [params.min_mq params.min_mq], 'k:')
    end
    bar(reshape([a; nan(4,size(a,2))],4,[])','grouped','LineStyle', 'none')
    legend('Am','Tm','Cm','Gm', 'Location', 'BestOutside')
    ylabel('Average Mapping Quality')
    set(map_qual,'Xtick',2:3:(3*numel(names)-1))
    set(map_qual,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
    xlim([0 (3*numel(names)+3)])
    set(map_qual, 'XTick', []);
    
    
    %Plot tail distance f
    tail_dist = subplot(4,1,4); hold on;
    title('Average tail distance');
    a=squeeze(c(25:32,:));
    if nargin>3
        plot([0 3*size(a,2)], [params.min_td params.min_td], 'k:')
        plot([0 3*size(a,2)], [params.max_td params.max_td], 'k:')
    end
    bar(reshape([a; nan(4,size(a,2))],4,[])','grouped', 'LineStyle', 'none')
    legend('Atd','Ttd','Ctd','Gtd', 'Location', 'BestOutside')
    ylabel('Average Tail Distance')
    set(tail_dist,'Xtick',2:3:(3*numel(names)-1))
    set(tail_dist,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
    xlim([0 (3*numel(names)+3)])
    set(tail_dist, 'XTick', []);
    
end

end


