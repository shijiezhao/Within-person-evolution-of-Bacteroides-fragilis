function show_coverage_window(start_location, coverageWindows,  names, sampleIndex)

%%Tami Lieberman, updated 2016 from a much older script
%Designed for plotting a region when table is clicked

%%
controlSample=1;

%initialize
figure(20); clf; hold on;
otherSamples=1:size(coverageWindows,1);
otherSamples(controlSample)=[];


l=[start_location:(start_location+size(coverageWindows,2)-1)];

for i=1:numel(otherSamples);
    plot(l,coverageWindows(otherSamples(i),:), '-', 'LineWidth', 1, 'Color', [i i i]*(.5/numel(otherSamples)) ); %plot other samples in black
end
c = plot(l,coverageWindows(controlSample,:)','k-', 'LineWidth',2); %plot control in black, wide


legendlist=[c]; legendnames={names{controlSample}};
if sampleIndex > 0
    s = plot(l,coverageWindows(sampleIndex,:),'r-', 'LineWidth',2); %plot clicked in red
    otherSamples(otherSamples==sampleIndex)=[];
    legendlist(end+1)=s;
    legendnames{end+1}=names{sampleIndex};
end

%legend, etc
xlim([l(1) l(end)])
ylim([0 inf])
ylabel('Coverage')
legend(legendlist, legendnames,'Location', 'BestOutside')


end

