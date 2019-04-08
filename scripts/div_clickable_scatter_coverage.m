
function div_clickable_scatter_coverage(strains, starts, ends, isdeletion, SampleNames, coverage_matrix)
figure; clf ; hold on;




for i=1:numel(strains)
    if isdeletion(i)>0
        plot([starts(i) ends(i)], [strains(i) strains(i)],'Color', rgb('Green'),'LineWidth',5)
    else
        plot([starts(i) ends(i)], [strains(i) strains(i)],'Color', rgb('Blue'),'LineWidth',5)
    end
    
end


p1=plot(starts,strains,'o');
set(p1,'MarkerSize', 4, 'MarkerFaceColor', rgb('Red'), 'MarkerEdgeColor', rgb('Red'),'ButtonDownFcn',@clicked);
ylim([0 numel(SampleNames)+1])

x=starts;
y=strains;

xrange=max(ends)-min(starts);
yrange=numel(SampleNames);




    function clicked(src,event)
        

        
        
        %find data point clicked, this region written by Roy Kishony March
        %2012
        ac=get(gca,'position') ;
        fc=get(gcf,'position') ;
        pc=get(gcf,'CurrentPoint') ;
        xl=xlim ;
        yl=ylim ;
        ax = [fc(3) * [ac(1), ac(1)+ac(3)]]  ;
        ay = [fc(4) * [ac(2), ac(2)+ac(4)]]  ;
        x0 = xl(1) + (xl(2)-xl(1))*(pc(1)-ax(1))/(ax(2)-ax(1)) ;
        y0 = yl(1) + (yl(2)-yl(1))*(pc(2)-ay(1))/(ay(2)-ay(1)) ;
        
        
        
        
        [~,ind] = min(((x-x0)/xrange).^2+((y-y0)/yrange).^2) ;
        
        

        extraonsides=5000;
        
        disp(strains(ind))
        disp(SampleNames(strains(ind)))
        disp(starts(ind))
        disp(ends(ind)-starts(ind))
      
        plotstart=starts(ind)-extraonsides;
        plotend=ends(ind)+extraonsides;
        
        if(isdeletion(ind)>0)
            disp('deletion')
        else
            disp('duplication')
        end
        

        
        %Absolute coverage plot
         figure(11); clf; hold on;
        %This plots each strain individually
        for j=1:numel(SampleNames)
  
            plot(plotstart:plotend,coverage_matrix(j,plotstart:plotend), 'LineWidth',.1,'Color', (.5*j/numel(SampleNames))*[1 1 1]) 
        end
        %This highlights the strain you clicked on
        plot(plotstart:plotend,coverage_matrix(strains(ind),plotstart:plotend),'LineWidth',10, 'Color', rgb('Red'))
        
        %This plots a control strain (randomly chosen, currently #2)
        plot(plotstart:plotend,coverage_matrix(2,plotstart:plotend),'LineWidth',5, 'Color', rgb('Black'))

        xlim([plotstart plotend])
        ylabel('Unnormalized coverage')
        xlabel('Position on genome')
        temp=ylim;
        plot([starts(ind) starts(ind)], [temp(1) temp(2)],'Color', rgb('Blue'));
        plot([ends(ind) ends(ind)],[temp(1) temp(2)], 'Color', rgb('Blue'));
        
  
    end

end
