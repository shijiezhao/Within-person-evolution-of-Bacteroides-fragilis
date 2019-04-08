function clickable_scatter_twocolor_oneoutput(x,y, colorred, output)

if sum(~colorred)>0
plot(x(~colorred),y(~colorred),'o','MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'ButtonDownFcn',@clicked) ;
end
if sum(colorred)>0
plot(x(colorred),y(colorred),'o','MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'ButtonDownFcn',@clicked) ;
end




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
        [~,ind] = min((x-x0).^2+(y-y0).^2) ;
        
        disp(output(ind))
        
        
    end

end
