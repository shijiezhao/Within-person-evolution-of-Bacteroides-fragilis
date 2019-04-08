figure; hold on;
p1 = 4797546	;
p2 = p1+1681;
chunk = all_coverage_per_bp(:,p1:p2);
    for j = 1:size(all_coverage_per_bp,1);
        normchunk = double(chunk(j,:))/double(Meancoverage);
%         if goodsamples(j)==1;
            plot(normchunk);

        
    end