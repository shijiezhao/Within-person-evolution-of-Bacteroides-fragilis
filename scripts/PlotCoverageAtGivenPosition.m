Chr = 15; Pos = 561;
P = ChrStarts(Chr)+Pos
chunk = all_coverage_per_bp(:,P-50:P+50);
close all force;
for j = 1:size(chunk,1)-12;
        normchunk = double(chunk(j,:))/double(meancoverage(j));
        if goodsamples(j)==1;
            plot(normchunk); hold on;
            sample=AllSampNames(j);
            coverage = mean(normchunk);
        end
end
    