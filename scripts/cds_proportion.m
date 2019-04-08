L=0
load([REFGENOMEFOLDER '/cds_sorted.mat'])
for i = 1:length(CDS);
    for j = 1:length(CDS{i});
        L = L + abs(CDS{i}(j).loc2-CDS{i}(j).loc1);
    end
end
G=0;I=0;
for i = 1:length(annotation_full);
    if annotation_full(i).type =='N' | annotation_full(i).type =='S';
        G=G+1;
    else
        I=I+1;
    end
end

r1 = L/size(all_coverage_per_bp,2)
r2 = G/(G+I)
sim=[];
for i = 1:1000;
    gs=0;
    for j = 1:length(annotation_full);
        tmp = rand(1);
        if tmp<r1;
            gs = gs+1;
        end
    end
    sim=[sim;gs/length(annotation_full)];
end
l1 = quantile(sim,0.05)
r1 = quantile(sim,0.95)


sim=[];
for i = 1:1000;
    gs=0;
    for j = 1:length(annotation_full);
        tmp = rand(1);
        if tmp<r2;
            gs = gs+1;
        end
    end
    sim=[sim;gs/length(annotation_full)];
end
l2 = quantile(sim,0.05)
r2 = quantile(sim,0.95)