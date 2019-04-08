function [final_dnds,sim2,x,y] = Prepare_dNdS(annotation_full,M2_ind);

SubM2=annotation_full(M2_ind);
NS_M2=expectedNS(SubM2);

%
N2=0;S2=0;
for i = 1:length(SubM2);
    if SubM2(i).type=='N';
        N2=N2+1;
    else if SubM2(i).type=='S';
        S2=S2+1;
        end
    end
end
sim2=[];
for i = 1:1000;
    r=N2/(N2+S2);
    n=0;
    for j = 1:N2+S2;
        if rand(1)<r;
            n=n+1;
        end
    end
    dNdS=(n/(N2+S2-n))/NS_M2;
    sim2=[sim2;dNdS];
end

final_dnds=(N2/S2)/NS_M2;

[x,y] = binomialCIdNdS(N2,S2,NS_M2);
end