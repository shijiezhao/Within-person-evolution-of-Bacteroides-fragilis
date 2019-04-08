function [LTG] =get_locustags_assembly(scaf,pos,CDS)


LTG = [];
for i = 1:length(scaf);
    mk = 0;
    s = scaf(i); p = pos(i);
    for j = 1:length(CDS{s});
        tmp = CDS{s}(j).indices-p;
        if tmp(1)*tmp(2)<0;
            mk = mk+1;
            if mk == 1;
                ltg = strsplit(CDS{s}(j).locustag,'_');
                ltg = str2double(ltg(2));
    %             ltg = str2num(ltg(2))

                LTG = [LTG; [ltg, abs(CDS{s}(j).indices(2)-CDS{s}(j).indices(1))]];
            end
        end
    end
    if mk == 0;
        LTG = [LTG; [0 0]];
    end
end

% %remove intragenic regions
% intergenic=cgenes~=floor(cgenes);
% intragenic=find(~intergenic);
% cgenes(intergenic)=[];
% 
% ltagnumbers=div_get_gene_numbers(codingsequences{1}(cgenes));
% 
% intergenic(intragenic(ltagnumbers==0))=1;  %mostly for RNA
% ltagnumbers=ltagnumbers(ltagnumbers>0);


end

