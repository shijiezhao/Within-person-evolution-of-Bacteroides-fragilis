function [expectedNS,ct] = expectedNS(annotation_full);
%% define categories
c=containers.Map; 
c('AT')=1;c('TA')=1;c('AC')=2;c('TG')=2;c('GC')=3;c('CG')=3;
c('GT')=4;c('CA')=4;c('AG')=5;c('TC')=5;c('GA')=6;c('CT')=6;

[dictN,dictS] = codon_dict();

%%
indiv_expectedNratio=[];
ct=[];
for i = 1:size(annotation_full,2);
    if mod(i,1000)==0; 
        i
    end
    SEQ = annotation_full(i).Sequence;
    if length(SEQ)>0 & length(annotation_full(i).nts)==2;
%         codonrank = rem(annotation_full(i).nt_pos+1,3);
        [~,~,eNS_spectrum] = findNS(SEQ,dictN,dictS);
        ancnt = annotation_full(i).anc;
        nts = annotation_full(i).nts;
        newnt = setdiff(nts,ancnt);
        chemical_cat = c([ancnt newnt]);
        ct=[ct;chemical_cat];
        eNS = eNS_spectrum(chemical_cat);
        indiv_expectedNratio = [indiv_expectedNratio; eNS];
    end
end
%%
expectedNratio = mean(indiv_expectedNratio);
expectedNS = expectedNratio/(1-expectedNratio);
end