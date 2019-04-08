function expectedNS_spectrum = expectedNS_CDS(CDS);
%% define categories
c=containers.Map; 
c('AT')=1;c('TA')=1;c('AC')=2;c('TG')=2;c('GC')=3;c('CG')=3;
c('GT')=4;c('CA')=4;c('AG')=5;c('TC')=5;c('GA')=6;c('CT')=6;

[dictN,dictS] = codon_dict();


%%
% [ChrStarts, glength, ~, ~, sequences]= genomestats(GenomeDir);
%%
Nspectrum=zeros(6,1); Sspectrum=zeros(6,1);
for i = 1:length(CDS);
    for j = 1:length(CDS{i});
        SEQ = CDS{i}(j).Sequence;
        if length(SEQ)>0;
            [Ns,Ss] = findNS(SEQ,dictN,dictS);
            Nspectrum = Nspectrum + Ns;
            Sspectrum = Sspectrum + Ss;
        end
        
        if rem(j,100)==0;
            j
        end
    end
end
expectedNS_spectrum=Nspectrum./(Nspectrum+Sspectrum);
end
        