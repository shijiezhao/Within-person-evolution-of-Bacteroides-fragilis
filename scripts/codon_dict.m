function [dictN,dictS] = codon_dict();
%% Build a dictionary
c=containers.Map; 
c('at')=1;c('ta')=1;c('ac')=2;c('tg')=2;c('gc')=3;c('cg')=3;
c('gt')=4;c('ca')=4;c('ag')=5;c('tc')=5;c('ga')=6;c('ct')=6;
c('AT')=1;c('TA')=1;c('AC')=2;c('TG')=2;c('GC')=3;c('CG')=3;
c('GT')=4;c('CA')=4;c('AG')=5;c('TC')=5;c('GA')=6;c('CT')=6;
dictN = containers.Map;
dictS = containers.Map;
NTs = 'atcg';
Ns = zeros(4,4);
Ss = zeros(4,4);
for n1 = 1:4;
    for n2 = 1:4;
        for n3 = 1:4;
            codon = NTs([n1,n2,n3]);
            aa1= double(nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11))
            N1 = [0;0;0;0;0;0]; % 1: at->ta; 2: at->cg; 3: gc->cg; 4: gc->ta; 5: at->gc; 6: gc->at;
            S1 = [0;0;0;0;0;0];
            for mut1 = 1:4;
                if mut1 ~= n1;
                    chemical_cat = c(NTs([n1,mut1]));
                    aa2 = double(nt2aa(NTs([mut1,n2,n3]),'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11));
                    if aa2 == aa1;
                        S1(chemical_cat) = S1(chemical_cat)+1;
%                         Ss(n1,mut1) = Ss(n1,mut1) + 1;
                    else;
                        N1(chemical_cat) = N1(chemical_cat)+1;
%                         Ns(n1,mut1) = Ns(n1,mut1) + 1;
                    end
                end
            end
            N2 = [0;0;0;0;0;0]; % 1: at->ta; 2: at->cg; 3: gc->cg; 4: gc->ta; 5: at->gc; 6: gc->at;
            S2 = [0;0;0;0;0;0];
            for mut2 = 1:4;
                if mut2 ~= n2;
                    chemical_cat = c(NTs([n2,mut2]));
                    aa2 = double(nt2aa(NTs([n1,mut2,n3]),'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11));
                    if aa2 == aa1;
                        S2(chemical_cat) = S2(chemical_cat)+1;
                    else;
                        N2(chemical_cat) = N2(chemical_cat)+1;
                    end
                end
            end
            N3 = [0;0;0;0;0;0]; % 1: at->ta; 2: at->cg; 3: gc->cg; 4: gc->ta; 5: at->gc; 6: gc->at;
            S3 = [0;0;0;0;0;0];
            for mut3 = 1:4;
                if mut3 ~= n3;
                    chemical_cat = c(NTs([n3,mut3]));
                    aa2 = double(nt2aa(NTs([n1,n2,mut3]),'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11));
                    if aa2 == aa1;
                        S3(chemical_cat) = S3(chemical_cat)+1;
                    else;
                        N3(chemical_cat) = N3(chemical_cat)+1;
                    end
                end
            end
            dictN(codon) = N1+N2+N3;%[N1;N2;N3];
            dictS(codon) = S1+S2+S3;%[S1;S2;S3];
        end
    end
    
end
end