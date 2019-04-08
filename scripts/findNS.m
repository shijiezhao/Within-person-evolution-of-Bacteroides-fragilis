function [Ns,Ss,NS] = findNS(SEQ,dictN,dictS);
tl = length(SEQ);

%%
Ns=zeros(6,1);
Ss=zeros(6,1);
for i = 1:tl/3;
    codon = SEQ(i*3-2:i*3);
    Ns = Ns + dictN(codon);
    Ss = Ss + dictS(codon);
end
NS = Ns./(Ns+Ss);

if rem(tl,3)>0;
    warning('Not 3x nucleotides!!!')
    Ns=zeros(6,1);
    Ss=zeros(6,1);
end
end