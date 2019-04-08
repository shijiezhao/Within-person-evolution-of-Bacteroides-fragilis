function [mutAF, mutantNT, all_mutAF]=mutant_frequency(cnts,hasmut, anc, calls)

%hasmut is a matrix with true/false
%ancestor is a vector of 0-4 representing unknown,A,T,C,G




[maf, maNT, minorNT] = div_major_allele_freq(cnts);

minormutation=(hasmut & (repmat(anc,1,size(hasmut,2))==maNT));


minorAF=div_minor_allele_freq(cnts);
mutAF=zeros(size(hasmut));
mutAF(hasmut>0)=maf(hasmut>0);
mutAF(minormutation)=minorAF(minormutation);

%if insertion or deletion, counts is uniformative
mutAF(hasmut>0 & ismember(calls,'ID'))=1;

%mutantNT
mutantNT=zeros(size(hasmut));
mutantNT(hasmut)=maNT(hasmut); 
mutantNT(minormutation)=minorNT(minormutation);

%gives mutation allele frequency only for rows with at least one called
%mutation
all_mutAF=zeros(size(mutAF));
all_mutAF(repmat(anc,1,size(hasmut,2))==maNT & minorAF > 0 )=minorAF(repmat(anc,1,size(hasmut,2))==maNT & minorAF > 0);
all_mutAF(repmat(anc,1,size(hasmut,2))~=maNT & maf > 0 )=maf(repmat(anc,1,size(hasmut,2))~=maNT & maf > 0);


