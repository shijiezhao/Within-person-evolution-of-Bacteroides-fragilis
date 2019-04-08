function [allthresholdsmet] = div_single_sample_test_thresholds(d, p, controlMAF, coveragethresholds, majorNT, minorNT)


%Changed December 2012 to require forward and reverse strand qualities
%Changed later in December 2012 to make coverage a percentile based thing


%Very similar but slightly different (single sample) then
%div_test_tresholds

num_reqs=12;


NPositions=size(d,2);

%unpack parameters
minorfreqthreshold=p.minorfreqthreshold;
minreads_perstrand=p.minreads_perstrand;
maxreads_perstrand=coveragethresholds(p.maxreads_perstrand_percentile);
min_bq=p.min_bq;
min_mq=p.min_mq;
min_td=p.min_td;
max_td=p.max_td;
max_sbp=p.max_sbp;
max_bqp=p.max_bqp;
max_mqp=p.max_mqp;
max_tdp=p.max_tdp;
max_percent_indels=p.max_percent_indels;
max_percent_ends=p.max_percent_ends;
min_control_MAF=p.min_control_MAF;
minreads_perstrand_per_allele=p.minreads_perstrand_per_allele;

if nargin<5
    %unpack relevant stats
    [~, majorNT, minorNT] = div_major_allele_freq(d);
end

allthresholdsmet=zeros(NPositions,1);
positionsv=(1:NPositions)';

n1=majorNT';
n2=minorNT';

minorfreq=double((d(sub2ind(size(d),n2,positionsv))+d(sub2ind(size(d),n2+4,positionsv))))'./sum(d(1:8,:));
readsf=sum(d(1:4,:));
readsr=sum(d(5:8,:));
f1 = d(sub2ind(size(d),n1,positionsv)); %major allele counts on forward strand
r1 = d(sub2ind(size(d),n1+4,positionsv)); %major allele counts on forward strand
f2 = d(sub2ind(size(d),n2,positionsv)); %minor allele counts on forward strand
r2 = d(sub2ind(size(d),n2+4,positionsv)); %minor allele counts on forward strand
majorbqf=d(sub2ind(size(d),n1+8,positionsv));
majorbqr=d(sub2ind(size(d),n1+12,positionsv));
minorbqf=d(sub2ind(size(d),n2+8,positionsv));
minorbqr=d(sub2ind(size(d),n2+12,positionsv));
majormqf=d(sub2ind(size(d),n1+16,positionsv));
majormqr=d(sub2ind(size(d),n1+20,positionsv));
minormqf=d(sub2ind(size(d),n2+16,positionsv));
minormqr=d(sub2ind(size(d),n2+20,positionsv));
majortdf=(d(sub2ind(size(d),n1+24,positionsv)));
minortdf=(d(sub2ind(size(d),n2+24,positionsv)));
majortdr=(d(sub2ind(size(d),n1+28,positionsv)));
minortdr=(d(sub2ind(size(d),n2+28,positionsv)));
percent_indels=double(d(end,:))./(sum(d(1:8,:))+double(d(end,:)));
percent_ends=double(d(end-1,:))./sum(d(1:8,:));

SBp=d(end-6,:);
BQp=d(end-5,:);
MQp=d(end-4,:);
TDFp=d(end-3,:);
TDRp=d(end-2,:);



%these lines were used when quality scores for forward and reverse were
%averaged
%majorbq=((f1.*d(sub2ind(size(d),n1+8,positionsv)))+(r1.*d(sub2ind(size(d),n1+12,positionsv))))./(f1+r1);
%majormq=((f1.*d(sub2ind(size(d),n1+16,positionsv)))+(r1.*d(sub2ind(size(d),n1+20,positionsv))))./(f1+r1);
%minorbq=((f2.*d(sub2ind(size(d),n2+8,positionsv)))+(r2.*d(sub2ind(size(d),n2+12,positionsv))))./(f2+r2);
%minormq=((f2.*d(sub2ind(size(d),n2+16,positionsv)))+(r2.*d(sub2ind(size(d),n2+20,positionsv))))./(f2+r2);
%Tbq= ((majorbq > min_bq) & (minorbq > min_bq))';
%Tmq = ((majormq > min_mq) & (minormq > min_mq))';


%Find true/false of meeting thresholds
Tminor = minorfreq > minorfreqthreshold;
Treads= (readsf > minreads_perstrand) & (readsr > minreads_perstrand) &...
    ((readsf +readsr) < maxreads_perstrand) & (f2' > minreads_perstrand_per_allele) ...
    & (r2' > minreads_perstrand_per_allele);
Tbq= ((majorbqf > min_bq) & (minorbqf > min_bq) & (majorbqr > min_bq) & (minorbqr > min_bq))';
Tmq = ((majormqf > min_mq) & (minormqf > min_mq) & (majormqr > min_mq) & (minormqr > min_mq))';
Ttd = ((majortdf > min_td) & (majortdf < max_td) & (majortdr < max_td) & (majortdr > min_td)...
    & (minortdf > min_td) & (minortdf < max_td) & (minortdr > min_td) & (minortdr < max_td))';
Tid = percent_indels < max_percent_indels;
Te = percent_ends < max_percent_ends;


TSBp = SBp < max_sbp;
TBQp = BQp < max_bqp;
TMQp = MQp < max_mqp;
TTDp = (TDFp < max_tdp) & (TDRp < max_tdp);

if numel(controlMAF > 1)
    control = controlMAF > min_control_MAF;
else
    control = ones(NPositions,1);
end


allreqs= Tminor + Treads + Tbq + Tmq + Ttd + Tid + Te + TSBp + TMQp+ TBQp +TTDp+ control; % TMQp 


allthresholdsmet((allreqs==num_reqs))=1;



end


