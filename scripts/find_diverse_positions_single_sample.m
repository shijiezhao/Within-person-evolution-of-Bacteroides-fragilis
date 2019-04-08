function [p, coveragethresholds, MAF_control, numfields] = find_diverse_positions_single_sample(SampleDir, SampleName, tempfolder)

load('for_finding_diverse_positions');


%unpack parameters
minorfreqthreshold=params.minorfreqthreshold;
%maxreads_perstrand=params.maxreads_perstrand;
minreads_perstrand=params.minreads_perstrand;
minreads_perstrand_per_allele=params.minreads_perstrand_per_allele;
min_bq=params.min_bq;
min_mq=params.min_mq;
max_td=params.max_td;
min_td=params.min_td;
max_sbp=params.max_sbp;
max_bqp=params.max_bqp;
%max_mqp=params.max_mqp;
max_tdp=params.max_tdp;
max_percent_indels=params.max_percent_indels;


num_reqs=9;


p=zeros(GenomeLength,1);


f=fopen(['find_diverse_positions_log/' SampleName '.txt'] ,'w');



coveragethresholds=zeros(100,1);



%load data
load([SampleDir '/diversity.mat']);


%coveragethresholds -- contains cdf cutoffs .01:.01:1.0.
%only counts positions where at least 1 read aligned
%purpose of this data structure is to allow downstream processes to
%remove positions with excess coverage in a way that accounts for
%variation in coverage between samples
cov=sum(double(data(1:8,:)));
cov(cov<1)=[];


if ~isempty(cov) & sum(cov)/length(data)> 0.5

    cov=sort(cov);
    cutoffs=.01:.01:1;
    for j=1:numel(cutoffs);
        coveragethresholds(j)=cov(floor(cutoffs(j)*numel(cov)));
    end
    
    
    
    
    
    %Parse data
    NPositions=size(data,2);
    [maf, majorNT, minorNT] = div_major_allele_freq(data);
    positionsv=(1:NPositions)';
    n1=majorNT';
    n2=minorNT';
    
    minorfreq=(double(data(sub2ind(size(data),n2,positionsv)))+double(data((sub2ind(size(data),n2+4,positionsv)))))'./sum(double(data(1:8,:)));
    readsf=sum(double(data(1:4,:)));
    readsr=sum(double(data(5:8,:)));
    f1 = data(sub2ind(size(data),n1,positionsv)); %major allele counts on forward strand
    r1 = data(sub2ind(size(data),n1+4,positionsv)); %major allele counts on forward strand
    f2 = data(sub2ind(size(data),n2,positionsv)); %major allele counts on forward strand
    r2 = data(sub2ind(size(data),n2+4,positionsv)); %major allele counts on forward strand
    majorbqf=data(sub2ind(size(data),n1+8,positionsv))';
    majorbqr=data(sub2ind(size(data),n1+12,positionsv))';
    minorbqf=data(sub2ind(size(data),n2+8,positionsv))';
    minorbqr=data(sub2ind(size(data),n2+12,positionsv))';
    majormqf=data(sub2ind(size(data),n1+16,positionsv))';
    majormqr=data(sub2ind(size(data),n1+20,positionsv))';
    minormqf=data(sub2ind(size(data),n2+16,positionsv))';
    minormqr=data(sub2ind(size(data),n2+20,positionsv))';
    
    %Changed December 2012 to require forward and reverse strand qualities
    %   majorbq=(((f1.*data(sub2ind(size(data),n1+8,positionsv)))+(r1.*data(sub2ind(size(data),n1+12,positionsv))))./(f1+r1))';
    % majormq=(((f1.*data(sub2ind(size(data),n1+16,positionsv)))+(r1.*data(sub2ind(size(data),n1+20,positionsv))))./(f1+r1))';
    %  minorbq=(((f2.*data(sub2ind(size(data),n2+8,positionsv)))+(r2.*data(sub2ind(size(data),n2+12,positionsv))))./(f2+r2))';
    % minormq=(((f2.*data(sub2ind(size(data),n2+16,positionsv)))+(r2.*data(sub2ind(size(data),n2+20,positionsv))))./(f2+r2))';
    majortdF=(data(sub2ind(size(data),n1+24,positionsv)))';
    minortdF=(data(sub2ind(size(data),n2+24,positionsv)))';
    majortdR=(data(sub2ind(size(data),n1+28,positionsv)))';
    minortdR=(data(sub2ind(size(data),n2+28,positionsv)))';
    percent_indels=double(data(end,:))./double(sum(data(1:8,:))+double(data(end,:)));
    SBp=data(end-6,:);
    BQp=data(end-5,:);
    MQp=data(end-4,:);
    TDFp=data(end-3,:);
    TDRp=data(end-2,:);
    
    
    
    
    
    %Find true/false of meeting thresholds
    Tminor = minorfreq > minorfreqthreshold;
    Treads= readsf > minreads_perstrand & readsr > minreads_perstrand & (f2' > minreads_perstrand_per_allele) & (r2' > minreads_perstrand_per_allele);
    % max reads per strand removed for now in favor of thresholding later
    % & readsf < maxreads_perstrand & readsr < maxreads_perstrand ;
    
    
    %Changed December 2012 to require forward and reverse strand qualities
    Tbq= ((majorbqf > min_bq) & (minorbqf > min_bq) & (majorbqr > min_bq) & (minorbqr > min_bq));
    Tmq = ((majormqf > min_mq) & (minormqf > min_mq) & (majormqr > min_mq) & (minormqr > min_mq));
    % Tbq= (majorbq > min_bq) & (minorbq > min_bq);
    % Tmq = (majormq > min_mq) & (minormq > min_mq);
    Ttd = (majortdF > min_td) & (majortdF < max_td) & (majortdR < max_td) & (majortdR > min_td)...
        & (minortdF > min_td) & (minortdF < max_td) & (minortdR > min_td) & (minortdR < max_td);
    Tid = percent_indels < max_percent_indels;
    TSBp = SBp < max_sbp;
    TBQp = BQp < max_bqp;
    %TMQp = MQp < max_mqp;
    TTDp = (TDFp < max_tdp) & (TDRp < max_tdp);
    
    
    
    %Report how many positions met each requirement
    fprintf(f,'MinorAlleleFreq: %g  \n',sum(Tminor)) ;
    fprintf(f,'Cov: %g  \n',sum(Treads)) ;
    fprintf(f,'minBQ: %g  \n',sum(Tbq)) ;
    fprintf(f,'minMQ: %g  \n',sum(Tmq)) ;
    fprintf(f,'SBp: %g  \n',sum(TSBp)) ;
    fprintf(f,'BQp: %g  \n',sum(TBQp)) ;
    %fprintf(f,'MQp: %g  \n',sum(TMQp)) ;
    fprintf(f,'TDp: %g  \n',sum(TTDp)) ;
    fprintf(f,'maxIndels: %g  \n',sum(Tid)) ;
    fprintf(f,'acceptableTD: %g  \n',sum(Ttd)) ;
    
    
    %Records positions that met all requirements
    allreqs= Tminor + Treads + Tbq + Tmq + Ttd + Tid + TSBp + TBQp + TTDp; %TMQp
    %fprintf(f,'Max requirements met: %g  \n',num_reqs) ;
    %num_reqs=max(allreqs);
    
    fprintf(f,'Positions meeting all requirements: %g  \n',sum(allreqs==num_reqs)) ;
    
    p(allreqs==num_reqs)=1;
    
    
    
    %Report how many positions are removed because of Isogenic control
    if numel(MAF_control)>1
        good=div_single_sample_test_thresholds(data, params, MAF_control, coveragethresholds);
        removedbycontrol=div_single_sample_test_thresholds(data, params, ones(size(MAF_control)), coveragethresholds);
        removedbycontrol(good>0)=0;
        fprintf(f,'Positions only removed from looseparameters because of Isogenic control: %g  \n',sum(removedbycontrol));
    else
        MAF_control=maf;
    end
    
    
else    
    MAF_control=zeros(length(data),1);
end

numfields=size(data,1);

p_sample=p;
coveragethresholds_sample=coveragethresholds;



save([tempfolder '/diverse_' char(SampleName)], 'p_sample', 'coveragethresholds_sample');

