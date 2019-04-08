%% Create intermidiary files for SNP and dMRCA summarization analysis
%dmrca_overtime,stddmrca_overtime: dmrca over time, with mean and std
%dsnp: identified dSNP over time, for individual dots
%t2r: tip-2-root for individual dots
%ty: major time points
%dMRCAs,SNPs: for collector curves

%% Get dMRCA and identified SNP overtime
dmrca_overtime=[];  stddmrca_overtime=[];
t2r_overtime=[];Tip2Root={};
%% Bin time points if they are within 5 days;
t=[];newTPS={};
for i = 1:length(st);
    if i == 1;
        t=[t;st(i)];
        newTPS{1}=TPS{1};
        mk = 1;
    else;
        if st(i)-t(length(t))<6;
            newTPS{mk} = [newTPS{mk};TPS{i}];
        else
            t = [t;st(i)];
            mk = mk+1;
            newTPS{mk} = TPS{i};
        end
    end
end
oldTPS=TPS; TPS=newTPS;

%% Number of major time points
snum = length(t);

%% Create the simulation results for collector curve
dMRCA_collector={};
SNP_collector={};
for ltp = 1:snum;
    % Initialization
    dMRCAs={};  antnts=ancnti(goodpos); stddMRCAs={};
    TempGoodsamples=TPS{ltp};   repeats=500;    SNPs={};
    % simulate use a subset of isolates
    for subsampling = 1:length(TempGoodsamples);
        dMRCAs{subsampling} = []; SNPs{subsampling} = [];stddMRCAs{subsampling}=[];
        for r = 1:repeats; % simulate repeats times
            dmrca = 0; totdmrca=[];
            antnts = ancnti(goodpos);
            Subsets = TempGoodsamples;
            SelectedSample=randsample(Subsets,subsampling); % randomly pick a subset of isolates
            SelMajorNT = majorNT(goodpos,SelectedSample);
            snpn=0;dmrca=0;
            for l = 1:length(antnts);
                if (min(SelMajorNT(l,:)) ~= max(SelMajorNT(l,:))); % Two alleles of this position in the subsets
                    dmrca = dmrca + sum(antnts(l) ~= SelMajorNT(l,:))/subsampling;
                    snpn = snpn+1;
                end
            end
            dMRCAs{subsampling}=[dMRCAs{subsampling};dmrca];
            SNPs{subsampling} = [SNPs{subsampling};snpn];
        end
    end
    
    tip2root=[];
    for s = 1:subsampling;      % When all isolates are included, calculate the dMRCA for each isolate
        tip2root = [tip2root;sum(antnts ~= SelMajorNT(:,s))];
        dmrca_within_tp=0;
        for l = 1:length(antnts);
            if (min(SelMajorNT(l,:)) ~= max(SelMajorNT(l,:))); % Two alleles of this position in the subsets
                dmrca_within_tp = dmrca_within_tp + sum(antnts(l) ~= SelMajorNT(l,s));
            end
        end
        totdmrca = [totdmrca; dmrca_within_tp];
    end
    dmrca_overtime = [dmrca_overtime; [mean(totdmrca) std(totdmrca)]];    % dMRCA overtime      
    t2r_overtime = [t2r_overtime;[mean(tip2root),std(tip2root)]];       % tip-2-root overtime
    Tip2Root{ltp} = tip2root;
    dMRCA_collector{ltp} = dMRCAs;
    SNP_collector{ltp} = SNPs;
end



%% Step 4: Calculate dSNP/isolate for each time point, other than the first one
% Initilization
ddSNP=[];     % All dSNP for all isolates
dSNP={};        % delta SNP

for ltp = 1:snum;      % Use time points > 1
    SubSetB = TPS{ltp};         % SubsetB is the set of isolates from a later time point
    SubSetA = TPS{1};          % SubsetA is what to compare with
    antnts = ancnti(goodpos);
    mdSNP=[];
    
    % Step 4.2
    % Basis, to make x->1, x>1, x is later time point
    antnts = ancnti(goodpos);
    mdSNP=[];
    dSNP{ltp}=[];
    % 1. Pick a subset of TPB; 2. Find how many new SNPs are there when
    % adding a new isolate

    BasisIso = SubSetA;
    Remains = SubSetB;

    BasisNT = majorNT(goodpos,BasisIso);    % Major NT from the isolates of the first time point
    
    for rpt2 = 1:length(Remains);   % Check each isolate from a later time point, and compare it to the first time point
        tmpdSNP = 0;
        rdmRemains = Remains(rpt2);
        newNT = majorNT(goodpos,rdmRemains);
        for pos = 1:length(goodpos);    % Go through each SNP position
            if length(unique(BasisNT(pos,:)))==1;   % This position is not polymorphic in the first time point
                if unique(BasisNT(pos,:)) ~= newNT(pos) & newNT(pos) ~=antnts(pos); % new mutation brings in by this isolate
                    tmpdSNP = tmpdSNP + 1;
                end
            end
        end
        dSNP{ltp} = [dSNP{ltp}; tmpdSNP];
    end
    ddSNP=[ddSNP mean(dSNP{ltp})]; 
end
