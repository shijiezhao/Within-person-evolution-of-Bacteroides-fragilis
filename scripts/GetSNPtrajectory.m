function [MIC,Days,MutCounts,TotCounts] = GetSNPtrajectory(i,counts,bin,N,Date,optc,smt)
    data = [];                  % Raw data corresponding to the SNP
    date = [];                  % Only include dates with positive read counts
    MutCounts = [];
    TotCounts = [];
    for j=1:size(counts,3);     % For each time point (206 in total)
        if sum(counts(1:8,i,j))>N;          % If sample has more than N counts
            date = [date Date(j)];
            data=[data counts(1:8,i,j)];
            As=sum(data(1,:))+sum(data(5,:));
            Ts=sum(data(2,:))+sum(data(6,:));
            Cs=sum(data(3,:))+sum(data(7,:));
            Gs=sum(data(4,:))+sum(data(8,:));
        end
    end
    
    % Find the indecies of major allele and minor allele
    [ii,ii] = sort([As Ts Cs Gs]);
    x11 = ii(4);
    x22 = ii(3);
    if x11==optc(i,3);
        mai=x11;mii=x22;
    else
        mai=x22;mii=x11;
    end
    MAC=[];MIC=[];              % Major allele counts and Minor allele counts
    Aseq=data(1,:)+data(5,:);
    Tseq=data(2,:)+data(6,:);
    Cseq=data(3,:)+data(7,:);
    Gseq=data(4,:)+data(8,:);
    Seqs=[Aseq;Tseq;Cseq;Gseq];
    Days=[];
    for j=1:length(data)/bin;                   % Group samples based on bin size
        mac = sum(Seqs(mai,(j-1)*bin+1:j*bin)); % Major allele counts
        mic = sum(Seqs(mii,(j-1)*bin+1:j*bin)); % Minor allele counts
        meanday = mean(date((j-1)*bin+1:j*bin));
        Days = [Days meanday];
        MAC=[MAC mac/(mac+mic)];
        MIC=[MIC mic/(mac+mic)];
        MutCounts = [MutCounts mic];
        TotCounts = [TotCounts mac+mic];
    end
%     for xyz = 1:smt+1;
%         if xyz>1;
%          MIC = smooth(MIC);
%         end
%     end
%     MIC=smooth(MIC);
%     MIC=smooth(MIC);
    return