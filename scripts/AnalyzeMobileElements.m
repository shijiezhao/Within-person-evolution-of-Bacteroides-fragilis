
%% 3. Create the files of sequences and patterns
TotalPattern=[]; TotalTYPE={}; SizeMED={};
load([REFGENOMEFOLDER '/cds_sorted.mat']);
% What to include in the files:
% Type
% Average coverage of good samples
% Length of the sequence
figure(2);hold on;
nestes=[];
NewPattern={};
TYPESEQ={};
for is = 1:size(newstartend,1); % For each mobile region
    % Close all the previous figures
    NewPattern{is}=[];              % The pattern of a particular sequence
    close all force;
    figure(1);hold on; 
    whichtype=0;                    %
    MC=[];                          % mean coverage for ?
    % Determine the start and end of the mobile elements
    p1=newstartend(is,1)+1;p2=newstartend(is,2)-1;
    chunk = all_coverage_per_bp(:,p1:p2);
    type=[];
    
    % Sequence:
    loc1 = p2chrpos(p1,ChrStarts); loc2 = p2chrpos(p2,ChrStarts);
    endloc = min(length(fr(loc1(1)).Sequence),loc2(2));
    SeqtoWrite = fr(loc1(1)).Sequence(loc1(2):loc2(2));         % Sequence to write in the output files
    
    
    % Write the pattern
    COVERAGE=[];
    cov_hist = []
    for j = 1:min(size(all_coverage_per_bp,1)-num_in_outgroup,157);      % For each isolate/genome
        normchunk = double(chunk(j,:))/double(meancoverage(j));
            plot(normchunk);
            sample=AllSampNames(j);
            coverage=mean(normchunk);
            NewPattern{is} = [NewPattern{is} coverage];
            cov_hist = [cov_hist;mean(normchunk)];  
            if coverage > 0.30 | (j>157 & coverage > 0.2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                type=[type 1];                                  % Assign 1 if average > 0.5
                COVERAGE=[COVERAGE;normchunk];                  % COVERAGE is saved for future plotting
                MC=[MC;coverage];                               % MC should be used for making the excel files
%                 if j<157;
%                     cov_hist = [cov_hist;mean(normchunk)];            
%                 end
            else;
                type=[type 0];
                AllSampNames(j)
            end        
    end
    
    MC=mean(MC);

    if is == 1 && MC > 0.0;
        TotalPattern=[TotalPattern;type];
        whichtype=1;
        TotalTYPE{1} = newstartend(is,:);
        SizeMED{1}=0;
    else if MC>0.0;
        for wis = 1:size(TotalPattern,1);
            if sum(abs(TotalPattern(wis,:)-type))/length(type)<0.05;    % Classification rules: (1). threshold: 0.5; (2). grouping: different within 5%
%               if max((TotalPattern(wis,:)-type).^2)<1;
                whichtype = wis;
                TotalTYPE{wis} = [TotalTYPE{wis};newstartend(is,:)];
            end
        end 
        if whichtype==0;
            TotalPattern=[TotalPattern;type];
            whichtype=size(TotalPattern,1);
            TotalTYPE{whichtype} = newstartend(is,:);
            SizeMED{whichtype}=0;
        end
        end
    end
    % Write the sequence infor into the fasta file
    
    % Find the sequence annotations;
    chrpos1 = p2chrpos(newstartend(is,1),ChrStarts);chrpos2 = p2chrpos(newstartend(is,2),ChrStarts);
    Annotation_MED={};
    mkk=0;
    for xyz = 1:length(CDS{chrpos1(1)});
        if (CDS{chrpos1(1)}(xyz).loc1>chrpos1(2) & CDS{chrpos1(1)}(xyz).loc1<chrpos2(2)) | (CDS{chrpos1(1)}(xyz).loc2>chrpos1(2) & CDS{chrpos1(1)}(xyz).loc2<chrpos2(2));
            Annotation_MED{mkk+1}=CDS{chrpos1(1)}(xyz).product;
            mkk=mkk+1;
        end
    end
    
    if MC>0.0;
        fprintf(fid2,['Region No.' num2str(is)]);
        fprintf(fid1,['>' Donor '_No.' num2str(is) '\t' 'Chromosome:' num2str(loc1(1)) ':' num2str(loc1(2)) '-' num2str(loc2(2)) '\t' 'Type' num2str(whichtype) '\t' 'Coverage:' num2str(MC) '\t' 'length of the chr:' num2str(length(fr(loc1(1)).Sequence)) '\n' SeqtoWrite '\n'])
        fprintf(fid2,[',' 'Type' num2str(whichtype) ',' num2str(MC) ',' num2str(loc2(2)-loc1(2)) ',' num2str(sum(SeqtoWrite=='C' | SeqtoWrite=='G')/length(SeqtoWrite)) ',' ])
        SizeMED{whichtype} = SizeMED{whichtype} + loc2(2)-loc1(2);
        if mkk>0;
            for xyz = 1:mkk;
                for lett = 1:length(Annotation_MED{xyz});
                    if sum(Annotation_MED{xyz}(lett)=='\n')>-1;
                        fprintf(fid2,[Annotation_MED{xyz}(lett)]);
                    end
                end
                fprintf(fid2,[';']);
            end
        end
        
        for metype = 1:length(type);
            fprintf(fid2,[',' num2str(type(metype))]);
        end
        saveas(gcf,['CoverageMap/' 'Type' num2str(whichtype) '_Region_No.' num2str(is)  '_noGC.pdf'],'pdf'); 
        fprintf(fid2,['\n']);
    end
    figure(444); hold on;
    map = rand(size(cov_hist,1),3);
    for kk = 1:size(cov_hist,1);
        histogram(cov_hist(kk,:),20,'Normalization','probability','facecolor',map(kk,:),'facealpha',.25,'edgecolor','none')
%         hold on;
    end
    saveas(gcf,['CoverageMap/' 'Type' num2str(whichtype) '_Region_No.' num2str(is)  '_Histogram.pdf'],'pdf'); 
    if whichtype>1;
        nestes=[nestes;newstartend(is,:)];
    end
    
    % Combine the sequences from a same mobile elements
    
end
tmpgoodpos=[];
for i = 1:length(goodpos);
    mark=0; pos = p(goodpos(i));
    for j = 1:size(newstartend,1);
        nst=newstartend(j,1);ned=newstartend(j,2);
        if pos<ned && pos>nst;
            mark=1;
        end
    end
    if mark==0;
        tmpgoodpos=[tmpgoodpos;goodpos(i)];
    end
end
goodpos=tmpgoodpos;