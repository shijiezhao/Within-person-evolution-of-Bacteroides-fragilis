%%
% This script has the following goals:
% 1. Write the sequences in the mobile element
% 2. Write the distribution of mobile elements across samples
% 3. Combine pieces of mobile elements that are close to each other

%%
% 0. Input the start/end information
mkdir 'CoverageMap';
fr = fastaread([REFGENOMEFOLDER '/genome.fasta']);
meancoverage=mean(all_coverage_per_bp,2);

f=fopen('startend.txt','r');
formatSpec='%d';
A=fscanf(f,formatSpec);
startend=[];
for i = 1:size(A);
    if mod(i,2)==1;
        B=[A(i)];
    else
        B = [B A(i)];
        startend=[startend;B];
    end
end
        
% 1. Create files to write the mobile element sequences, and the pattern of
% mobile element distribution
fid1=fopen([Donor 'MobileSequences.fasta'],'w');
fid2=fopen([Donor 'MobileElementPattern.csv'],'w');
fprintf(fid2,['MobileElement'])
% 1.1 Write the first row of the pattern file
fprintf(fid2,[',' 'Mobile Element Linkage Group' ',' 'Coverage' ',' 'Region Length' ])
for j = size(all_coverage_per_bp,1)-11:size(all_coverage_per_bp,1);
    fprintf(fid2,[',' AllSampNames{j}])
end
fprintf(fid2,['\n'])

%
% 2. Combine end start if they are close
newstartend=[];
for i = 1:size(startend,1);
    se = startend(i,:);
    if i==1;
        newstartend=[newstartend;se];
    else
        x1=p2chrpos(newstartend(size(newstartend,1),2),ChrStarts);x2=p2chrpos(se(1),ChrStarts);x3=p2chrpos(se(2),ChrStarts)
        if strcmp(Donor, 'D44') & newstartend(size(newstartend,1),2) + 20 > se(1) && x1(1)==x2(1);
            newstartend(size(newstartend,1),2) = se(2);
        else if newstartend(size(newstartend,1),2) + 500 > se(1) && x1(1)==x2(1);
            newstartend(size(newstartend,1),2) = se(2);
        else
%             if x2(1)==x3(1) & se(2)-se(1)>1000;
                
%             if loc1(1) < loc2(1);
%                 se11=se(1);se12=se(1)+length(fr(loc1(1)).Sequence)-loc1(2)-1;se21=se12+5;se22=se(2);
%                 newstartend=[newstartend;[se11,se12;se21,se22]];
%             else
                newstartend=[newstartend;se];
%             end
%             end
            end
        end
    end
end

%% 3. Create the files of sequences and patterns
TotalPattern=[];

% What to include in the files:
% Type
% Average coverage of good samples
% Length of the sequence

for is = 1:size(newstartend,1); % For each mobile region
    % Close all the previous figures
    close all force;
    figure(1);hold on; 
    whichtype=0;
    MC=[];
    % Determine the start and end of the mobile elements
    p1=newstartend(is,1)+1;p2=newstartend(is,2)-1;
    chunk = all_coverage_per_bp(:,p1:p2);
    type=[];
    
    % Sequence:
    loc1 = p2chrpos(p1,ChrStarts); loc2 = p2chrpos(p2,ChrStarts);
    endloc = min(length(fr(loc1(1)).Sequence),loc2(2))
    SeqtoWrite = fr(loc1(1)).Sequence(loc1(2):loc2(2));
    
    
    % Write the pattern
%     normchunk=double(chunk)./repmat(meancoverage,1,(size(chunk,2))); 
    
    for j = 1:size(chunk,1)-11;
        normchunk = double(chunk(j,:))/double(meancoverage(j));
%         if goodsamples(j)==1;
%             plot(normchunk);
            sample=AllSampNames(j);
            coverage = mean(normchunk);
            if coverage > 0.40; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                type=[type 1];
%                 fprintf(fid2,[',' num2str(1)]);
                MC=[MC;coverage];
            else;
                type=[type 0];
%                 fprintf(fid2,[',' num2str(0)]);
                AllSampNames(j)
            end
%         end
    end
    
    
    
    MC=mean(MC);

    if is ==1 && MC > 0.0;
        TotalPattern=[TotalPattern;type];
        whichtype=1;
    else if MC>0.0;
        for wis = 1:size(TotalPattern,1);
            if sum(abs(TotalPattern(wis,:)-type))/length(type)<0.05;
                whichtype = wis;
            end
        end 
        if whichtype==0;
            TotalPattern=[TotalPattern;type];
            whichtype=size(TotalPattern,1);
        end
        end
    end
    % Write the sequence infor into the fasta file
    type=[];
    for j = size(chunk,1)-11:size(chunk,1);
        normchunk = double(chunk(j,:))/double(meancoverage(j));
%         if goodsamples(j)==1;
            plot(normchunk);
            sample=AllSampNames(j);
            coverage = mean(normchunk);
            if coverage > 0.40; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                type=[type coverage];
%                 fprintf(fid2,[',' num2str(1)]);
            else;
                type=[type coverage];
%                 fprintf(fid2,[',' num2str(0)]);
                AllSampNames(j)
            end
%         end
    end
    
    if MC>0.0;
        fprintf(fid2,['Region No.' num2str(is)]);
        fprintf(fid1,['>Region No.' num2str(is) '\t' 'Chromosome:' num2str(loc1(1)) ':' num2str(loc1(2)) '-' num2str(loc2(2)) '\t' 'Linkage Group ' num2str(whichtype) '\t' 'Coverage:' num2str(MC) '\t' 'length of the chr:' num2str(length(fr(loc1(1)).Sequence)) '\n' SeqtoWrite '\n'])
        fprintf(fid2,[',' 'L.G.' num2str(whichtype) ',' num2str(MC) ',' num2str(loc2(2)-loc1(2))])
        for metype = 1:length(type);
            fprintf(fid2,[',' num2str(type(metype))]);
        end
        saveas(gcf,['CoverageMap/' 'L.G.' num2str(whichtype) '_Region_No.' num2str(is)  '.pdf'],'pdf'); 
        fprintf(fid2,['\n']);
    end
%      pause
    
end