%%
% This script has the following goals:
% 1. Write the sequences in the mobile element
% 2. Write the distribution of mobile elements across samples
% 3. Combine pieces of mobile elements that are close to each other

%% Input the start/end information
mkdir 'CoverageMap';                                % Make a directory to save the outputs
fr = fastaread([REFGENOMEFOLDER '/genome.fasta']);  % fr: the genome sequence information
meancoverage=mean(all_coverage_per_bp,2);           % extract the mean coverage of each sample/genome/isolate
f=fopen('startend.txt','r');                        % open the file that contains the starts and ends information
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
mkdir('MobileElementsInformation')
fid1=fopen('MobileElementsInformation/MobileSequences.fasta','w');
fid2=fopen('MobileElementsInformation/MobileElementPattern.csv','w');
fprintf(fid2,['MobileElement'])
% 1.1 Write the first row of the pattern file
fprintf(fid2,[',' 'Type' ',' 'Coverage' ',' 'Region Length' ',' 'GC content' ',' 'Gene Annotations' ])
KK=find(~in_outgroup);
for j = 1:size(KK);
    fprintf(fid2,[',' AllSampNames{j}])
end
fprintf(fid2,['\n'])

% 2. Combine end start if they are close
newstartend=[];                          % Stitch starts/ends if they are close
for i = 1:size(startend,1);
    se = startend(i,:);
    if i==1;
        newstartend=[newstartend;se];
    else
        x1=p2chrpos(newstartend(size(newstartend,1),2),ChrStarts);x2=p2chrpos(se(1),ChrStarts);x3=p2chrpos(se(2),ChrStarts)
        if strcmp(Donor, 'D44') & newstartend(size(newstartend,1),2) + 20 > se(1) && x1(1)==x2(1);
            newstartend(size(newstartend,1),2) = se(2);
        else if newstartend(size(newstartend,1),2) + 500 > se(1) && x1(1)==x3(1);
            newstartend(size(newstartend,1),2) = se(2);
            else if x2(1)==x3(1);
                newstartend=[newstartend;se];
                end
            end
        end
    end
end
