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
        

% 2. Combine end start if they are close
newstartend=[];
for i = 1:size(startend,1);
    se = startend(i,:);
    if i==1;
        newstartend=[newstartend;se];
    else
        x1=p2chrpos(newstartend(size(newstartend,1),2),ChrStarts);x2=p2chrpos(se(1),ChrStarts);x3=p2chrpos(se(2),ChrStarts);
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

