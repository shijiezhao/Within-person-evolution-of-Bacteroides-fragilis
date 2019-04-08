function  annotations = append_annotations(annotations, ancnti, calls, cnts, hasmutation, promoterdistance)


NTs='ATCG';

[~, majorNT, minorNT] = div_major_allele_freq(cnts);

Npositions=numel(annotations);
positions_to_annotate = 1:Npositions;
Nsamples=size(calls,2);


% Add in ancestor to each
for i = 1:Npositions
    annotations(i).anc=NTs(ancnti(i));
end

% Impose quality cut off
goodpositions=sum(hasmutation,2)>0;
positions_to_annotate = positions_to_annotate(goodpositions);




%Modify annotations such that it is a stand-alone record
for i = positions_to_annotate
    % Fill in annotations with AA, and whether promoter
    % or intergenic
    annotations(i).AApos=floor(((double(annotations(i).nt_pos)-1)/3+1)*10)/10;
    if floor(annotations(i).gene_num)==annotations(i).gene_num
        %intragenic
        protein='';
        for j=1:size(annotations(i).protein,1)
            protein=[protein '' annotations(i).protein(j,:)];
        end
        annotations(i).annotation=protein;
        if numel(annotations(i).AA)~=4
            annotations(i).type='U';
        end
        
    else
        %intergenic
        annotations(i).annotation=[num2str(annotations(i).distance1) ...
            '-' annotations(i).locustag1 ' , ' num2str(annotations(i).distance2) ...
            '-' annotations(i).locustag2 ];
        if (~isempty(annotations(i).distance1) && annotations(i).distance1 > -1*promoterdistance && annotations(i).distance1 < 0) || ...
                (~isempty(annotations(i).distance2) && annotations(i).distance2 > -1*promoterdistance && annotations(i).distance2 < 0)
            annotations(i).type='P';
        else
            annotations(i).type='I';
        end
    end
    
    % Fill in annotations with whether NS or Syn mutation
    annotations(i).AAs=''; % actual different AA's found across all samples
    
    % record ancestral
    if ancnti(i) > 0
        annotations(i).nts=[NTs(ancnti(i))];
        if numel(annotations(i).AA)==4
            annotations(i).AAs(end+1) = annotations(i).AA(ancnti(i));
        end
    else
        annotations(i).nts=char([]);
    end
    
    % Iterate through all samples
    for j = 1:Nsamples
        
        if hasmutation(i,j)>0
            
            %fixed mutation
            if calls(i,j)>0
                annotations(i).nts(end+1)=NTs(calls(i,j));
                if numel(annotations(i).AA)==4
                  annotations(i).AAs(end+1)=annotations(i).AA(majorNT(i,j));
                end
            elseif calls(i,j)==-1
            %if diverse (calls not a mutation), add minor and major call
                annotations(i).nts(end+1)=NTs(majorNT(i,j));
                annotations(i).nts(end+1)=NTs(minorNT(i,j));
                if numel(annotations(i).AA)==4
                    annotations(i).AAs(end+1)=annotations(i).AA(majorNT(i,j));
                    annotations(i).AAs(end+1)=annotations(i).AA(minorNT(i,j));
                end
            end
            
            
        end
    end
    
    if numel(annotations(i).AA)==4
        annotations(i).type='S';
    end
    
    % remove duplicates
    annotations(i).nts=unique(annotations(i).nts);
    annotations(i).AAs=unique(annotations(i).AAs);
    
%     if numel(unique(annotations(i).AAs)) > 2
%         disp(i)
%         warning('This pipeline as uploaded in 2016 is not suited for handling SNPs with 3 nucleotides')
%     end
    
    % record if nonsynonymous mutation
    if numel(annotations(i).AAs)>1
        annotations(i).type='N';
    end
    
    % Record all observed mutations across all isolates
    % E.g. K134Y, W47A, etc.
    if numel(annotations(i).AA) > 0
        annotations(i).muts={};
        rAA=annotations(i).AA(NTs==annotations(i).anc);
        nrAA=annotations(i).AAs(annotations(i).AAs~=rAA);
        for j=1:numel(nrAA)
            annotations(i).muts{j}=[rAA num2str(floor(annotations(i).AApos)) nrAA(j)];
        end
    end
    

end





end