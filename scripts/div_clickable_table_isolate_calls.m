function [annotations, tabledata, annotations_all] = div_clickable_table_isolate_calls(muts, calls, allp, ancnti, cnts, fwindows, cwindows, hasmutation, MutQual, MutQualIsolates, RefGenome, ScafNames, SampleInfo, ChrStarts, promoterdistance, QualSort)
% for isolates!!

global BARCHARTYPE

% set optional variables
if nargin < 17
    QualSort = 0;
end

fprintf(1,'\nGenerating table...\n')

% params
scrsz = get(0,'ScreenSize');
strand=['g-','k-'];
window_size=floor(size(fwindows,1)/2);
IsGenomeLoaded = false ;
NTs='ATCG';
show_alignment=1;
d = hasmutation;

[maf, majorNT, minorNT] = div_major_allele_freq(cnts);

Npositions=numel(muts);
positions_to_iterate = 1:Npositions;
goodpositions = zeros(Npositions,1);
Nsamples=size(calls,2);


goodmaf=zeros(size(d));
goodmaf(d>0)=maf(d>0);


% Impose quality cut off
goodpositions=sum(hasmutation,2)>0;
positions_to_iterate = positions_to_iterate(goodpositions);


fprintf('\nGood positions length %i\n', length(goodpositions));
fprintf('\nLength of positions to iterate %i\n', length(positions_to_iterate));

%Modify annotations such that it is a stand-alone record
annotations=muts;
for i = positions_to_iterate
    % Fill in annotations with quality score, AA, and whether promoter
    % or intergenic
    annotations(i).qual = MutQual(i);
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
        if d(i,j)>0
            %major mutation
            if calls(i,j)~='N'
                annotations(i).nts(end+1)=calls(i,j);
                if numel(annotations(i).AA)==4
                  annotations(i).AAs(end+1)=annotations(i).AA(majorNT(i,j));
                end
            end
            
                        
            %if diverse (calls not a mutation), add minor and major call
            if calls(i,j)=='N' | calls(i,j)==NTs(ancnti(i))
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
    
    % convert to N mutation
    if numel(annotations(i).AAs)>1
        annotations(i).type='N';
    end
    
    % Record all observed mutations across all isolates
    % E.g. W47W, W47A, etc.
    if numel(annotations(i).AA) > 0
        annotations(i).muts={};
        rAA=annotations(i).AA(NTs==annotations(i).ref);
        nrAA=annotations(i).AAs(find(annotations(i).AAs~=rAA));
        for j=1:numel(nrAA)
            annotations(i).muts{j}=[rAA num2str(floor(annotations(i).AApos)) nrAA(j)];
        end
    end
    
    % This is for indels
    %         if sum(calls(i,:)=='I' & d(i,:),2) > 0
    %             annotations(i).nts= [annotations(i).nts 'I'];
    %             annotations(i).type='D';
    %         elseif sum(calls(i,:)=='D' & d(i,:),2) > 0
    %             annotations(i).nts= [annotations(i).nts 'D'];
    %             annotations(i).type='D';
    %         end
end

fprintf('\nFinished generating all annotations\n')

%     %for returning
%     allannotations=annotations;

%actual table
annotations_all=annotations;
annotations=annotations(goodpositions>0);
MutQualIsolates = MutQualIsolates(goodpositions>0,:);
allp=allp(goodpositions>0);
cnts(:,~goodpositions,:)=[];

if numel(fwindows)>1
    fwindows(:,~goodpositions,:)=[];
    cwindows(:,~goodpositions,:)=[];
end


%generate table
if numel(ScafNames)>1
    colnames={'Qual', 'Type','Chr','Pos', 'Locustag', 'Gene','Annotation', 'AApos', 'NTs', 'AAs'};
    widths=num2cell([48, 16, 12, 58, 50, 38, 250, 42, 40, 40, ones(1,numel(SampleInfo))*26]);
else
    colnames={'Qual', 'Type','Pos', 'Locustag','Gene','Annotation', 'AApos', 'NTs', 'AAs'};
    widths=num2cell([48, 16, 58, 50, 38, 250, 45, 40, 40, ones(1,numel(SampleInfo))*26]);
end

nonsamplecols=numel(colnames);
for i=1:Nsamples
    colnames{end+1}=SampleInfo(i).Sample;
end

tabledata=cell(sum(goodpositions),numel(colnames));




for i=1:numel(annotations)
    if numel(annotations(i).locustag)>0
        locustag=annotations(i).locustag(end-4:end);
    else
        locustag=0;
    end
    
    % ___ generate table ___ %
    
    if numel(ScafNames)>1
        tabledata(i,1:nonsamplecols)=[{[annotations(i).qual]} {annotations(i).type} {annotations(i).scaffold} {annotations(i).pos} ...
            {locustag} {annotations(i).gene} {annotations(i).annotation} {annotations(i).AApos} {[annotations(i).nts]} ...
            {[annotations(i).AAs]}];
    else
        tabledata(i,1:nonsamplecols)=[ {[annotations(i).qual]} {annotations(i).type} {annotations(i).pos} ...
            {locustag} {annotations(i).gene} {annotations(i).annotation} {annotations(i).AApos} {[annotations(i).nts]} ...
            {[annotations(i).AAs]} ];
    end
    
    % FOR SCRAPES
    % ___ fill in mutAF for each isolate ___ %
    %         for j=1:Nsamples
    %             if (annotations(i).mutAF(j) > 0) && (annotations(i).mutAF(j) < 1)
    %                 n=[num2str(annotations(i).mutAF(j)) '0' '0'];
    %                 tabledata(i,nonsamplecols+j)=cellstr(n(2:4));
    %             elseif (annotations(i).mutAF(j) == 1)
    %                 tabledata(i,nonsamplecols+j)={'1.0'};
    %             elseif (annotations(i).mutAF(j) == -1)
    %                 tabledata(i,nonsamplecols+j)={'I'};
    %             elseif (annotations(i).mutAF(j) == -2)
    %                 tabledata(i,nonsamplecols+j)={'D'};
    %             else
    %                 tabledata(i,nonsamplecols+j)={'0'};
    %             end
    %         end
    
    % FOR ISOLATES
    % ___ fill in maNT for each isolate ___ %
    for j = 1:Nsamples
        tabledata(i,nonsamplecols+j) = cellstr(calls(positions_to_iterate(i),j));
    end
end


positions = 1:size(tabledata,1);
oldtable=tabledata;
if QualSort==1
    [tabledata, positions] = sortrows(oldtable,1);
end

%display table
figure();clf;hold on;
set(gcf, 'Position',[10         50        1250         550]);
uicontrol('Style','text','Position',[400 45 120 20],'String','Vertical Exaggeration')
t = uitable('Units','normalized','Position',[0 0 1 .97], 'Data', tabledata,...
    'ColumnName', colnames,...
    'RowName',[], ...
    'CellSelectionCallback',@mut_matix_clicked, ...
    'ColumnWidth',widths);
h.checkbox1 = uicontrol('Units','normalized','Style','checkbox','String','Show Alignment in IGV when clicked (must have IGV viewer open already)','Min',0,'Max',1, 'Value',0, 'Position',[0 .97 1 .03]);

%intialize divbarchartswindow
figure(660);clf;
set(660,'Position',[scrsz(3)*2/3 scrsz(4)/20 scrsz(3)/3 scrsz(4)/2]);clf;hold on;


    function mut_matix_clicked(src, event)
        % ___ MODIFY TO SHOW ONLY ISOLATES IN MUTQUALISOLATES!! ___ %
        
        scrsz = get(0,'ScreenSize');
        strand=['g-','k-'];
        window_size=floor(size(fwindows,1)/2);
        
        rc = event.Indices ;
        dt = get(src,'data') ;
        
        ind=positions(rc(1));
        
        if isfield(annotations(ind),'scafold')
            chr=annotations(ind).scafold;
        else
            chr=1;
        end
        position= annotations(ind).pos;
        
        disp(ind);
        disp(annotations(ind));
        
        % get the two isolates used for MutQual in this position
        calledisolates = MutQualIsolates(ind,:);
        fprintf('\nIsolates used to make call is %i and %i\n', calledisolates(1), calledisolates(2));
        
        if rc(2) > nonsamplecols
            sample=rc(2)-nonsamplecols;
            disp(sample)
        else
            sample=1;
        end
        
        %Bar charts of counts
        if BARCHARTYPE==1
            div_bar_charts(squeeze(cnts(:,ind,:)), sample, {SampleInfo.Sample})
        elseif BARCHARTYPE==2
            div_bar_charts(squeeze(cnts(:,ind,calledisolates)), sample, {SampleInfo(calledisolates).Sample})
        elseif BARCHARTYPE==3
            div_bar_charts(squeeze(cnts(:,ind,[sample calledisolates])), sample, {SampleInfo([sample calledisolates]).Sample})
        end
        
        %Plot MAF in region neighboring locus
        region=(find(allp>allp(ind)-window_size,1):find(allp<allp(ind)+window_size,1,'last'));
        
        if numel(fwindows)>1
            div_maf_window(annotations(ind), allp(ind)-ChrStarts(chr), window_size, [], squeeze(fwindows(:,ind,:)), allp(region)-ChrStarts(chr), goodmaf(region,:), {SampleInfo.Sample},sample,0)
            div_cov_window(annotations(ind), allp(ind)-ChrStarts(chr), window_size, squeeze(cwindows(:,ind,:)), {SampleInfo.Sample},sample, 0)
        end
        
        show_alignment=get(h.checkbox1, 'Value');
        if show_alignment==1
            
            t = tcpip('localhost', 60152) ;
            fopen(t) ;
            if rc(2) > nonsamplecols
                bai_name = ['/Volumes/sysbio/KISHONY LAB/illumina_pipeline/' SampleInfo(sample).ExperimentFolder '/' SampleInfo(sample).Sample '/' SampleInfo(sample).AlignmentFolder '/aligned.sorted.bam.bai' ]
                
                
                if ~exist(bai_name,'file')
                    error('Create bai files for viewing alignment')
                end
                
                
                if ~IsGenomeLoaded
                    run_cmd(['genome  /Volumes/sysbio/kishonylab/illumina_pipeline/Reference_Genomes/' RefGenome '/genome.fasta' ])
                    IsGenomeLoaded = true ;
                end
                
                run_cmd(['load "' bai_name(1:end-4) '"']) ;
                
                run_cmd(['goto "' ScafNames{chr} ':' num2str(position) '"'])
                
            end
            
            
            fclose(t);
            delete(t);
            clear t
            
        end
        
        function run_cmd(c)
            disp(c)
            fprintf(t, c);
            response = fgetl(t);
        end
        
        
        
    end


end