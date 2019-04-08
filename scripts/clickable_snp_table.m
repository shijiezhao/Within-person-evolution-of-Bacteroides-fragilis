function clickable_snp_table(annotations, callsi, cnts, SampleNames, ScafNames, MutQual, QualSort, coverageMatrix, p)

if nargin > 7
    plot_coverage_windows=1;
else
    plot_coverage_windows=0;
end
    
NTs='ATCG';

calls=repmat('N',size(callsi));
calls(callsi>0)=char(NTs(callsi(callsi>0)));


%generate table
if numel(ScafNames)>1
    colnames={'Qual', 'Type','Scaf','Pos', 'Locustag', 'Gene','Annotation', 'AApos', 'NTs', 'AAs'};
    widths=num2cell([48, 16, 12, 58, 50, 38, 250, 42, 40, 40, ones(1,numel(SampleNames))*26]);
else
    colnames={'Qual', 'Type','Pos', 'Locustag','Gene','Annotation', 'AApos', 'NTs', 'AAs'};
    widths=num2cell([48, 16, 58, 50, 38, 250, 45, 40, 40, ones(1,numel(SampleNames))*26]);
end

nonsamplecols=numel(colnames);
for i=1:numel(SampleNames)
    colnames{end+1}=SampleNames{i};
end

tabledata=cell(numel(annotations),numel(colnames));




for i=1:numel(annotations)
    if numel(annotations(i).locustag)>0
        locustag=annotations(i).locustag();
    else
        locustag=0;
    end
    
       
    % ___ generate table ___ %
    
    if numel(ScafNames)>1
        tabledata(i,1:nonsamplecols)=[{MutQual(i)} {annotations(i).type} {annotations(i).scaffold} {annotations(i).pos} ...
            {locustag} {annotations(i).gene} {annotations(i).annotation} {annotations(i).AApos} {[annotations(i).nts]} ...
            {[annotations(i).AAs]}];
    else
        tabledata(i,1:nonsamplecols)=[ {MutQual(i)} {annotations(i).type} {annotations(i).pos} ...
            {locustag} {annotations(i).gene} {annotations(i).annotation} {annotations(i).AApos} {[annotations(i).nts]} ...
            {[annotations(i).AAs]} ];
    end
    
    for j = 1:numel(SampleNames)
        tabledata(i,nonsamplecols+j) = cellstr(calls(i,j));
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


%Write to csv outfile
fid=fopen('snptable.csv','w');
fprintf(fid,'Contig,Pos,Locus,Type,Gene,Protein,');
for i=1:numel(SampleNames)
    fprintf(fid,[SampleNames{i} ',']);
end
fprintf(fid,'\n');
for i=1:numel(annotations)
    fprintf(fid,[num2str(annotations(i).scaffold) ',' num2str(annotations(i).pos) ',' ...
        annotations(i).locustag ',' ...
        annotations(i).type ',' annotations(i).gene ',' annotations(i).protein ',']);
    for j=1:numel(SampleNames)
        if callsi(i,j)>0
            fprintf(fid,[NTs(callsi(i,j)) ',']);
        else
            fprintf(fid,'N,');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);


%intialize divbarchartswindow

figure(660);clf;

    function mut_matix_clicked(~, event)
                
        rc = event.Indices ;
        
        ind=positions(rc(1));

        %Bar charts of counts
        div_bar_charts(squeeze(cnts(:,ind,:)), SampleNames)
        
        
        
        if plot_coverage_windows 
            if rc(2) > nonsamplecols
                sample=rc(2)-nonsamplecols;
                disp(sample)
            else
                sample=1;
            end
            
            if annotations(ind).pos < 500
                show_coverage_window(1, coverageMatrix(:,p(ind)-annotations(ind).pos:p(ind)+500), SampleNames, sample)
            else
                disp(ind)
                disp(p(ind)-500)
                %disp(coverageMatrix(:,p(ind)-500:p(ind)+500))
                show_coverage_window(annotations(ind).pos-500, coverageMatrix(:,p(ind)-500:p(ind)+500), SampleNames, sample)
            end
        end
        
        disp(annotations(ind))
        disp(ind)
        disp(annotations(ind).anc)

    end


end