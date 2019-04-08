function Positions = generate_positions_single_sample(StrainDir, SampleName, maxFQ)


%% Tami Lieberman 2016

% An earlier version of this was not readable, re-wrote

%%


load for_matlab_scripts
[ChrStarts,GenomeLength,~,ScafNames]=genomestats(REF_GENOME_DIRECTORY);
path(SCRIPTSDIRECTORY,path);



include = zeros(GenomeLength,1) ;


fname_in=[StrainDir '/variant.vcf'];
fid=fopen(fname_in,'r');
line = fgetl(fid);

while ischar(line)
    
    if line(1)~='#'
        
        %parse line
        lt=textscan(line,'%s','Delimiter', '\t');
        l=lt{1};
        
        position_on_chr= sscanf(l{2},'%f', inf); %faster than str2double
        position=ChrStarts(strcmp(l{1},ScafNames)) + position_on_chr;
        
        alt=l{5};
        ref=l{4};
        
        %only consider for simple calls (not indel, not ambigious)
        if ~isempty(alt) & ~any(alt==',') & length(alt)==length(ref) & length(ref)==1
            
            %find and parse quality score
            xt=textscan(l{8},'%s','Delimiter', ';');
            x=xt{1};
            entrywithFQ=find(strncmp(x,'FQ',2));
            fq=sscanf(x{entrywithFQ}((find(x{entrywithFQ}=='=')+1):end),'%f', inf);  %faster than str2double
            
            if int16(fq) < maxFQ; %better than maxFQ??
                include(position)=1;
            end
        end
        
    end
    line = fgetl(fid);
    
end

fclose(fid);
Positions=p2chrpos(find(include),ChrStarts);
save([TEMPORARYFOLDER '/vcf_' char(SampleName)], 'Positions');

end

