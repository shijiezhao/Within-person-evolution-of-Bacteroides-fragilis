%% Define the subjects that will be used here
TotSim={};      % Record the simulated mutations
TotReal={};     % Record the real mutations

%% Get the SNP information from each subject/lineage and combine them
Names={'L01','L02','L03','L04','L05','L06','L07','L08','L09','L10','L11','L12'};
workingdir=char(pwd);       % set up current directory as working dir.
cd('..')                    % move up two directories
masterdir=char(pwd);        % set up the parent dir. as masterdir
SCRIPTSDIRECTORY = [masterdir '/scripts']; path(SCRIPTSDIRECTORY,path);    % Add the srcipt folder to the search path
A=[];       % keep record of oall mutation in one variable
for name = 1:12;
    LINEAGE = Names{name};      % the lineage to look at   
    LINEAGE
    REFGENOMEFOLDER=[masterdir '/references/S' LINEAGE(2:3)];  % Reference genome directory
    cd(workingdir); cd(LINEAGE); 
    load([REFGENOMEFOLDER '/cds_sorted.mat'])
    [ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);
    % Read in the mutation information, ONE lineage after one lieage
    load(['S' LINEAGE(2:3) 'prokka.mat']);
    G=annotation_full;
    % Correct for intra-lineage nucleotide level parallel evolution
    for k = 1:numel(G);
        if G(k).pos==265509 & G(k).scaffold==5 & strcmp(LINEAGE,'L01');
            G=[G G(k)];
        end
        if G(k).pos==36564 & G(k).scaffold==17 & strcmp(LINEAGE,'L08');
            G=[G G(k)];
        end
        if G(k).pos==86601 & G(k).scaffold==14 & strcmp(LINEAGE,'L10');
            G=[G G(k) G(k)];
        end
        if G(k).pos==44770 & G(k).scaffold==12 & strcmp(LINEAGE,'L04');
            mark = k; % manually delete the double mutation, as they are likely not emerge in parallel
        end
    end
    if strcmp(LINEAGE,'L04');
        G=[G(1:mark-1) G(mark+1:length(G))];
    end

    for k = 1:numel(G);
        G(k).Subject=LINEAGE; % Assign lineage information
    end
    if ~ismember('locustag1',fieldnames(G))
        G(1).locustag1=[];
        G(1).locustag2=[];
        G(1).protein1=[];
        G(1).protein2=[];
        G(1).gene1=[];
        G(1).gene2=[];
        G(1).distance1=nan;
        G(1).distance2=nan;
    end

    A=[A G];        % keep record of all mutations in A
    annotation_full = G;
    
    % Below starts the simulation process
    % Threshold length to include genes: > 1 SNP/genelengththreshold
    genelengththreshold=2000;
    % Get the observed locus information
    locustags=get_locustags_assembly([annotation_full(:).scaffold],[annotation_full(:).pos],CDS);
    % ultg: unique locus tag; rp: the first time a unique locustag appears
    [ultg,rp] = unique(locustags(:,1));
    % look for genes with parallel mutation within subject
    multi = [];
    for i = 1:length(ultg);
        % 1. gene exists; 2. more than 1 SNP; 3. pass the length threshold
        if ultg(i)>0 & sum(ultg(i)==locustags(:,1))>1 & sum(ultg(i)==locustags(:,1))/locustags(rp(i),2) > 1/genelengththreshold;
            multi = [multi; ultg(i)];
        end
    end
    % Count for the total multiple mutation genes
    % 1. For each subject; 2. records of the locustags; 3. number of
    % mutations; 4. length of the genes.
    TotReal{name} = {{locustags},{multi},{length(multi)}};
    % Simulate for 1000 times for a particular subject
    numtrials = 1000;
    Sim_Storage={};
    for nt = 1:numtrials;
        % Randomly draw N positions in the particular genome
        % N=#mutation in the lineage
        % Draw N random positions from the genome
        sim_pos = datasample(1:GenomeLength,length(annotation_full));
        % Get the locus information
        scaf_pos = p2chrpos(sim_pos',ChrStarts);
        sim_locustags=get_locustags_assembly(scaf_pos(:,1),scaf_pos(:,2),CDS);
        % sultg: simulated unique locustag
        [sultg,srp] = unique(sim_locustags(:,1));
        smulti = [];
        % Find how many genes with in person parallel evolution
        for i = 1:length(sultg);
            if sultg(i)>0 & sum(sultg(i)==sim_locustags(:,1))>1 & sum(sultg(i)==sim_locustags(:,1))/sim_locustags(srp(i),2) > 1/genelengththreshold;
                smulti = [smulti; sultg(i)];
            end
        end
        % Store simulation results
        Sim_Storage{nt,1} = sim_locustags; 
        Sim_Storage{nt,2} = smulti;
        Sim_Storage{nt,3} = length(smulti);
    end
    TotSim{name} = {Sim_Storage};
end
annotation_full = A;
