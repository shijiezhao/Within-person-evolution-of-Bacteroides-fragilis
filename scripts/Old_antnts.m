% 6.2. Determine the major allele in the outgroup
% ancnti=mode([Calls(:,in_outgroup>0) refnti],2);     % if empty, use the reference

% Rule 1: use closest 
for i = 1:size(SampleNames);
    if strcmp(closest,SampleNames(i));
        ancnti=Calls(:,i);
        fprintf 'The closest strain is', SampleNames(i)
    end
end

% Rule 2: use
for i = 1:size(Calls,1);
    if ancnti(i) == 0;  
        tmp = Calls(i,in_outgroup>0);
        ancnti(i)=mode(tmp(tmp>0));
    end
    if ~(ancnti(i)>-1) & ( strcmp('D55',Donor) | strcmp('D52',Donor) | strcmp('D77',Donor) | strcmp('D440',Donor)); % We know which one is the ancestor, by looking at the tree
        tmp = Calls(i,in_outgroup==0);
        ancnti(i) = 0;
    end
    if strcmp('D77',Donor) & p2chrpos(p(i),ChrStarts)==[197 11111];
        ancnti(i)=1;
    end
%     if strcmp('D440',Donor);
%         if p2chrpos(p(i),ChrStarts)==[10 3120];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[13 24527];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[13 32611];ancnti(i) = 3;end
%         if p2chrpos(p(i),ChrStarts)==[79 5551];ancnti(i) = 3;end
%         if p2chrpos(p(i),ChrStarts)==[82 8870];ancnti(i) = 2;end
%         if p2chrpos(p(i),ChrStarts)==[86 11852];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[111 7245];ancnti(i) = 1;end
%         if p2chrpos(p(i),ChrStarts)==[113 23078];ancnti(i) = 3;end
%         if p2chrpos(p(i),ChrStarts)==[207 1818];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[207 5085];ancnti(i) = 3;end
%         if p2chrpos(p(i),ChrStarts)==[255 5515];ancnti(i) = 2;end % unclear
%         if p2chrpos(p(i),ChrStarts)==[259 1345];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[259 6575];ancnti(i) = 3;end
%         if p2chrpos(p(i),ChrStarts)==[268 1357];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[268 2492];ancnti(i) = 1;end
%         if p2chrpos(p(i),ChrStarts)==[277 163];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[333 9148];ancnti(i) = 3;end
%         if p2chrpos(p(i),ChrStarts)==[333 10069];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[344 5096];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[395 273];ancnti(i) = 4;end
%         if p2chrpos(p(i),ChrStarts)==[486 198];ancnti(i) = 4;end
%     end
%     if strcmp('D66',Donor);
%         if p2chrpos(p(i),ChrStarts)==[63 10969];ancnti(i) = 3;end
%         if p2chrpos(p(i),ChrStarts)==[66 1198];ancnti(i) = 2;end
%         if p2chrpos(p(i),ChrStarts)==[140 1586];ancnti(i) = 4;end
%     end
    if strcmp('D66_201707',Donor);
        if p2chrpos(p(i),ChrStarts)==[43 21573];ancnti(i) = 3;end
        if p2chrpos(p(i),ChrStarts)==[118 9227];ancnti(i) = 2;end
        if p2chrpos(p(i),ChrStarts)==[258 5266];ancnti(i) = 3;end
    end
    if strcmp('D131',Donor);
        if p2chrpos(p(i),ChrStarts)==[105 2044];ancnti(i) = 3;end
        if p2chrpos(p(i),ChrStarts)==[150 35349];ancnti(i) = 1;end
        if p2chrpos(p(i),ChrStarts)==[150 38941];ancnti(i) = 2;end
        if p2chrpos(p(i),ChrStarts)==[334 11220];ancnti(i) = 4;end % This is a tricky one... Unclear because of low mapping quality
    end
    if strcmp('D44',Donor);
        if p2chrpos(p(i),ChrStarts)==[1 1198543];ancnti(i) = 3;end
        if p2chrpos(p(i),ChrStarts)==[1 1365866];ancnti(i) = 3;end
        if p2chrpos(p(i),ChrStarts)==[4 234330];ancnti(i) = 1;end %% Wired... 
        if p2chrpos(p(i),ChrStarts)==[16 1641];ancnti(i) = 2;end
        if p2chrpos(p(i),ChrStarts)==[12 64179];ancnti(i) = 1;end
        if p2chrpos(p(i),ChrStarts)==[12 61341];ancnti(i) = 2;end
        if p2chrpos(p(i),ChrStarts)==[12 25392];ancnti(i) = 2;end
    end
end