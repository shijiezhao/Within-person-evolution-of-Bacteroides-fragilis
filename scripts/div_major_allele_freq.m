function [MAF, majorNT,minorNT, minorAF] = div_major_allele_freq(cnts)

c=cnts(1:4,:,:)+cnts(5:8,:,:);

[sorted, sortedpositions] = sort(c,1);
maxcount = sorted(end,:,:);
minorcount = sorted(end-1,:,:);

MAF=double(maxcount)./sum(c,1);
minorAF=double(minorcount)./sum(c,1);

majorNT = squeeze(sortedpositions(end,:,:));
minorNT = squeeze(sortedpositions(end-1,:,:));

MAF=squeeze(MAF);
MAF(isnan(MAF))=0; %set to 0 to indicate no data

minorAF=squeeze(minorAF); %leave no datat as 'nan' because 0 here could mean a pure position


return

