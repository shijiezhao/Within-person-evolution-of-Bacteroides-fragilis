function features = parse_location_gb(s)

s(s=='<') = [] ;

if s(1)=='c'       % if 'true', gene is found on complement strand. ('features.strand = false'). 
    loc = s(12:end-1) ;     % used to skip 'complement('
    features.strand=false ;         
else
    loc = s ;
    features.strand=true ;
end


features.loc1 = str2num(loc(1:find(loc=='.',1)-1)) ;
features.loc2 = str2num(loc(find(loc=='.',1)+2:end)) ;

return
