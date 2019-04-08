function [ind,flag] = featurelocation(str)
%FEATURELOCATION Translates strings with feature locations into indices.
%
%  [IND,UFB] = FEATURELOCATION(STR) Translates strings with feature
%  locations into indices IND accordingly with the DDBJ/EMBL/GenBank
%  Feature Table Definition (http://www.ncbi.nlm.nih.gov/collab/FT/).
%  Pointers to other records, sites between symbols and inconclusive sites
%  are all parsed as NaNs. Unknown Feature Boundaries are parsed, but UFB
%  will be TRUE if it occurs in at least one section of the descriptor.

% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.1.8.3 $   $Date: 2009/01/30 14:40:46 $

if ~isempty(regexp(str,'^complement(','once'))
    [ind,flag] = featurelocation(str(12:end-1));
    ind = fliplr(ind);
elseif ~isempty(regexp(str,'^join(','once'))
    str = str(6:end-1);
    % find subsections of the string (comma separated)
    h = find([1 ~cumsum((str=='(') - (str==')')) & str==',' 1]);
    ind = [];
    flag = false;
    for i = 1:numel(h)-1
        [indT,flagT] = featurelocation(str(h(i):h(i+1)-2));
        ind = [ind indT];
        flag = flag || flagT;
    end
elseif ~isempty(regexp(str,'^.*(?=:)','once')) % points to another record
    ind = [NaN NaN];
    flag = false;
else
    flag = false;
    % sites between symbols cannot be represented with indices
    str = regexprep(str,'\d+\^\d+','NaN');
    % inconclusive locations neither can be represented with indices
    str = regexprep(str,'\(\d+\.\d+\)','NaN');
    % boundaries of unknown feature limits will be parsed, but flagged
    h = regexp(str,'[<>]');
    if ~isempty(h)
        str(h) = '';
        flag = true;
    end
    % is it a range or a single position ?
    if isempty(regexp(str,'\.\.','once'))
        ind = [str2double(str) str2double(str)];
    else
        ind = str2num(strrep(str,'.',' ')); %#ok is a vec
    end
end
    
% REFERENCE: http://www.ncbi.nlm.nih.gov/collab/FT/#3.2
%
% The following is a list of common location descriptors with their meanings: 
% Location                  Description   
% 
% 467                       Points to a single base in the presented sequence 
% 
% 340..565                  Points to a continuous range of bases bounded by and 
%                           including the starting and ending bases
% 
% <345..500                 Indicates that the exact lower boundary point of a 
%                           feature is unknown.  The location begins at some 
%                           base previous to the first base specified (which need 
%                           not be contained in the presented sequence) and con-
%                           tinues to and includes the ending base 
% 
% <1..888                   The feature starts before the first sequenced base 
%                           and continues to and includes base 888
% 
% (102.110)                 Indicates that the exact location is unknown but that 
%                           it is one of the bases between bases 102 and 110, in-
%                           clusive
% 
% (23.45)..600              Specifies that the starting point is one of the bases 
%                           between bases 23 and 45, inclusive, and the end point 
%                           is base 600 
% 
% (122.133)..(204.221)      The feature starts at a base between 122 and 133, 
%                           inclusive, and ends at a base between 204 and 221, 
%                           inclusive
% 
% 123^124                   Points to a site between bases 123 and 124
% 
% 145^177                   Points to a site between two adjacent bases anywhere 
%                           between bases 145 and 177 
% 
% join(12..78,134..202)     Regions 12 to 78 and 134 to 202 should be joined to 
%                           form one contiguous sequence
% 
% complement(join(2691..4571,4918..5163)
%                           Joins regions 2691 to 4571 and 4918 to 5163, then 
%                           complements the joined segments (the feature is 
%                           on the strand complementary to the presented strand)
%  
% join(complement(4918..5163),complement(2691..4571))
%                           Complements regions 4918 to 5163 and 2691 to 4571, 
%                           then joins the complemented segments (the feature is 
%                           on the strand complementary to the presented strand)
%   
% complement(34..(122.126)) Start at one of the bases complementary to those  
%                           between 122 and 126 on the presented strand and finish
%                           at the base complementary to base 34 (the feature is 
%                           on the strand complementary to the presented strand)
% 
% J00194:100..202           Points to bases 100 to 202, inclusive, in the entry 
%                           (in this database) with primary accession number 
%                           'J00194'
%  
