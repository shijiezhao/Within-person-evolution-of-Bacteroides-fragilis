function num=locustag2num(tag)

%This function finds the location in locustags where it starts to be only
%numeric. Sometimes, there are letters too, so this finds the first case of
%4 numbers in a row. Edited November 2016 by Tami Lieberman.

x=char(tag);
if numel(x) > 1
    isnum=x>=48 & x<58;
    num=find(isnum(1:end-3) & isnum(2:end-2) & isnum(3:end-1) & isnum(4:end),1,'last');
else
    num=0;
end

end