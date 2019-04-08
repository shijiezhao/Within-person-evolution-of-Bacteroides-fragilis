function p = chrpos2index(pos, cstarts)

if (size(pos,1) < size(pos,2)) & (size(pos,2)>2)
    pos=pos';
  %  fprintf(1,'reversed orientation of chrpos');
end


if numel(cstarts) > 1
    p=(cstarts(pos(:,1))'+pos(:,2));
else
    p=pos(:,2);
end