function pos = p2chrpos(p, cstarts)

chr=ones(size(p));
if numel(cstarts) > 1
    for i=2:numel(cstarts)
        chr=chr+(p>cstarts(i));
    end
    positions= p - cstarts(chr)'; %Position on Chromosome
    pos=[chr positions];
else
    pos=[chr p];
end