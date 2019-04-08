function color_tree(oldfilename, newname, names, colors)

fout=fopen(newname,'w');

fprintf(fout, '#NEXUS\n');
fprintf(fout, 'begin taxa;\n');
fprintf(fout, ['\tdimensions ntax=' num2str(numel(names)) ';\n']);
fprintf(fout, '\ttaxlabels\n');

for i=1:numel(names)
    if numel(names{i})<11
        fprintf(fout, ['\t' names{i} '[&!color=#-' colors{i} ']\n']);
    else
        fprintf(fout, ['\t' names{i}(1:10) '[&!color=#-' colors{i} ']\n']);
    end
end

fprintf(fout, ';\n');
fprintf(fout, 'end;\n\n');

fprintf(fout, 'begin trees;\n');

fprintf(fout, '\t tree tree_1 = [&R] ');

treenumber=1;

%copy old tree over
fid=fopen(oldfilename,'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fprintf(fout, tline);
end
fclose(fid);

fprintf(fout, '\nend;\n\n');

fprintf(fout, 'begin figtree;\n');
fprintf(fout, '\tset tipLabels.fontName="Helvetica";\n');
fprintf(fout, '\tset tipLabels.fontSize=11;\n');
fprintf(fout, '\tset tipLabels.fontStyle=1;\n');
fprintf(fout, '\tset layout.layoutType="POLAR";\n');
fprintf(fout, '\tset polarLayout.showRoot=false;\n');
fprintf(fout, 'end;\n');


