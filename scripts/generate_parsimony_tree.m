function treefile = generate_parsimony_tree(calls, names, postfix)




timestamp=datestr(now, 'yyyy-mm-dd-HH-MM-SS');

%write input file
generate_phylip_input(calls, names, [timestamp '_' postfix '_infile.txt'])

%write parameter file
fid = fopen([timestamp '_' postfix '_optionfile.txt'],'w');
fprintf(fid, [timestamp '_' postfix '_infile.txt\n']);
fprintf(fid, 'f\n');
fprintf(fid, [timestamp '_' postfix '_out.txt\n']);
fprintf(fid, '5\n');
fprintf(fid, 'V\n');
fprintf(fid, '1\n');
fprintf(fid, 'y\n');
fprintf(fid, 'f\n');
fprintf(fid, [timestamp '_' postfix '_out.tree\n\n']);
treefile=[timestamp '_' postfix '_out.tree'];
fclose(fid);

%run
fid = fopen('temp.sh','w');
fprintf(fid, ['! ../dnapars < ' timestamp '_' postfix '_optionfile.txt > ' timestamp '_' postfix '_outfile.txt\n']);
fclose(fid);

if ~exist('outfile','file')
    fid = fopen('outfile','w');
    fprintf(fid, 'foo');
    fclose(fid);
    fid = fopen('outtree','w');
    fprintf(fid, 'foo');
    fclose(fid);
end


! chmod +x temp.sh
! ./temp.sh

