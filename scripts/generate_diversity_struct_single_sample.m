function  generate_diversity_struct_single_sample(SampleDir, SampleName)


load('for_generate_diversity_struct_single_sample');


fprintf(1,'Creating counts 3 dimensional matrix \n') ;


load([SampleDir '/diversity.mat']);
genome_size=size(data,2);


Countsi=data(:,p');

save([TEMPORARYFOLDER '/countsatp_' SampleName], 'Countsi');


end

