function  Counts = generate_diversity_struct_simple(SampleDirs, SampleNames, p, numfields)



Counts=zeros(numfields, length(p), length(SampleDirs),'uint16');


for i=1:length(SampleDirs)
    
    fprintf(1,'Loading counts matrix for sample: %g  \n',i) ;
    load([SampleDirs{i} '/diversity.mat']);
    Counts(:,:,i)=data(:,p);
   
end


end