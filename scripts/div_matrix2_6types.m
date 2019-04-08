function x = div_matrix2_6types(matrix)

x=zeros(6,size(matrix,3));

%compress into 6 elements

%AT, TA
%AC, TG
%AG, TC
%GC, CG
%GT, CA
%GA, CT

x(1,:)=squeeze(matrix(1,2,:)+matrix(2,1,:));
x(2,:)=squeeze(matrix(1,3,:)+matrix(2,4,:));
x(3,:)=squeeze(matrix(1,4,:)+matrix(2,3,:));
x(4,:)=squeeze(matrix(4,3,:)+matrix(3,4,:));
x(5,:)=squeeze(matrix(4,2,:)+matrix(3,1,:));
x(6,:)=squeeze(matrix(4,1,:)+matrix(3,2,:));

end