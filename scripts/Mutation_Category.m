ct1=[];ct2=[];ct3=[];
[~,ct1]=expectedNS(annotation_11L);
[~,ct2]=expectedNS(annotation_Norm);
[~,ct3]=expectedNS(annotation_Hyper);

%%
CT=zeros(6,3);
for i =1:6;
    CT(i,1) = sum(ct1==i);
    CT(i,2) = sum(ct2==i);
    CT(i,3) = sum(ct3==i);
end
CT(:,1)=CT(:,1)/sum(CT(:,1));
CT(:,2)=CT(:,2)/sum(CT(:,2));
CT(:,3)=CT(:,3)/sum(CT(:,3));
