%% F: normal direction measurements of ODs
%% R: measurements of ODs with 180 degree rotations
load('competition_exp/ODs_for_FigureS5.mat')
recrd = 28;
A=zeros(48,12);
for t = 0:5;
    for i = 1:8;
        for j = 1:12;
            ifo = t*8 + i;
            ire = t*8 + 9 - i;
            jr = 13-j;
            A(ifo,j) = 0.5*(F(ifo,j)+R(ire,jr));
        end
    end
end
%%
ODs=[];
STDs=[];
mk=1;
ODS={};
colors=rand(5,3);
h = zeros(5, 1);

for i = 1:3:10;
    for j = 1:8;
        mk
        od=[mean(A(8,4:6))];
        std=[0];
        tm=[0;3;6;9;12;15;22];
        for t = 0:5;
            od = [od; mean(A(t*8+j,i:i+2))];
            std = [std; var(A(t*8+j,i:i+2))];
        end
        ODS{mk} = od([1 3  4 5 6 7])-0.047;

        if mk > 22 & mk <27;
            figure(2);hold on;
            h(1)=errorbar(tm,od,std,'r','LineWidth',2);
        end
        if mk > 26 & mk <31;
            figure(2);hold on;
            h(2)=errorbar(tm,od,std,'b','LineWidth',2);
        end
        
        if mk > 3 & mk <16 & mod(mk,3)==1;
            figure(2);hold on;
            h(3)=errorbar(tm,od,std,'y','LineWidth',2);
            recrd=recrd+1
        end
        if mk > 3 & mk <16 & mod(mk,3)==2;
            figure(2);hold on;
            h(4)=errorbar(tm,od,std,'k','LineWidth',2);
        end
        if mk > 3 & mk <16 & mod(mk,3)==0;
            figure(2);hold on;
            h(5)=errorbar(tm,od,std,'g','LineWidth',2);
        end
        mk = mk+1;
    end
end
box on;
ax = gca;
ax.FontSize = 16; 
ax.LineWidth = 1;
xlim([0 22])

legend(h,'SL1 only','SL2 only','SL1:SL2=1:9','SL1:SL2=5:5','SL1:SL2=9:1');