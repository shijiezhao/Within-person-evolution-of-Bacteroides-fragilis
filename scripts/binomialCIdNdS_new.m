
function [u, l, dnds]= binomialCIdNdS_new(N, S, expected) 
    
%Takes number of N observed and number of S observed and an expected N/S
%ratio and gives the 95% confidence interval for dN/dS
%Edited to take vectors 11/12/2013


u=[];l=[];dnds=[];

for i = 1:length(N);
    i
    II=binornd(N(i)+S(i),N(i)/(N(i)+S(i)),2000,1);
    l=[l;quantile(II,0.025)/(N(i)+S(i)-quantile(II,0.025))/expected(i)];
    u=[u;min(quantile(II,0.975)/(N(i)+S(i)-quantile(II,0.975))/expected(i),128)];
    dnds=[dnds;N(i)/S(i)/expected(i)];
end

end
