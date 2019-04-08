
function [u, l, dnds]= binomialCIdNdS(N, S, expected) 
    
%Takes number of N observed and number of S observed and an expected N/S
%ratio and gives the 95% confidence interval for dN/dS
%Edited to take vectors 11/12/2013


t=N+S; %number of trials

e=N;
p=N./t;



[~,b] = binofit(N, t);

l=b(:,1).*t;
u=b(:,2).*t;

u=u./(t-u)./expected;
l=l./(t-l)./expected;

dnds=(N./S)./expected;


end
