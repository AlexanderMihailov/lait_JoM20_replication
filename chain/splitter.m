function [res_1,res_2] = splitter(draws)
%[res_1,res_2]
%Splits draws into two. For later implementation of Geweke's chi-squared
%test
[ndraw nvar] = size(draws);
p1 = 0.5; % AM160725 0.1 = 1st 10% of the sample replaced by 0.5 = 1st 50%
p2 = 0.5;
nobs1 = round(p1*ndraw);
nobs2 = round(p2*ndraw);

draws1 = draws(1:nobs1,:);
draws2 = trimr(draws,nobs2,0);

res_1 = momentg(draws1);
res_2 = momentg(draws2);

end

