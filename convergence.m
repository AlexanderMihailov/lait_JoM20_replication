function [nse4, nse8, pvalue, shrink] = convergence(Thetapost)
% AM150810: use convergence.m as a function called from the main program
% AM160223: addding the output matrix from this function to be printed into a LaTeX table in the main program

%================================================== 
% convergence.m
%
% This program computes some convergence statistics.
% All the functions used here are either from Dynare,
% or heavily based on a Dynare code.
%
% Antonio Pompa Rangel - June, 2014
%
%==================================================
% Function calls:
%     momentg -> calculates nse and rse statistics
%
%==================================================

addpath(genpath('chain/'));
load('chain/mh_wip');

nse = momentg(Thetapost);

%NSE using 4% autocovariance tapered estimate
nse4 = [];
for i=1:28
    nse4 = [nse4,nse(i).nse1];
end
nse4 = nse4';

%NSE using 8% autocovariance tapered estimate
nse8 = [];
for i=1:28
    nse8 = [nse8,nse(i).nse2];
end
nse8 = nse8';

%Geweke's chi-squared test
[draws1, draws2] = splitter(Thetapost);
geweke = apm(draws1, draws2);

pvalue = [];
for i=1:28
    pvalue = [pvalue,geweke(i).prob(2)];
end
pvalue = pvalue';

%Brooks and Gelman shrink factor

shrink = psrf(Thetapost)';

convmat = [nse4, nse8, pvalue, shrink]; % AM160223 added ";" to save the matrix rather than print it