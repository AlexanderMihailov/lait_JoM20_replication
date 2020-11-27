function [newdata,dtitle,date_vec] = bra_datadoc_qifspq(COUNTRY,data,datadate,interval,sampldate,startdate,PLOT_DATA)

% bra_datadoc_qifspq.m
% Last modified: AM150808; AM151219 (just the comment on line 131)
%
% This script manipulates data imported into bra_estimate_qifspq.m. It also plots
% raw and "doctored" data.
% 
% T.Kam and K.Lees, Feb 2, 2006.
% AM150726: Gap variables expressed as quarterly % natural log-deviation from HP trend (i.e., *100), to be comparable to %
% variables; then ALL variables DEMEANED.
%
% Log; Feb 2, 2006: Cruel and unnatural things done to the data - HP
% filtered GDP/person, ToT, RER. Also removed GST spike in inflation data
% using least squares regression.
%
% Functions called:
%   multiplot2Draw
%   bra_multiplot2Ddoc - AM150807: country-specific time-series sample, starting for Brazil in 2000:01
%   suptitle2
%   hpfilter
%   ols - AM150804: not really used; but kept for POTENTIAL use (in the commented out block near the end
%         of the program) ONLY IF OUTLIERS ARE DUMMIED OUT (as for Australia in the original KLL code)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Option
fontsizeplus = 8; % increase default font size in function suptitle2

tq_raw = [1999.00 1999.25 1999.50 1999.75 2000.00 2000.25 2000.50 2000.75 2001.00 2001.25 2001.50 2001.75 2002.00 2002.25 2002.50 2002.75 2003.00 2003.25 2003.50 2003.75 2004.00 2004.25 2004.50 2004.75 2005.00 2005.25 2005.50 2005.75 2006.00 2006.25 2006.50 2006.75 2007.00 2007.25 2007.50 2007.75 2008.00 2008.25 2008.50 2008.75 2009.00 2009.25 2009.50 2009.75 2010.00 2010.25 2010.50 2010.75 2011.00 2011.25 2011.50 2011.75 2012.00 2012.25 2012.50 2012.75 2013.00 2013.25 2013.50 2013.75 2014.00 2014.25 2014.50 2014.75]';
      
tq_doc = [2004.00 2004.25 2004.50 2004.75 2005.00 2005.25 2005.50 2005.75 2006.00 2006.25 2006.50 2006.75 2007.00 2007.25 2007.50 2007.75 2008.00 2008.25 2008.50 2008.75 2009.00 2009.25 2009.50 2009.75 2010.00 2010.25 2010.50 2010.75 2011.00 2011.25 2011.50 2011.75 2012.00 2012.25 2012.50 2012.75 2013.00 2013.25 2013.50 2013.75 2014.00 2014.25 2014.50 2014.75]';

% 1999.00 1999.25 1999.50 1999.75 2000.00 2000.25 2000.50 2000.75 2001.00 2001.25 2001.50 2001.75 2002.00 2002.25 2002.50 2002.75 2003.00 2003.25 2003.50 2003.75

date = sampldate; % AM150806 sampldate marks the correct date labelling here
frequency = rows(interval);

date_vec = tq_raw(:); % AM150806 tq_raw indexes the correct date labelling here

i = 1;
while i < rows(datadate)-1
    date_vec = [date_vec; date*ones(rows(interval),1) + interval];
    date = date+1;
i = i + frequency;
end

date_vec = date_vec(1:end); % AM150806: was (1:end-2)
dtitle = datadate(1,2:end);

%AM150731: plot raw data - in natural log-levels or percentage (per annum, i.e, quarterly rates in annualised terms)
if PLOT_DATA
    multiplot2Draw(data,dtitle,date_vec);
    suptitle2([COUNTRY,'Raw Data ', num2str(date_vec(1)),' to ', num2str(date_vec(end))],fontsizeplus)
end

% Pass data through HP filter to get s, smoothed series for those that are
% in natural log-levels in the data input file

lambda = 1600; % smoothing parameter

[nertrend,desvabs] = hpfilter(data(:,1),lambda);
[qtrend,desvabs] = hpfilter(data(:,3),lambda);
[strend,desvabs] = hpfilter(data(:,4),lambda);
[ytrend,desvabs] = hpfilter(data(:,5),lambda);
[usytrend,desvab] = hpfilter(data(:,9),lambda);
[ctrend,desvabs] = hpfilter(data(:,11),lambda);
[dtrend,desvabs] = hpfilter(data(:,12),lambda);

% AM150729: expressing all gap variables in % terms from trend (per quarter, by *100)
nergap = (data(:,1)-nertrend)*100;
qgap = (data(:,3)-qtrend)*100;
sgap = (data(:,4)-strend)*100;
ygap = (data(:,5)-ytrend)*100;
usygap = (data(:,9)-usytrend)*100;
cgap = (data(:,11)-ctrend)*100;
dgap = ((data(:,12)-dtrend)-ytrend)*100; % AM150730: proxy for d(t) using international reserves scaled to ytrend=y_SS

newdata = data;

% AM150731: gap variables
newdata(:,1) = nergap;
newdata(:,3) = qgap;
newdata(:,4) = sgap;
newdata(:,5) = ygap;
newdata(:,9) = usygap;
newdata(:,11) = cgap;
newdata(:,12) = dgap;

% AM150731: percentage variables enter expressed in quarterly terms, and are next transformed into annualized terms
  % AM150731: inflation rates
newdata(:,2) = data(:,2)*4;
newdata(:,6) = data(:,6)*4;
newdata(:,8) = data(:,8)*4;
  % AM150731: money market (interest) rates
newdata(:,7) = data(:,7)*4;
newdata(:,10) = data(:,10)*4;

nd = (startdate-sampldate)*rows(interval)+1;
newdata = newdata(nd:end,:);

% AM150731: demeaning of the above gap and percentage variables (to proxy for zero steady state values implied by the model)
meandata = ones(size(newdata,1),1)*mean(newdata);
newdata = newdata-meandata;

% AM150731: The following block from the original KLL codes for Australia
% was commented out by Antonio, as we do not eliminate spikes in inflation
% (at least so far) - HOWEVER, the block is kept because it MAY BE USED to
% eliminate outliers in the data for a particular country

% remove GST effect on CPI inflation in 2001
%inf = newdata(:,6);
%GSTindex = find(inf==3.6566);
%inf(GSTindex) = (inf(GSTindex-1)+inf(GSTindex+1))/2;
%newdata(:,6) = %inf;
% infc=ones(rows(inf),1);
% gstdum=zeros(rows(inf),1);
% gstdum(83)=1; % we know that inf is highest on this GST date
% y=inf;
% x=[infc gstdum];
% results = ols(y,x);  % Need to install LeSage's jplv7 Econometrics Toolbox
% inf_gst = results.resid;
% newdata(:,6) = inf_gst;

date = startdate; % AM150806 startdate marks the correct date labelling here
frequency = rows(interval);

date_vec = tq_doc(:); % AM150806 bra_tq_doc indexes the correct date labelling here

i = 1;
while i < rows(datadate)-20 % AM151219 -20 = 5 years by 4 quarters, so as to start with index i=1=2004:1 is the correct date labelling here
    date_vec = [date_vec; date*ones(rows(interval),1) + interval];
    date = date+1;
i = i + frequency;
end

date_vec = date_vec(1:end)
dtitle = datadate(1,2:end);

%AM150731: plot HP-smoothed and demeaned gap data and demeaned percentage (per annum) data
if PLOT_DATA
    bra_multiplot2Ddoc(newdata,dtitle,date_vec);
    suptitle2([COUNTRY,'HP-Smoothed Data ', num2str(date_vec(1)),' to ', num2str(date_vec(end))],fontsizeplus)
end