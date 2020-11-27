function bra_multiplot2Ddoc(ydata,ydataname,tq_doc) 
% Last modified: AM150810; AM151219 (just the font and comment on line 44)
% AM150719 varargin replaced by bra_tq_doc
% AM150808: country-specific time-series data sample allowed, and indicated in the file name.

% Based on: multiplot2D.m  (c) T.Kam 2006
%
% This function plots multiple data series in subplots typically suitable
% for journal layout sizes
%
% data is T x n matrix, where T is #points, n is #variables
% dataname is n x 1 string vector
% varargin is "xdata" for x-axis labels

% AM150804: bra_tq_doc needs to be modified to accommodate a different sample
% of "doctored" data, by moving/commenting in or out certain dates - as here below

tq_doc = [2004.00 2004.25 2004.50 2004.75 2005.00 2005.25 2005.50 2005.75 2006.00 2006.25 2006.50 2006.75 2007.00 2007.25 2007.50 2007.75 2008.00 2008.25 2008.50 2008.75 2009.00 2009.25 2009.50 2009.75 2010.00 2010.25 2010.50 2010.75 2011.00 2011.25 2011.50 2011.75 2012.00 2012.25 2012.50 2012.75 2013.00 2013.25 2013.50 2013.75 2014.00 2014.25 2014.50 2014.75]';
 
% 1999.00 1999.25 1999.50 1999.75 2000.00 2000.25 2000.50 2000.75 2001.00 2001.25 2001.50 2001.75 2002.00 2002.25 2002.50 2002.75 2003.00 2003.25 2003.50 2003.75  

nwin = size(ydata,2); % #subplots, nwin = n
yscale = 0.2; % create buffer with size of +/- 0 <= yscale <= 1

figure
    if nwin <= 6
        wc = 2;                           % window cols
    elseif nwin > 6 & nwin <= 9
        wc = 3;
    else
        wc = 4;
    end
    
        wr = ceil(nwin/wc); % window rows
    
    for i = 1:nwin     
        subplot(wr,wc,i)
            
            if nargin <= 2
                plot(tq_doc, ydata(:,i));
            elseif nargin > 2
                plot(tq_doc,ydata(:,i));
            end
            ylabel(ydataname(i),'FontSize',10) % AM151219: smaller now (due to my long names of series), was 16
        grid on
       
        axis([min(tq_doc), max(tq_doc), ...
              (min(ydata(:,i)) - (max(ydata(:,i))-min(ydata(:,i)))*yscale), ...
              (max(ydata(:,i)) + (max(ydata(:,i))-min(ydata(:,i)))*yscale)])
    end