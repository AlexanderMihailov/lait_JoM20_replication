%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per_est_qifspq_099_035_mh02_rng14twister_kll.m   % _kll appended to name, AM160518 
% Last modified: AM160822 - SEE genTheta_per099_035_kll.m in the "Functions calls" subsection below
% for the PARTICULAR calibration used (as "embedded" in the appended part of the function name:
% NEEDS TO BE CHANGED FOR ANY ALTERNATIVE CALIBRATION!)
% 
% This program estimates the original version of the KLL (2009: JMCB) model
% but with country-specific calibration of the import share in consumption alpha (see below) - 
% that is, excluding the preference shock and the assumption of incomplete asset markets,
% hence excluding also model-adjusted international reserves as the observable to proxy the foreign bond,
% as well as the parameters chi, rhog, sigg and the preference shock innovation epsilon_{g},
% for the Peruvian IMF/IFS (and OECD/National Accounts) quarterly data.
%
% It is based on Kam, Lees and Liu's (2009: JMCB) model for Australian
% data, aus_estimate.m
%
% Estimation using the Random-Walk Metropolis-Hasting Markov Chain
% Monte Carlo (RWMH-MCMC) method to simulate posterior densities
% of the parameter vector.
%
% Original author: Philip Liu  (gm_mcmc.m code)
% Last tweak: T. Kam, Feb 7, 2006.
% Final tweak by original author: PL and TK April 4, 2006.
%
% Modified by Antonio Pompa Rangel, March 1, 2014.
% Final modification and check (plus uniform formatting: max ~110 columns): Alexander Mihailov, 22 August 2016
%
% AM150731: NOW CORRECTED AND AMENDED model_likelihood_kll.m AND model_solve_kll.m as of
% 150730->160519.
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script calls:
%     
% Function calls:
%     per_datadoc_qifspq.m to manipulate and plot data series
%AM150805: added next 
%     genTheta_per099_035_kll.m - fixes the 2 calibrated parameters, beta
%          and alpha, respectively, for Peru as indicated by the filename
%          NEEDS TO BE CHANGED FOR ANY ALTERNATIVE CALIBRATION!
%     model_likelihood_kll.m for evaluating the likelihood; _kll appended to name, AM160518
%     priorln.m for evaluating the prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc
clear
tic
warning off
delete output.txt
diary('per_est_qifspq_099_035_mh02_rng14twister_kll_output.txt') 
     %AM150809: country- and calibration-specific output file name 
disp('Estimating for Peru_kll: NK-foreign VAR model & muq not zero')

LOAD_OLD = 1;

%% Policy Choice
POLICY = 0;         % 0 = discretion; 1 = commitment

% Local path options
addpath(genpath('func/'))
% addpath(genpath('stats/')) AM150724 (with MATLAB R2014a)

if LOAD_OLD == 1
    if POLICY==0
       load_path='chain/mh_dis';
    else
       load_path='chain/mh_com';
    end
end
wip_path='chain/mh_wip';

%% Settings Options
COUNTRY = 'Peru - ';
PLOT_DATA = 1;      % plot raw and filtered/dummied data
NLAGPOLICY = 1;     % = nlp if nlp #lag of r(t) appears in the state vector. Input: 0, 1, 2 .. etc
mh_scale = 0.2;     % scale variance to obtain optimal acceptance rate (20%)
N = 2*1000000;      % number of MCMC draws
MH_TRIM = 0.5;      % trimming
Nburn = MH_TRIM*N;  % number of burn in periods
loadmh = 0;         % 1 load previously saved chain, 0 start new one
ny = 21; %AM160518! % number of state variables, rows(ynames)
         %(23-2 NOW with complete markets: foreign bond to SS-ouput d_t and preference shock g_t NOT used, as in KLL)
nx = 1;             % number of policy (control/instrument) variables
PER_SAVE=1;         % periodic save option
Xsave=5000;         % saving after Xsave draws

% rand('state',34)

% AM151224: added block for seed and PRNG function for full replicability
% Generate Random Numbers That Are Repeatable - Save and Restore the Generator Settings: from Matlab documentation, see rng

% First, initialize the random number generator to make the results in this example repeatable.

rng(14,'twister'); % AM151224: here the seed and the PRNG function may be changed for EACH run

% Save the generator settings in a structure, s.

s = rng;

% Now, return the generator to the original state stored in s, and reproduce the SAME random draws (sample)

rng(s);

disp(['Number of simulations: ', num2str(N), ' draws']);
disp(['mh scale: ', num2str(mh_scale)]);
disp(['rng (14,twister)']); % AM160113 addition to print out in output file

%% Import COUNTRY data and setting data frequency and estimation startdate
[data,datadate,raw] = xlsread('data/perdata_qifspq.xlsx'); % AM150730: NEW IFS-based data set including a proxy for d(t)
% save raw_data data datadate raw
%load raw_data
%data = raw;
%data(:,3) = -data(:,3); % innvert RER series to log(eP*/P)

interval = [0.1, 0.2, 0.3, 0.4]'; % set to == 0 if annual, else, %[0.1, 0.2, 0.3, 0.4]'
sampldate = 1999; % AM150810: 2003 rather than 1999, to ignore perumns with missing data, does NOT work since the raw matrix has different dimensions
startdate = 2005;

% Call function to manipulate and plot data
[newdata,dtitle,date_vec] = per_datadoc_qifspq(COUNTRY,data,datadate,interval,sampldate,startdate,PLOT_DATA);
newdata(:,1)=[];            % get rid of the NER [AM150729: not consumption, as in KLL] series
newdata(:,10)=[];           % get rid of the consumption series c_t
                            % (was 11th column in the xlsx data, now has become 10th!!! - that's why repeated next!)
newdata(:,10)=[];           % AM160518: get rid of the model-adjusted international reserves series d_t
                            % (was 12th, now 10th) in the "_kll" complete markets version here!

% AM150805 THE NEXT BLOCK, CANNOT BE DELETED: while it seems repetitive (with regard to similar commands in
% per_datadoc_qifspq, THE REDUCTION BY 3 COLUMNS JUST ABOVE of newdata AFFECTS its SIZE, AND HENCE
% THE RE-NUMBERING of the variables by the index i!

percent_index = [2,3,4,8]; % AM150729-IMPORTANT! 10 added - was MISSING (and the other index numbers cross-checked)!
                           % AM160518: 10 deleted, as now there 9 observables
for i = 1:length(percent_index)
    newdata(:,percent_index(i)) = newdata(:,percent_index(i));
end

meandata = ones(size(newdata,1),1)*mean(newdata);
newdata = newdata-meandata;

para_latex_kll;              % Create matrix of paramater names - text/LaTeX strings
                             % -> AM160518: comes from para_latex_kll.m, a separate file that needs to be created
                             % (as well as the next command, referring to 'paraname')
npara = rows(paraname);

%% Prior Specifcations:
%--------------------------------------------------------------------------
% pshape: 0 is point mass, both para and p2 are ignored
%         1 is BETA(mean,stdd)
%         2 is GAMMA(mean,stdd)
%         3 is NORMAL(mean,stdd)
%         4 is INVGAMMA(s^2,nu)
%         5 is UNIFORM [p1,p2]
% para: Theta draws
% p1, p2 are the mean and std. of prior distribution except for uniform
%     pshape    p1      p2

pri=[   0       0.99  0.001;  % AM160528 -> added Theta vector numbering: Theta(1); beta ~ calibrated
        0       0.3495400689   0.1;   % Theta(2); alpha ~ calibrated to the country-specific mean
        1       0.6   0.2;    % Theta(3) but (1) estimated (i.e., non-zero entry in 1st column!); h ~ BETA
        2       1     0.5;    % Theta(4) but (2) estimated; sigma ~ GAMMA
        2       1.5   0.25;   % Theta(5) but (3) estimated; phi ~ NORMAL
        2       1     0.5;    % Theta(6) but (4) estimated; eta ~ GAMMA
        1       0.7   0.2;    % Theta(7) but (5) estimated; deltah ~ BETA
        1       0.7   0.2;    % Theta(8) but (6) estimated; deltaf ~ BETA
        1       0.5   0.2;    % Theta(9) but (7) estimated; thetaH ~ BETA
        1       0.5   0.2;    % Theta(10) but (8) estimated; thetaF ~ BETA
        2       0.5   0.2;    % Theta(11) but (9) estimated; a1 ~ BETA
        0       0     0.1;    % Theta(12); a2 ~ BETA
        0       0     0.1;    % Theta(13); a3 ~ NORMAL
        0       0     0.1;    % Theta(14); b1 ~ NORMAL
        2       0.5   0.2;    % Theta(15) but (10) estimated; b2 ~ BETA
        0       0     0.1;    % Theta(16); b3 ~ NORMAL
        0       0     0.1;    % Theta(17); c1 ~ NORMAL
        0       0     0.1;    % Theta(18); c2 ~ NORMAL
        2       0.5   0.2;    % Theta(19) but (11) estimated; c3 ~ BETA
        0       0     0.2;    % Theta(20); rhoh ~ BETA
        0       0     0.2;    % Theta(21); rhof ~ BETA
        1       0.5   0.2;    % Theta(22) but (12) estimated; rhoa ~ BETA
        1       0.9   0.2;    % Theta(23) but (13) estimated; rhoq ~ BETA
        1       0.25  0.2;    % Theta(24) but (14) estimated; rhos ~ BETA
        0       0     0.2;    % Theta(25); rhor ~ BETA
        2       0.5   0.3;    % Theta(26) but (15) estimated; mu_q ~ GAMMA %check genTheta_[...].m, set newTheta(26)=0 if mu_q=0
        2       0.5   0.3;    % Theta(27) but (16) estimated; mu_y ~ GAMMA
        2       0.5   0.3;    % Theta(28) but (17) estimated; mu_r ~ GAMMA
        4       0.5   0.25;   % Theta(29) but (18) estimated; sigmah ~ INVGAMMA
        4       0.5   0.25;   % Theta(30) but (19) estimated; sigmaf ~ INVGAMMA
        4       1     0.4;    % Theta(31) but (20) estimated; sigmaa ~ INVGAMMA
        4       2     0.5;    % Theta(32) but (21) estimated; sigmaq ~ INVGAMMA
        4       1     0.4;    % Theta(33) but (22) estimated; sigmas ~ INVGAMMA
        4       1     0.4;    % Theta(34) but (23) estimated; sigmapi* ~ INVGAMMA
        4       1     0.4;    % Theta(35) but (24) estimated; sigmay* ~ INVGAMMA
        4       1     0.4;    % Theta(36) but (25) estimated; sigmar* ~ INVGAMMA
        4       1     0.4;]   % Theta(37) but (26) estimated; sigmar ~ INVGAMMA
        %0       0.05    0.004;  % Theta(38); chi ~ calibrated to the value of JP(2010) in JAEtrics--->Added by Antonio 2014
                                 % AM150805 MODIFIED VALUES AND CORRECTED genTheta.m - WAS WRONG BEFORE!
        %1       0.5     0.2;    % Theta(39); rhog ~ BETA       --->Added by Antonio July 2015 (preference shock)
        %4       1       0.4];   % Theta(40); sigmag ~ INVGAMMA --->Added by Antonio July 2015

%AM160528: These last 5 parameters (BELOW) were ALREADY commented out in the KLL code - NOT USED BY US EITHER!        
%         3       2.5      1;    % c_pif
%         3       2.5     .5;    % c_pi
%         2       6       1;     % c_r
%         3       2.5     .5;    % c_pist
%         2       5       1];    % c_rst

TPOL = [15:17]; % pick index of policy parameters AM160528: numbering as the ESTIMATED parameters, NOT Theta(i) (see above)!
TPRI = [1:8];   % pick index of private deep parameters AM160528: numbering as the ESTIMATED parameters, NOT Theta(i)!
                          
% Truncations and bounds
point_b = [-100 100];  % bounds for point mass
bet_b = [0 1];         % bounds for BETA density
gam_b = [0 100];       % bounds for GAMMA density
nor_b = [-10 10];      % bounds for NORMAL density
invgam_b = [0 25];     % bounds for INVERSE GAMMA density

bounds=zeros(length(pri),2);
for i=1:length(pri)
    if pri(i,1)==0
        bounds(i,:)=point_b;
    elseif pri(i,1)==1
        bounds(i,:)=bet_b;
    elseif pri(i,1)==2
        bounds(i,:)=gam_b;
    elseif pri(i,1)==3
        bounds(i,:)=nor_b;
    else
        bounds(i,:)=invgam_b;
    end
end
pshape = pri(:,1);
p1 = pri(:,2);
p2 = pri(:,3);
lb = bounds(:,1);
ub = bounds(:,2);

%% Jumping Distribution
%%
% Load starting values else use prior mean

% AM160519 Generate an array, 'primeanTheta', containing the middle column (i.e., the means) of the prior, 'pri', matrix
% entered above and save it to a MAT-file called 'mh_init0.mat'.
%primeanTheta = pri(:,2);
%save('mh_init0.mat','primeanTheta');
% AM160519: The above line was commented out, as it did not seem to work - and was replaced by the definition of mh_init
% in the next line (taking the prior means as above) and showing them as mh_init.mat in the Workspace
%mh_init = (primeanTheta)';
load mh_init0
% AM160527 (commented out): load mh_init0
% AM160519: The SAVED KLL mh_init0.mat was used EARLIER, with 37 or 40 (as extended by Antonio) parameters - BUT these had
%           SLIGHTLY DIFFERENT VALUES! Now - see 3rd line below - code takes the means of the Theta vector (last 2 commands)
if loadmh==0,
    mh_old=0;
    Theta = mh_init ;
    Theta_s =[];
    loglike_s = [];
    logpri_s = [];
else
    load(load_path);
    nn = min(find(loglike_s==0));
    if isempty(nn)
        nn = length(loglike_s)+1;
    end
    Theta_s = Theta_s(1:nn-1,:);
    loglike_s = loglike_s(1:nn-1,:);
    logpri_s = logpri_s(1:nn-1,:);
    mh_old=length(Theta_s);
    Theta=Theta_s(mh_old,:);
    disp(['Number of previous draws loaded: ', num2str(nn-1), ' draws']);
end

%%
% Increment random vector, z ~ N(0,variance_normal) for RW-MH algorithm MCMC: 
% Theta(s) = Theta(s-1) + z
if POLICY == 0,
    variance_normal(1)=0.01;
    variance_normal(2)=0.01;
    variance_normal(3)=0.01;
    variance_normal(4)=0.01;
    variance_normal(5)=0.1;
    variance_normal(6)=0.1;
    variance_normal(7)=0.01;
    variance_normal(8)=0.01;
    variance_normal(9)=0.01;
    variance_normal(10)=0.01;
    variance_normal(11)=0.05;
    variance_normal(12)=0.01;
    variance_normal(13)=0.05;
    variance_normal(14)=0.05;
    variance_normal(15)=0.05;
    variance_normal(16)=0.05;
    variance_normal(17)=0.01;
    variance_normal(18)=0.01;
    variance_normal(19)=0.05;
    variance_normal(20)=0.01;
    variance_normal(21)=0.01;
    variance_normal(22)=0.01;
    variance_normal(23)=0.01;
    variance_normal(24)=0.01;
    variance_normal(25)=0.01;
    variance_normal(26)=0.01;
    variance_normal(27)=0.01;
    variance_normal(28)=0.01;
    variance_normal(29)=0.05;
    variance_normal(30)=0.05;
    variance_normal(31)=0.05;
    variance_normal(32)=0.05;
    variance_normal(33)=0.05;
    variance_normal(34)=0.05;
    variance_normal(35)=0.05;
    variance_normal(36)=0.05;
    variance_normal(37)=0.05;
    %variance_normal(38)=0.01;   %AM160518: commented out - and chi variance corrected from 0.05 to 0.01
    %variance_normal(39)=0.05;   %Antonio July 2015 -> AM160518: commented out - rhog
    %variance_normal(40)=0.05;   %Antonio July 2015 -> AM160518: commented out - epsilon g
%     variance_normal(41)=0.05;
%     variance_normal(42)=0.05;
else 
    variance_normal(1)=0.01;
    variance_normal(2)=0.01;
    variance_normal(3)=0.01;
    variance_normal(4)=0.01;
    variance_normal(5)=0.1;
    variance_normal(6)=0.1;
    variance_normal(7)=0.01;
    variance_normal(8)=0.01;
    variance_normal(9)=0.01;
    variance_normal(10)=0.01;
    variance_normal(11)=0.05;
    variance_normal(12)=0.01;
    variance_normal(13)=0.05;
    variance_normal(14)=0.05;
    variance_normal(15)=0.05;
    variance_normal(16)=0.05;
    variance_normal(17)=0.01;
    variance_normal(18)=0.01;
    variance_normal(19)=0.05;
    variance_normal(20)=0.01;
    variance_normal(21)=0.01;
    variance_normal(22)=0.01;
    variance_normal(23)=0.01;
    variance_normal(24)=0.01;
    variance_normal(25)=0.01;
    variance_normal(26)=0.01;
    variance_normal(27)=0.01;
    variance_normal(28)=0.01;
    variance_normal(29)=0.05;
    variance_normal(30)=0.05;
    variance_normal(31)=0.05;
    variance_normal(32)=0.05;
    variance_normal(33)=0.05;
    variance_normal(34)=0.05;
    variance_normal(35)=0.05;
    variance_normal(36)=0.05;
    variance_normal(37)=0.05;
    %variance_normal(38)=0.01;   %AM160518: commented out - chi
    %variance_normal(39)=0.05;   %Antonio July 2015 -> AM160518: commented out - rhog
    %variance_normal(40)=0.05;   %Antonio July 2015 -> AM160518: commented out - epsilon g
%     variance_normal(41)=0.05;
%     variance_normal(42)=0.05;
end
variance_normal=variance_normal*mh_scale; 

%% Metropolis-Hastings MCMC algorithm

[loglike,PROBLEM1,PROBLEM2,PROBLEM3] = ...
        model_likelihood_kll(Theta,newdata,ny,nx,POLICY,NLAGPOLICY); % _kll appended to name, AM160518

disp('Initial model Log likelihood')
disp(loglike)

logpri = priorln(Theta, pshape, p1, p2);

loglikelogpri=loglike+logpri;

save initiallikepriNK_per_est_qifspq_099_035_mh02_rng14twister_kll loglikelogpri %AM160222

pau=0;      % initial value for the acceptance rate, pau is acceptance rate 
test1=0;
probs1=0;
probs2=0;
hh = waitbar(0,'Starting RW Metropolis-Hasting');

set(hh,'Name','RWMH-MCMC: Please wait.')
Theta_s = [Theta_s; zeros(N,npara)];
loglike_s = [loglike_s; zeros(N,1)];
logpri_s = [logpri_s; zeros(N,1)];

for j=1:N,
         
        [newTheta] = genTheta_per099_035_kll(Theta, 0, variance_normal, npara, bounds);    %generate new Theta
	     
        newlogpri = priorln(newTheta, pshape, p1, p2);
        
        if newlogpri>-Inf 
                [newloglike,PROBLEM1,PROBLEM2,PROBLEM3] = ...
                        model_likelihood_kll(newTheta,newdata,ny,nx,POLICY,NLAGPOLICY);
                
                probs1=probs1+PROBLEM1;             % indeterminacy
                probs2=probs2+PROBLEM2+PROBLEM3;    % problem in evaluating the likelihood
                
                if newloglike==-Inf;
                        ratio=0;
                else
                        ratio=exp(newloglike+newlogpri-(loglike+logpri)); %alpha in notes
                        
                        if rand<=ratio              % we accept with prob ratio
                                loglike=newloglike;
                                logpri=newlogpri;
                                Theta=newTheta;
                                pau=pau+1;
                        end
                end
        end
    
	Theta_s(mh_old+j,:)=Theta;
	loglike_s(mh_old+j,1)=loglike;
	logpri_s(mh_old+j,1)=logpri;
    
    % saving per X draws
%    if rem(length(loglike_s),Xsave)==0 & PER_SAVE==1;
%        save(load_path, 'loglike_s', 'logpri_s', 'Theta_s')
%    else
%    end

    prtfrc = j/N;
    if rem(100*j,N)==0;
      if PER_SAVE==1;
        save(load_path, 'loglike_s', 'logpri_s', 'Theta_s')
        disp([num2str(mh_old+j), ' draws saved']);
      else
      end
    toc;
    disp(['completion rate (%): ', num2str(prtfrc*100)]);              % AM160204: all 3 running statistics harmonised as %
    disp(['acceptance rate (%): ', num2str(pau/j*100)]);               % AM160204: all 3 running statistics harmonised as %
    disp(['indeterminancy rate (%): ' num2str((probs1)/j*100)]);       % AM160527: "probs1+probs2" split to distinguish both
    disp(['rate of invalid likelihood (%): ' num2str((probs2)/j*100)]);% AM160527: "probs2" ONLY here, "probs1" ONLYt above
    time=fix(clock)  
    end     
    
    waitbar(prtfrc,hh,sprintf('percent done: %f; acceptance rate: %f',prtfrc*100,pau/j*100)); 
    % AM160527 (line just above) re-worded; "*100" two times inserted, to give %s ; "NK-AR-q" deleted
end  
toc;  
close(hh)

pau=pau/N*100;
probs1=probs1/N*100;
probs2=probs2/N*100; %AM160527 was wrong in the original KLL codes, now corrected, "probs2" in the RHS too, not "probs1"
ThetaMn=mean(Theta_s);

fid4=fopen('per_est_qifspq_099_035_mh02_rng14twister_kll_NKsummary.txt','w+');
fprintf(fid4,'%s ','# of draws:');
fprintf(fid4,'%d\n',N);
fprintf(fid4,'%s ','rate of acceptance (in %):');
fprintf(fid4,'%f\n',pau);
fprintf(fid4,'%s ','% of indeterminancies:');
fprintf(fid4,'%f\n',probs1);
fprintf(fid4,'%s ','% of invalid likelihood:');
fprintf(fid4,'%f\n',probs2); 
fclose(fid4);
%end

%matlabpool close
%% Marginal density
nn = min(find(loglike_s==0));
if isempty(nn)
   nn = length(loglike_s)+1;
end
Theta_s = Theta_s(1:nn-1,:);
loglike_s = loglike_s(1:nn-1,:);
logpri_s = logpri_s(1:nn-1,:);

% getting rid of fix parameters
Theta_marginal = [];
Thetapost = [];
paran = [];

parfor i=1:npara
     if pshape(i)~=0
        Theta_marginal = [Theta_marginal, Theta_s(:,i)];
        Thetapost  = [Thetapost, Theta_s(Nburn+1:end,i)];
        paran = [paran; paraname(i,:)];
     else
     end
end
marginal = marginal_density(Theta_marginal,logpri_s,loglike_s,MH_TRIM)

save(load_path, 'loglike_s', 'logpri_s', 'Theta_s', 'marginal')
save(wip_path);

%% Plot PRIOR and POSTERIOR DENSITIES of the estimated parameters

% You need to modify paran, prior, if you are calibrating some of the
% parameters. Delete the ones from the list of paraname and pri which are
% calibrated.

priordens = priorsim(pshape,p1,p2,N-Nburn);

optim = 0;  % kernel estimation bandwidth parameter

plot_density(Thetapost,priordens,optim,paran,TPOL,TPRI);

%% Plot POSTERIOR MCMC CHAINS INDIVIDUALLY for the estimated parameters

%AM160528: commented out (see next comment block)
%figure
%TPLOT = length(Thetapost);
%for i = 1:cols(Thetapost)
%    subplot(5,6,i)
%    plot(Thetapost(end-TPLOT+1:end,i))
%    title(paran(i,:),'FontSize',18)
%end
%print -depsc MCMCq
%hgsave('MCMCq')

%AM160528: Splitting the trace plot figure into two figures (for better visibility), like the two posterior figures

Thetapost01TPOL = Thetapost(:,15:17);   %AM160528: Picking up the policy parameters (numbering as estimated - see above)
Thetapost01TDEEP = Thetapost(:,1:8);    %AM160528: Picking up the deep parameters (numbering as estimated)
Thetapost01 = horzcat(Thetapost01TPOL,Thetapost01TDEEP);   %AM160528: concatenating horizontally the above 2 matrices

Thetapost02Trho = Thetapost(:,9:14);    %AM160528: Picking up the persistence parameters (numbering as estimated)
Thetapost02Tsigma = Thetapost(:,18:26); %AM160528: Picking up the st.dev. parameters (numbering as estimated)
Thetapost02 = horzcat(Thetapost02Trho,Thetapost02Tsigma);  %AM160528: concatenating horizontally the above 2 matrices

paran01TPOL = paran(15:17,:);  %AM160528: Picking up the LABELS for the policy parameters (numbering as estimated)
paran01TDEEP = paran(1:8,:);   %AM160528: Picking up the LABELS for the deep parameters (numbering as estimated)
paran01 = vertcat(paran01TPOL,paran01TDEEP);               %AM160528: concatenating verticaly the above 2 matrices

paran02Trho = paran(9:14,:);   %AM160528: Picking up the LABELS for the persistence parameters (numbering as estimated)
paran02Tsigma = paran(18:26,:);%AM160528: Picking up the LABELS for the st.dev. parameters (numbering as estimated)
paran02 = vertcat(paran02Trho,paran02Tsigma);              %AM160528: concatenating verticaly the above 2 matrices

figure %AM160528: 1st fig collecting the MCMC chanins for the policy and deep parameters
TPLOT01 = length(Thetapost01);
for i = 1:cols(Thetapost01)
    subplot(4,3,i)
    plot(Thetapost01(end-TPLOT01+1:end,i))
    title(paran01(i,:),'FontSize',16) %AM160528: Picking up the correct LABELS (paran01); fontsize reduced from 18 to 16
    end
    print -depsc 01MCMCqDeep
    hgsave('01MCMCqDeep')
    
figure %AM160528: 2nd fig collecting the MCMC chanins for the persistence and st.dev. parameters (of exogenous shocks)
TPLOT02 = length(Thetapost02);
for i = 1:cols(Thetapost02)
    subplot(4,5,i)                    %AM160528: Subplot dimesions changed from (5,4,i) to (4,5,i) (for better visibility)
    plot(Thetapost02(end-TPLOT02+1:end,i))
    title(paran02(i,:),'FontSize',16) %AM160528: Picking up the correct LABELS (paran02); fontsize reduced from 18 to 16
    end
    print -depsc 02MCMCqExog
    hgsave('02MCMCqExog')

%% Compute CONVERGENCE STATISTICS for the estimated parameters
% This section of the program computes some convergence statistics.
% All the functions used here are either from Dynare, or heavily based on a Dynare code.
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

% NSE using 4% autocovariance tapered estimate
nse4 = [];
for i=1:cols(Thetapost) % AM160519: upper bound of range linked to size of parameter vector
    nse4 = [nse4,nse(i).nse1];
end
nse4 = nse4';

% NSE using 8% autocovariance tapered estimate
nse8 = [];
for i=1:cols(Thetapost) % AM160519
    nse8 = [nse8,nse(i).nse2];
end
nse8 = nse8';

% Geweke's chi-squared test
[draws1, draws2] = splitter(Thetapost);
geweke = apm(draws1, draws2);

pvalue = [];
for i=1:cols(Thetapost) % AM160519
    pvalue = [pvalue,geweke(i).prob(2)];
end
pvalue = pvalue';

% Brooks and Gelman shrink factor

shrink = psrf(Thetapost)';

convmat = [nse8, pvalue, shrink];

%% Construct TABLES of posterior estimator statistics

priormat = [];   % AM160225: our version (of prior table)
newstatmat = []; % AM160222: our version (of statmat)

i = 1;
while i <= cols(Thetapost)
    
    prior_up = sort(priordens(:,i));
        thetapri_l = prior_up(round(0.025*(N-Nburn)));
        thetapri_u = prior_up(round(0.975*(N-Nburn)));
        thetapri_mean = mean(prior_up);
        thetapri_std = std(prior_up);
        thetapri_med = median(prior_up);
    
    post_up = sort(Thetapost(:,i));
        thetapost_l = post_up(round(0.025*(N-Nburn)));
        thetapost_u = post_up(round(0.975*(N-Nburn)));
        thetapost_mean = mean(post_up);
        thetapost_std = std(post_up);
        thetapost_med = median(post_up);

priormat(i,:) = [ thetapri_mean, thetapri_med, thetapri_std, thetapri_l, thetapri_u]; % AM160225: our version
newstatmat(i,:) = [ thetapost_mean, thetapost_std, thetapri_l, thetapri_u];           % AM160222: our version             

i = i+1;   
end

% Now write newstatmat into a TeX table

% AM160222: Generate IDENTICAL row labels for ALL table versions - we also use conv further below
for i=1:rows(paran) 
    rowlabel{:,i} = paran(i,:); 
end

collabel = {'Prior mean','median','Std','2.5%','97.5%'};
matrix2latex2(priormat, 'per_est_qifspq_099_035_kll_PERprior.tex', 'rowLabels', rowlabel, ...
            'columnLabels', collabel, 'alignment', 'c', 'format', '%-6.2f','size', 'tiny');

collabel = {'Post Mean','Post SD','2.5%','97.5%'};        
matrix2latex2(newstatmat, 'per_est_qifspq_099_035_mh02_rng14twister_kll_PERstat.tex', 'rowLabels', rowlabel, ...
            'columnLabels', collabel, 'alignment', 'c', 'format', '%-6.2f','size', 'tiny'); 
        
% AM160222: adding a table for convmat
collabel = {'NSE(8%)','p-value','B-G'};        
matrix2latex2(convmat, 'per_est_qifspq_099_035_mh02_rng14twister_kll_PERconv.tex', 'rowLabels', rowlabel, ...
            'columnLabels', collabel, 'alignment', 'c', 'format', '%-6.2f','size', 'tiny');

% AM160222: adding a final table concatenating horizontally newstatmat and convmat (as our preferred presentation)
statconvmat = horzcat(newstatmat,convmat);
collabel = {'Post Mean','Post SD','2.5%','97.5%','NSE(4%)','NSE(8%)','p-value','B-G'};        
matrix2latex2(statconvmat, 'per_est_qifspq_099_035_mh02_rng14twister_kll_PERstatconv.tex', 'rowLabels', rowlabel, ...
            'columnLabels', collabel, 'alignment', 'c', 'format', '%-6.2f','size', 'tiny');
        
%%
toc;
diary off;

disp('DONE estimating for Peru_kll: NK-foreign VAR model & muq not zero')