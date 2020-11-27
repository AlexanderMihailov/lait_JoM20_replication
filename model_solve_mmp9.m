function [AA,BB,PROBLEM] = model_solve_mmp9(Theta,POLICY)

% model_solve_mmp9.m <=> AM190604 -> mmp9 version AM190611 znames last para rhog added - to work?!...
% Last modified: AM190605 -> mmp9 version; AM150804; AM150730: MINIMAL CHANGES in some equations (comments entered);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: policy=0 or 1; 1 is commitment and 0 is discretionary
%
% AM150729: 0 and 1 were swapped in this description above and below
%(inconsequential for running the codes), but I corrected it (4 places)
%
%        n = number of endogenous variables
%        THETA = vector of the parameters 
%
% output: under commitment (where policy=1)
%               Z(t) = AAZ(t-1) + BB v(t);
%               AA=(n+1)x(n+1), BB=(n+1)xns
%               Z(t)=[lambdat yt xt]'; v(t)=[eq, eg, eu, ers]';                
%         under discretion (where policy=0)
%               Z(t) = AAZ(t-1) + BB v(t);    
%               AA=(2n+1)x(2n+1), BB=(2n+1)x ns
%               Z(t)=[yt xt]'; v(t)=[eq, eg, eu, ers]';
%         prob=1 indicate convergence problem else 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calls:
%       quadpref.m      Produces loss function Q, W matrices

% Domestic economy
beta = Theta(1);                  % household discount rate
alpha = Theta(2);                 % degree of openness
h = Theta(3);                     % habit parameter
sigma = Theta(4);                 % inverse elasticity of substitution
phi = Theta(5);                   % labour elasticity
eta = Theta(6);                   % home and foreign good substitution
deltah = Theta(7);                % degree of home inflation indexation
deltaf = Theta(8);                % degree of import inflation indexation
thetah = Theta(9);                % calvo home inflation
thetaf = Theta(10);               % calvo import inflation

chi = Theta(38);                  % debt elastic parameter ---->Added by Antonio 22 Sep 2014

% Foreign VAR parameters
a1 = Theta(11);      % pi-pi
a2 = Theta(12);      % pi-output
a3 = Theta(13);      % pi-interest

b1 = Theta(14);      % output-pi
b2 = Theta(15);      % output-output
b3 = Theta(16);      % output-interest

c1 = Theta(17);      % interest-pi
c2 = Theta(18);      % interest-output
c3 = Theta(19);      % interest-interest

% Exogenous AR(1) processes
rhoh = Theta(20);                 % persistence for domestic mark-up shock
rhof = Theta(21);                 % persistence for import mark-up shock
rhoa = Theta(22);                 % persistence for productitivity shock
rhoq = Theta(23);                 % persistence for ER shock
rhos = Theta(24);                 % persistence for TOT shock
rhor = Theta(25);                 % persistence of MP shock
%rhog = Theta(39);                 % persistence of preference shock. % Added by Antonio July 2015
                                                                       % Commented out for mmp9 version by Alex 20190604
% Central bank preferences

mu_q = Theta(26);         % weight on ER
mu_y = Theta(27);         % weight on output
mu_pi = 1;                % numeraire weight on inflation
mu_r = Theta(28);         % weight on interest rate smoothing

mu_vec = [mu_q; mu_y; mu_pi; mu_r];

% Std deviation of shocks
sigh = Theta(29);       % epsilon H
sigf = Theta(30);       % epsilon F
siga = Theta(31);       % epsilon a
sigq = Theta(32);       % epsilon q
sigs = Theta(33);       % epsilon s
a4 = Theta(34);         % pi* shock
b4 = Theta(35);         % output* shock
c4 = Theta(36);         % interest* shock
sigr = Theta(37);       % epsilon r
% sigg = Theta(40);       % epsilon g (preference shock)---->Added by Antonio July 2015
                                                             % Commented out for mmp9 version by Alex 20190604
% List variables string vector

ynames = [  'Consumption     '; 
            'Dom. Inflation  ';
            'For. Inflation  ';
            'Real Ex. Rate   ';
            'ToT             ';
            'Output          ';
            'CPI Inflation   ';
            'r_{t} - r_{t-1} ';
            'r_{t}           ';
            '\epsilon_{a}    ';
            '\epsilon_{H}    ';
            '\epsilon_{F}    ';
            '\epsilon_{q}    ';
            '\epsilon_{s}    ';
            '\epsilon_{\pi*} ';
            '\epsilon_{y*}   ';
            '\epsilon_{r*}   ';
            '\epsilon_{r}    ';
			'CPI Inflation -1';      % added Tk, Feb 7, 2006
			'CPI Inflation -2';		 % added Tk, Feb 7, 2006
			'Annual CPI Inf  ';      % added Tk, Feb 7, 2006
            'For. Bond Hold  ' ] ;      % added Antonio, Sept 22, 2014
            %'\epsilon_{g}    '];	 % added Antonio, July 2015 % Commented out for mmp9 version by Alex 20190604
        
xnames = [  'r_{t}      ' ] ;

znames = [  '\epsilon_{a}    ';
            '\epsilon_{H}    ';
            '\epsilon_{F}    ';
            '\epsilon_{q}    ';
            '\epsilon_{s}    ';
            '\epsilon_{\pi*} ';
            '\epsilon_{y*}   ';
            '\epsilon_{r*}   ';
            '\epsilon_{r}    '] ;
          % '\epsilon_{g}    '] ; %(preference shock) added Antonio, July 2015
                                         % Commented out for mmp9 version by Alex 20190604    

% List of variables names as arguments in loss function: (q,ygap,rt-rt-1,pibar)
        
targets_y = [ynames(4,:); ynames(6,:); ynames(8,:); ynames(21,:)];

% Counting number of endogenous, control and exogenous variables
ny = size(ynames,1);
nx = size(xnames,1);
nz = size(znames,1);

% Shock covariance matrix
OMEGA = diag([siga,sigh,sigf,sigq,sigs,a4,b4,c4,sigr,zeros(1,ny-nz)],0); %added sigg for shock preference Antonio July 2015
%AM190604 deleted sigg

% Call central bank preference setup
[Q,W] = quadpref(mu_vec,targets_y,ynames,ny,nx);

% Setup equilibrium conditions
% Composite deep parameters
lambdah=(1-beta*thetah)*(1-thetah)/thetah;
lambdaf=(1-beta*thetaf)*(1-thetaf)/thetaf;

% STATE SPACE FORM
A0 = zeros(ny,ny);  % y_t
A1 = zeros(ny,ny);  % y_t-1
A2 = zeros(ny,ny);  % Et y_t+1
A3 = zeros(ny,nx);  % x_t
A4 = zeros(ny,nx);  % Et x_t+1
A5 = zeros(ny,ny);  % z_t


% Model setup HERE 
% Consumption Euler (as written in the paper)----> Added Antonio Sep, 2014
A0(1,1) = 1 + h;                % c(t)
A0(1,9) = (1-h)/sigma;          % r(t) 
% A0(1,23) = -(1-h)/sigma;        % epsilon_g(t) % Commented out for mmp9 version by Alex 20190604
A1(1,1) = h;                    % c(t-1)        
A2(1,1) = 1;                    % Et_c(t+1)
A2(1,7) = (1-h)/sigma;          % Et_pi(t+1)
% A2(1,23) = -(1-h)/sigma;        % Et_epsilon_g(t+1) % Commented out for mmp9 version by Alex 20190604

% Domestic goods inflation NKPC
A0(2,2) = 1 + beta*deltah;
% A0(2,4) = lambdah;
A0(2,5) = -lambdah*alpha;
A0(2,6) = -lambdah*phi;
A0(2,1) = -lambdah*sigma/(1-h);
A1(2,1) = -lambdah*sigma*h/(1-h);  % Corrected typo: missing h, Mar 30
% A0(2,17) = -lambdah*sigma/(1-h);
A1(2,2) = deltah;
% A1(2,17) = -lambdah*sigma*h/(1-h);
A2(2,2) = beta;
% A5(2,1) = lambdah;
%A5(2,2) = -lambdah*(1+phi);
A0(2,10) = lambdah*(1+phi);         % epsilon_a(t)
%A5(2,3) = lambdah;
A0(2,11) = -lambdah;                % epsilon_H(t)

% Imports inflation NKPC
A0(3,3) = 1 + beta*deltaf;
A0(3,4) = -lambdaf;
A0(3,5) = lambdaf*(1-alpha); % previously, incorrect -ve. changed typo: Mar 06
A1(3,3) = deltaf;
A2(3,3) = beta;
%A5(3,4) = lambdaf;
A0(3,12) = -lambdaf;            % epsilon_F(t)

% RUIP condition (+ chi d_t)  --> Modified by Antonio 22-09-2014;
% AM150730: a_t renamed to d_t in the comment above (inconsequential), as in the model
A0(4,4) = 1;
A0(4,9) = 1;                   
A0(4,17) = -1;                  % A0(4,18)=1 Changed by Antonio, September 2014. Possible initial mistake as it refered to epsilon_r and not epsilon_r*
A2(4,4) = 1;
A2(4,7) = 1;                   
A2(4,15) = -1;
%A5(4,5) = 1;
A0(4,13) = chi;                 % epsilon_q(t)
A0(4,22) = chi;                 % Added by Antonio 23-09-2014

% Terms of trade growth
A0(5,2) = 1;
A0(5,3) = -1;
A0(5,5) = 1;
A1(5,5)  = 1;
%A5(5,6) = 1;   % Deleted TOT shock from model: mar 30, 2006
A0(5,14) = -1;

% Market clearing using LOP gap definition
A0(6,6) = 1;
A0(6,1) = -(1-alpha);
A0(6,4) = -alpha*eta;
A0(6,5) = -alpha*eta; %incorrect: -(1-alpha*(eta-1)); changed typo: Mar 06
A0(6,16) = -alpha;

% CPI inflation index
A0(7,7) = 1;
A0(7,2) = -(1-alpha);
A0(7,3) = -alpha;

% Dummy identities: rt-rt-1 = rt - rt-1 + epsilon_r
% CORRECTED AM150729
A0(8,8) = 1;           %----> Commented out Antonio July 2015
A0(8,9) = -1;
A1(8,9) = -1;
A0(8,18) = -1;         %----> Commented out Antonio July 2015; AM150729: commented back in, but not quite sure?

% Dummy indentity rt = rt
A0(9,9) = 1;
A3(9,1) = 1;

% AR(1) technology shock
A0(10,10) = 1;
A1(10,10) = rhoa;
A5(10,1) = 1;

% AR(1) H-good markup shock
A0(11,11) = 1;
A1(11,11) = rhoh;
A5(11,2) = 1;

% AR(1) F-good markup shock
A0(12,12) = 1;
A1(12,12) = rhof;
A5(12,3) = 1;

% AR(1) RUIP shock
A0(13,13) = 1;
A1(13,13) = rhoq;
A5(13,4) = 1;

% AR(1) TOT shock
A0(14,14) = 1;
A1(14,14) = rhos;
A5(14,5) = 1;

% Foreign VAR(1) block
A0(15:17,15:17) = eye(3);
A1(15:17,15:17) = [a1,a2,a3; b1,b2,b3; c1,c2,c3];
A5(15:17,6:8) = eye(3);

% Monetary policy shock
A0(18,18) = 1;
A1(18,18) = rhor;
A5(18,9) = 1;

% Dummy; pi(-1) = pi(-1)
A0(19,19) = 1;
A1(19,7) = 1;

% Dummy; pi(-2) = pi(-2)
A0(20,20) = 1;
A1(20,19) = 1;

% Dummy; pibar(1) = [pi + pi(-1) + pi(-2) + pi(-3)]/4
A0(21,21) = 1;
A0(21,7) = -1/4;
A1(21,7) = 1/4;
A1(21,19) = 1/4;
A1(21,20) = 1/4;

% New Equations added by Antonio

% FLOW BUDGET CONSTRAINT; c + a + alpha^2 * s +alpha*q - y = 1/beta*a(-1)
% AM150729: RE-WRITTEN with d in place of a (inconsequential): d = 1/beta*d(-1) - c - alpha^2 * s - alpha*q + y
A0(22,1) = 1; %AM190604 ALL 6 ROWS HERE COMMENTED OUT in this mmp9 version where dt is not used as observable
A0(22,22) = 1;
A0(22,4) = alpha;
A0(22,5) = alpha^2;
A0(22,6) = -1;
A1(22,22) = 1/beta;

% AR(1) Preference shock
%A0(23,23) = 1; %AM190604 ALL 3 ROWS HERE COMMENTED OUT in this mmp9 version
%A1(23,23) = rhog;
%A5(23,10) = 1;

if POLICY == 0
    nlm = 0; % number of dynamic lagrange multipliers
    
    % Pass A0, A1, ...,A5, Q, W, etc to disc.m
    [F1,F2,H1,H2,PROBLEM] = discretion(A0,A1,A2,A3,A4,A5,ny,beta,Q,W);
        
    if PROBLEM == 1;
        AA=zeros(rows(H1)+rows(F1), cols(H1)+nx);
        BB=zeros(rows(H2)+rows(F2), cols(F2));
    else
        AA = [H1, zeros(ny,nx) ; 
              F1, zeros(nx,nx)];
        BB = [H2(:,1:nz); 
              F2(:,1:nz)];
    end

elseif POLICY == 1 % Commitment case
    nlm = ny; % number of dynamic lagrange multipliers
    
    [AA,BB,PROBLEM] = commitment(A0,A1,A2,A3,A4,A5,ny,beta,Q,W);
     
    if PROBLEM == 1;
        AA=zeros(2*ny+1,2*ny+1);
        BB=zeros(2*ny+1,nz);
    else
        BB=BB(:,1:nz);
    end
    
end