% TVP_VAR_DPS_DMA_h1.m - Forecasting with Large TVP-VAR using forgetting factors
% MULTIPLE MODEL CASE / DYNAMIC PRIOR SELECTION (DPS) AND DYNAMIC MODEL
% AVERAGING (DMA)
%-------------------------------------------------------------------------------
    % The model is:
    %
%	 y[t] = theta[t] x[t] + e[t]
%	 theta[t] = theta[t-1] + u[t]
%
% where x[t] = I x (y[t-1],...,y[t-p]) (Kronecker product), and e[t]~N(0,V[t])
% and u[t]~N(0,Q[t]).
%
% Additionally:
    %
%  V[t] = kappa V[t-1] + (1-kappa) e[t-1]e[t-1]'
%  Q[t] = (1 - 1/lambda) S[t-1|t-1]
%
% This code estimates lambda and allows it to be time-varying. The specification is:
%
%  lambda[t] = lambda[min] + (1-lambda[min]) LL^(e[t]e[t]')
%
%-------------------------------------------------------------------------------
    %  - This code allows to calculate ONLY iterated forecasts
%  - This new version does analytical forecasting only for h=1
%  - This code "online" forecasting, i.e. the Minnesota prior should not be
%  dependent on the data, so that the Kalman filter runs once for 1:T.
%-------------------------------------------------------------------------------
    % Written by Dimitris Korobilis and Gary Koop
% University of Glasgow and University of Strathclyde
% This version: 17 January, 2012
%-------------------------------------------------------------------------------
    
    clear all;
clc;

% Add path of data and functions
addpath('data');
addpath('functions');

%-------------------------------PRELIMINARIES--------------------------------------
    forgetting = 1;    % 1: use constant factor; 2: use variable factor

lambda = 0.99;        % Forgetting factor for the state equation variance
kappa = 0.96;      % Decay factor for measurement error variance

eta = 0.99;   % Forgetting factor for DPS (dynamic prior selection) and DMA

% Please choose:
    p = 4;             % p is number of lags in the VAR part
nos = 1;           % number of subsets to consider (default is 3, i.e. 3, 7, and 25 variable VARs)
% if nos=1 you might want a single model. Which one is this?
single = 3;        % 1: 3 variable VAR
% 2: 7 variable VAR
% 3: 25 variable VAR

prior = 1;         % 1: Use Koop-type Minnesota prior
% 2: Use Litterman-type Minnesota prior

% Forecasting
first_sample_ends = 1974.75; % The end date of the first sample in the 
% recursive exercise (default value: 1969:Q4)
% NO CHOICE OF FORECAST HORIZON, h=1 IN ALL INSTANCES OF THIS CODE

% Choose which results to print
% NOTE: CHOOSE ONLY 0/1 (FOR NO/YES) VALUES!
    print_fore = 1;           % summary of forecasting results
print_coefficients = 0;   % plot volatilities and lambda_t (but not theta_t which is huge)
print_pred = 0;           % plot predictive likelihoods over time
print_Min = 0;            % print the Minnesota prior over time
%----------------------------------LOAD DATA----------------------------------------
    load ydata.dat;
load ynames.mat;
load tcode.dat;
load vars.mat;
%load yearlab.dat;

% Create dates variable
start_date = 1959.00; %1959.Q1
end_date = 2010.25;   %2010.Q2
yearlab = (1959:0.25:2010.25)';
T_thres = find(yearlab == first_sample_ends); % find tau_0 (first sample)

% Transform data to stationarity
% Y: standard transformations (for iterated forecasts, and RHS of direct forecasts)
[Y,yearlab] = transform(ydata,tcode,yearlab);

% Select a subset of the data to be used for the VAR
if nos>3
error('DMA over too many models, memory concerns...')
end

Y1=cell(nos,1);
Ytemp = standardize1(Y,T_thres);
M = zeros(nos,1);
for ss = 1:nos
if nos ~= 1
single = ss;
end
select_subset = vars{single,1};
Y1{ss,1} = Ytemp(:,select_subset);
M(ss,1) = max(size(select_subset)); % M is the dimensionality of Y
end
t = size(Y1{1,1},1);

% The first nfocus variables are the variables of interest for forecasting
nfocus = 3;

% ===================================| VAR EQUATION |==============================
% Generate lagged Y matrix. This will be part of the X matrix
x_t = cell(nos,1);
x_f = cell(nos,1);
y_t = cell(nos,1);
K = zeros(nos,1);
for ss=1:nos
ylag = mlag2(Y1{ss,1},p); 
ylag = ylag(p+1:end,:);
[temp,kk] = create_RHS(ylag,M(ss),p,t);
x_t{ss,1} = temp;
K(ss,1) = kk;
x_f{ss,1} = ylag;
y_t{ss,1} = Y1{ss,1}(p+1:end,:);
end
yearlab = yearlab(p+1:end);
% Time series observations
t=size(y_t{1,1},1); %#ok<*NASGU>

%----------------------------PRELIMINARIES---------------------------------
%========= PRIORS:
% Set the alpha_bar and the set of gamma values
alpha_bar = 10;
gamma = [1e-10,1e-5,0.001,0.005,0.01,0.05,0.1];
nom = max(size(gamma));  % This variable defines the number of DPS models
%-------- Now set prior means and variances (_prmean / _prvar)
theta_0_prmean = cell(nos,1);
theta_0_prvar = cell(nos,1);
for ss=1:nos
if prior == 1            % 1) "No dependence" prior
for i=1:nom
[prior_mean,prior_var] = Minn_prior_KOOP(alpha_bar,gamma(i),M(ss),p,K(ss));   
theta_0_prmean{ss,1}(:,i) = prior_mean;
theta_0_prvar{ss,1}(:,:,i) = prior_var;        
end
Sigma_0{ss,1} = cov(y_t{ss,1}(1:T_thres,:)); %#ok<*SAGROW> % Initialize the measurement covariance matrix (Important!)
elseif prior == 2        % 2) Full Minnesota prior
for i=1:nom
[prior_mean,prior_var,sigma_var] = Minn_prior_LITT(y_t{ss,1}(1:T_thres,:),x_f{ss,1}(1:T_thres,:),alpha_bar,gamma(i),M(ss),p,K(ss),T_thres);   
theta_0_prmean{ss,1}(:,i) = prior_mean;
theta_0_prvar{ss,1}(:,:,i) = prior_var;       
end
%Sigma_0{ss,1} = sigma_var; % Initialize the measurement covariance matrix (Important!)
Sigma_0{ss,1} = cov(y_t{ss,1}(1:T_thres,:));
end
end

% Define forgetting factor lambda:
lambda_t = cell(nos,1);
for ss=1:nos
if forgetting == 1
% CASE 1: Choose the forgetting factor   
inv_lambda = 1./lambda;
lambda_t{ss,1}= lambda*ones(t,nom);
elseif forgetting == 2
% CASE 2: Use a variable (estimated) forgetting factor
lambda_min = 0.97;
inv_lambda = 1./0.99;
alpha = 1;
LL = 1.1;
lambda_t{ss,1} = zeros(t,nom);
else
error('Wrong specification of forgetting procedure')
end
end

% Initialize matrices
theta_pred = cell(nos,1);   
theta_update = cell(nos,1);
R_t = cell(nos,1);
S_t = cell(nos,1);
y_t_pred = cell(nos,1);
e_t =cell(nos,1);
A_t = cell(nos,1);
V_t = cell(nos,1);
y_fore = cell(nos,1);
omega_update = cell(nos,1);
omega_predict = cell(nos,1);
ksi_update = zeros(t,nos);
ksi_predict = zeros(t,nos);
w_t = cell(nos,1);
w2_t = zeros(t,nos);
f_l = zeros(nom,1);
max_prob_DMS = zeros(t,1);
index_best = zeros(t,1);
index_DMA = zeros(t,nos);

anumber = t-T_thres+1;
nfore = 1;
y_t_DMA = zeros(nfore,nfocus,anumber);
y_t_DMS = zeros(nfore,nfocus,anumber);
LOG_PL_DMA = zeros(anumber,nfore);
MSFE_DMA = zeros(anumber,nfocus,nfore);
MAFE_DMA = zeros(anumber,nfocus,nfore);
LOG_PL_DMS = zeros(anumber,nfore);
MSFE_DMS = zeros(anumber,nfocus,nfore);
MAFE_DMS = zeros(anumber,nfocus,nfore);
logpl_DMA = zeros(anumber,nfocus,nfore);
logpl_DMS = zeros(anumber,nfocus,nfore);
offset = 1e-8;  % just a constant for numerical stability

%----------------------------- END OF PRELIMINARIES ---------------------------