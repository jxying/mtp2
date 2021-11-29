%% Fast Projected Newton-like Method for Precision Matrix Estimation
%% with Nonnegative Partial Correlations

close all;
clear;
clc;
%%
addpath(genpath('functions'));
addpath(genpath('data'));

%% load synthetic data: BA model with degree one consisting of 1000 nodes
data=load('BA1000.mat');
Mtrue=data.Mtrue;

%% set sample size
p=size(Mtrue,1);
n= p/2;

%% Set the underlying multivariate Gaussian distribution
Mu    = zeros(1,p);     % mean
Sigma = inv(Mtrue);     % covariance matrix

%% Generate sample covariance matrix
Y = mvnrnd(Mu,Sigma,n);
S = Y'*Y/n;             % sample covariance matrix

%% Generate weighted regularization parameters
alpha = 0.008;
lmbd  = generate_weighted_lambda(S, alpha);

%% Estimate precision matrices
opts.max_iter = 1e3;
opts.tol      = 1e-6;
out           = solver_fpn(S, lmbd, opts);
X_est         = out.X_est;  % estimated precision matrix
run_time      = out.time;   % run time
