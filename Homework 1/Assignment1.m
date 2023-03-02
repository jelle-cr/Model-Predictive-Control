% Model predictive control Homework assingment 1
% Jelle Cruijsen & Michiel Wind
close all;
clear all;
clc

Ts = 0.1;
N = 20;

%% Question 4.1 (discrete-time state-space prediction model)
Ac = [-0.003 0.039 0 -0.322; -0.065 -0.319 7.74 0; 0.02 -0.101 -0.429 0; 0 0 1 0];
Bc = [0.01 1; -0.18 -0.04; -1.16 0.598; 0 0];
Cc = [1 0 0 0; 0 -1 0 7.74];
Dc = 0;

ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[n,m] = size(B);

[phi, gamma] = predictionModel(A,B,N,n,m);

%% Cost function (Likely incorrect due to the cost function being different from the instructions)
Q = 10*eye(n);
R = 0.1*eye(m);
[~,P,~] = dlqr(A,B,Q,R);                    % P is the solution of the ARE in discrete-time

omega = blkdiag(kron(eye(N-1),Q), P);       % Kronecker product
psi = kron(eye(N),R);

G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*phi;

%% Question 4.2 (design a model predictive controller)
h0 = 0;
v0 = 10;
x0 = [v0; -2.2541; -2.23e-17; -0.2912];
u0 = [0.3792; 0.0203];
