% Model predictive control Homework assingment 1
% Jelle Cruijsen & Michiel Wind
close all;
clear all;
clc

Ts = 0.1;
N = 3;

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

%% Cost function
Q = 1*eye(n);
R = 1*eye(m);
[K,P,~] = dlqr(A,B,Q,R);                    % In MATLAB documentation, u(k)=-Kx(k)

omega = blkdiag(kron(eye(N-1),Q), P);       % Kronecker product
psi = kron(eye(N),R);

G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*phi;