% Model predictive control Homework assingment 1
% Jelle Cruijsen & Michiel Wind
close all;
clear all;
clc

t = 300; 
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
p = height(C);

% Initial conditions
h0 = 0;
v0 = 10;
x0 = [v0; -2.2541; -2.23e-17; -0.2912];
u0 = [0.3792; 0.0203];
y0 = [v0; h0];

ref = zeros(2*(t+N),1);
k = 0;
for i = 1:2:length(ref)
    if(k<100)
        ref(i) = 10;
        ref(i+1) = 0;
    elseif(k<200)
        ref(i) = 10;
        ref(i+1) = 5;
    else
        ref(i) = 8;
        ref(i+1) = 5;
    end
    k = k+1;
end
ref1 = zeros((t+N),1);
k = 0;
for i = 1:2:length(ref)
    k = k+1;
    ref1(k) = ref(i);
end
ref2 = zeros((t+N),1);
k = 0;
for i = 1:2:length(ref)
    k = k+1;
    ref2(k) = ref(i+1);
end

%% Constraints and matrices
umin = [-20; -20];
umax = [20; 20];
xmin = [-25; -25; -pi/20; -pi/2];
xmax = [25; 25; pi/20; pi/2];
ymin = -15;
ymax = 15;
dumin = [-15; -15];
dumax = [15; 15];

Q = 10*eye(p);
R = 0.1*eye(m);

[phi, gamma] = predictionModel(A,B,N,n,m);
omega = kron(eye(N),Q);       % Kronecker product
psi = kron(eye(N),R);
C_bar = kron(eye(N),C);

T = diag(ones(N*m,1));
nT = size(T,1);
T(2:nT+1:end) = -1;

%% Question 4.2 (design a model predictive controller)
k_sim = t;
%% unconstrained MPC
% constraint = 'unconstrained';           %changes plot title
% uk = [u0 zeros(m,t)];
% xk = [x0 A*x0+B*u0 zeros(n,t)];
% yk = [y0 C*xk(:,2) zeros(p,t-1)];
% G = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% i = 0;
% for k = 2:(k_sim+1)               %time starts at k = 0, xk is known for k = 1
%     i = i+2;                      %used for indexing due to reference being twice as long;
%     Rk = ref((i+1):(i+p*N));
%     v = [uk(:,k-1); zeros(2*(N-1),1)];
%     Uk = -2*G^-1*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
%     uk(:,k) = Uk(1:m);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k) = C*xk(:,k);
% end
%% constrained MPC
constraint = 'constrained';           %changes plot title
%[Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m);  %no delta u and y constraints
[Ccal, Dcal, Ecal, Mcal,Ebar] = caligraphicMatricesExtended(umin,umax,xmin,xmax,ymin,ymax,dumin,dumax,N,p,n,m);  
uk = [u0 zeros(m,t)];
xk = [x0 A*x0+B*u0 zeros(n,t)];
yk = [y0 C*xk(:,2) zeros(p,t-1)];
H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
%L = Mcal*gamma + Ecal;           %no delta u and y constraints
L = Mcal*gamma + Ecal+Ebar*T;
W = -Dcal-Mcal*phi;
i = 0;
for k = 2:(k_sim+1)               %time starts at k = 0, xk is known for k = 1
    i = i+2;                      %used for indexing due to reference being twice as long;
    Rk = ref((i+1):(i+p*N));
    v = [uk(:,k-1); zeros(2*(N-1),1)];
    f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
    %[Uk,fval,exitflag] = quadprog(H,f,L,Ccal+W*xk(:,k),[],[],[],[],[],[]);  %no delta u and y constraints
    [Uk,fval,exitflag] = quadprog(H,f,L,Ccal+W*xk(:,k)+Ebar*v,[],[],[],[],[],[]);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
            break;                %optimization failed, break loop then plot results
        end
    end
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k) = C*xk(:,k);
end

%% plotting
close all;
%constraint = 'Altconstrained';
font = 18;
thickness = 2;
plottingFunction(constraint,font,thickness,t,xk,uk,yk,ref1,ref2);

