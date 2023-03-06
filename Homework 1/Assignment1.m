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

ref = zeros(2,t+2*N);
for i = 1:width(ref)
    if(i<201)
        ref(1,i) = 10;
    else
        ref(1,i) = 8;
    end
    if(i<101)
        ref(2,i) = 0;
    else
        ref(2,i) = 5;
    end
end
% ref(1,1:200) = 10;
% ref(1,201:end) = 8;
% ref(2,1:100) = 0;
% ref(2,101:end) = 5;
% ref = reshape(ref,2*(t+2*N),1);

%% Constraints and matrices
constr.statelb = [-25;-25;-pi/20;-pi/2;-15]; % v,w,q,theta,h
constr.stateub =  [25;25;pi/20;pi/2;15]; % v,w,q,theta,h
constr.initialstatelb = -abs([x0; 0]); % v,w,q,theta,h
constr.initialstateub = abs([x0; 0]);
constr.terminalstatelb = [-25;-25;-pi/20;-pi/2;-15]; % Change if needed after monday
constr.terminalstateub = [25;25;pi/20;pi/2;15];
constr.inputlb = [-20;-20];
constr.inputub =  [20;20];
constr.deltainputlb = [-15;-15];
constr.deltainputub =  [15;15];

Q = 10*eye(p);
R = 0.1*eye(m);

[phi, gamma] = predictionModel(A,B,N,n,m);
omega = kron(eye(N),Q);       % Kronecker product
psi = kron(eye(N),R);
C_bar = kron(eye(N),C);

T = diag(ones(N*m,1));
nT = size(T,1);
T(2:nT+1:end) = -1;

[W,L,c,S] = getWLcS(constr,N,B,gamma,phi,T);

%% Question 4.2 (design a model predictive controller)
k_sim = t;
%% unconstrained
xk = [x0 zeros(n,t)];
yk = [y0 zeros(p,t)];
uk = [zeros(m,1) u0 zeros(m,t-1)];
xk(:,2) = A*xk(:,1)+B*uk(:,1); % Calculate x1 (because we know u0)
yk(:,2) = C*xk(:,1);
v = [u0; zeros(2*(N-1),1)];
Rk = double.empty;
for i = 0:(N-1)
   Rk = [Rk; ref(1,1+i)];
   Rk = [Rk; ref(2,1+i)];
end
G = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
F = 2*(gamma'*C_bar'*omega*C_bar*phi*x0-gamma'*C_bar'*omega*Rk-T'*psi*v);
M = zeros(m,2*N);
M(1:m,1:m) = eye(m);
for k = 1:k_sim
    Rk = double.empty;
    for i = 0:(N-1)
       Rk = [Rk; ref(1,k+i)];
       Rk = [Rk; ref(2,k+i)];
    end
    uk(:,k) = -M*inv(G)*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v); 
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k+1) = C*xk(:,k);
    v = [uk(:,k) ; zeros(2*(N-1),1)];
end
%% constrained
% xk = [x0 zeros(n,k_sim)];       %state at each time index k
% yk = [y0 zeros(p,k_sim)];       %output at each time index k
% uk = [zeros(size(u0)),u0,zeros(m,t-1)];
% xk(:,2) = A*xk(:,1)+B*uk(:,1); % Calculate x1 (because we know u0)
% yk(:,2) = C*xk(:,1);
% v = [u0 ; zeros(2*(N-1),1)];
% H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
% c = c+S*v; % update c to new constraints
% check = 0;
% for k = 2:t
%     [Uk,fval,exitflag] = quadprog(H,f,L,c+W*xk(:,k),[],[],[],[],[],[]);
%     if exitflag ~= 1
%         warning('exitflag quadprog = %d\n', exitflag)
%         if exitflag == -2
%             sprintf('Optimization problem is infeasible.')
%         end
%     end
%     uk(:,k) = Uk(1:m);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k+1) = C*xk(:,k);
%     v = [uk(:,k) ; zeros(2*(N-1),1)];
%     c = c+S*v; 
%     f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(k:(k-1+2*N))-2*T'*psi*v);
%     check = check+1
% end
% % 
figure()
subplot(1,2,1)
stairs(0:t,yk(1,:))
xlabel('$k$','Interpreter','latex');
ylabel('$v [ft/sec]$','Interpreter','latex');

subplot(1,2,2)
stairs(0:t,yk(2,:))
xlabel('$k$','Interpreter','latex');
ylabel('$h [ft/sec]$','Interpreter','latex');
sgtitle('Outputs')

figure()
subplot(1,2,1)
stairs(0:t,uk(1,:))
xlabel('$k$','Interpreter','latex');
ylabel('$e$','Interpreter','latex');

subplot(1,2,2)
stairs(0:t,uk(2,:))
xlabel('$k$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex');
sgtitle('Inputs')

ref = reshape(ref,[2,t+2*N]);
figure()
subplot(1,2,1)
stairs(0:t-1,ref(1,1:t)')
xlabel('$k$','Interpreter','latex');
ylabel('$v ref$','Interpreter','latex');
ylim([min(ref(1,:))-1 max(ref(1,:))+1]);

subplot(1,2,2)
stairs(0:t-1,ref(2,1:t)')
xlabel('$k$','Interpreter','latex');
ylabel('$h ref$','Interpreter','latex');
sgtitle('reference');
ylim([min(ref(2,:))-1 max(ref(2,:))+1]);