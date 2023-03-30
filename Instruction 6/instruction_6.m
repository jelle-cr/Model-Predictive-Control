% 5LMB0 Model predictive control
% Instruction 6
clear all
close all
%% (a) Linearized discrete model around origin
M = 1;
k0 = 0.33;
hd = 1.1;
Ts = 0.4;
A0 = [1,Ts;-Ts*k0/M,1-Ts*hd/M];
B0 = [0;Ts/M];
[nx,nu] = size(B0);

xmin = [-2.65;-2];
xmax = [+2.65;+2];
umin = -4.5;
umax = +4.5;
N = 5;

x0 = [-1;1];
% x0 = [-2;2];

% Q,R,P,Terminal Set
Q = 1*eye(nx);
R = 0.1*eye(nu);
[K,P,~] = dlqr(A0,B0,Q,R);
K_lqr = -K;

X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);
CA_set = Polyhedron(U_set.A*K_lqr,U_set.b);
IA_set = CA_set&X_set;
model = LTISystem('A',A0+B0*K_lqr);
INVset = model.invariantSet('X',IA_set); % invariant set with CA set
M_N = INVset.A;
b_N = INVset.b;

%% (b) Linear MPC - Yalmip
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
objective = 0;
constraints = [];
for i = 1:N
    objective = objective + x(:,i)'*Q*x(:,i)+u(:,i)'*R*u(:,i);
    constraints = [constraints, x(:,i+1)==A0*x(:,i)+B0*u(:,i)];
    constraints = [constraints, xmin<=x(:,i)<=xmax, umin<=u(:,i)<=umax];
end
objective = objective + x(:,N+1)'*P*x(:,N+1);
constraints = [constraints, M_N*x(:,N+1)<=b_N];
MPC_Linear = optimizer(constraints,objective,[],x(:,1),u);

k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic
    Uk = MPC_Linear(xk(:,k));
    tk(:,k) = toc;
    uk(:,k) = Uk(1:nu);
    % Use discrete nonlinear model
    xk(:,k+1) = [1,Ts;-Ts*k0/M*exp(-xk(1,k)),1-Ts*hd/M]*xk(:,k)+B0*uk(:,k);
end
data(1).x = xk;
data(1).u = uk;
data(1).t = tk;

%% (c) Nonlinear MPC - YALMIP (solver:fmincon)
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
objective = 0;
constraints = [];
for i = 1:N
    objective = objective + x(:,i)'*Q*x(:,i)+u(:,i)'*R*u(:,i);
    constraints = [constraints, x(:,i+1)==[1,Ts;-Ts*k0/M*exp(-x(1,i)),1-Ts*hd/M]*x(:,i)+B0*u(:,i)];
    constraints = [constraints, xmin<=x(:,i)<=xmax, umin<=u(:,i)<=umax];
end
objective = objective + x(:,N+1)'*P*x(:,N+1);
constraints = [constraints, M_N*x(:,N+1)<=b_N];
MPC_Nonlinear = optimizer(constraints,objective,[],x(:,1),u);
% Check the solver of "MPC_Nonlinear": "Solver: FMINCON-STANDARD"

k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic
    Uk = MPC_Nonlinear(xk(:,k));
    tk(:,k) = toc;
    uk(:,k) = Uk(1:nu);
    % Use discrete nonlinear model
    xk(:,k+1) = [1,Ts;-Ts*k0/M*exp(-xk(1,k)),1-Ts*hd/M]*xk(:,k)+B0*uk(:,k);
end
data(2).x = xk;
data(2).u = uk;
data(2).t = tk;

%% (c) Nonlinear MPC - MATLAB fmincon
U_pre = zeros(nu,N);
k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic
    [Uk, fval, exitflag] = fmincon(@(u)myfun(u,xk(:,k)),U_pre,[],[],[],[],[],[],@(u)mycon(u,xk(:,k),M_N,b_N));
    tk(:,k) = toc;
    uk(:,k) = Uk(1:nu);
    U_pre = Uk;
    % Use discrete nonlinear model
    xk(:,k+1) = [1,Ts;-Ts*k0/M*exp(-xk(1,k)),1-Ts*hd/M]*xk(:,k)+B0*uk(:,k);
end
data(3).x = xk;
data(3).u = uk;
data(3).t = tk;

%% (d) iterative QP: Graded assignment work
clc
n = 2;
m = 1;

xmin = [-2.65; -2];
xmax = [2.65; 2];
umin = -4.5;
umax = 4.5;

omega = blkdiag(kron(eye(N-1),Q), P);       % Kronecker product
psi = kron(eye(N),R);

[Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m);

B = [0; 1/M];
t=20;
k_sim = t;
xk = [x0 zeros(n,k_sim+1)];     %state at each time index k
uk = zeros(m,k_sim);            %input at each time index k
pk = zeros(N,k_sim);
Pk = exp(-x0(1))*ones(N,1);
pkall = Pk;
Xk = zeros(n*(N),1);
Uk_prev = zeros(m*(N),1);
tolerance = 0.0000005;
Loop = 100;

A1 = [0 1; -k0/M*Pk(1) -hd/M];
A2 = [0 1; -k0/M*Pk(2) -hd/M];
A3 = [0 1; -k0/M*Pk(3) -hd/M];
A4 = [0 1; -k0/M*Pk(4) -hd/M];
A5 = [0 1; -k0/M*Pk(5) -hd/M];

[phi, gamma] = predictionModel(A1,A2,A3,A4,A5,B,N,n,m);
G = 2*gamma'*omega*gamma;
F = 2*gamma'*omega*phi;

L = Mcal*gamma + Ecal;
W = -Dcal-Mcal*phi;

opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym');
for k = 1:k_sim+1
    for l = 1:Loop
        [Uk,~,exitflag] = quadprog(G,F*xk(:,k),L,Ccal+W*xk(:,k),[],[],[],[],[],opt);
        if exitflag ~= 1
            warning('exitflag quadprog = %d\n', exitflag)
            if exitflag == -2
                sprintf('Optimization problem is infeasible.')
            end
        end
        Xk=phi*xk(:,k)+gamma*Uk;
        Pk = rhoFunction(Xk,n);
        pkall = [pkall Pk];
        A1 = [0 1; -k0/M*Pk(1) -hd/M];
        A2 = [0 1; -k0/M*Pk(2) -hd/M];
        A3 = [0 1; -k0/M*Pk(3) -hd/M];
        A4 = [0 1; -k0/M*Pk(4) -hd/M];
        A5 = [0 1; -k0/M*Pk(5) -hd/M];
        
        [phi, gamma] = predictionModel(A1,A2,A3,A4,A5,B,N,n,m);
        G = 2*gamma'*omega*gamma;
        F = 2*gamma'*omega*phi;
        
        L = Mcal*gamma + Ecal;
        W = -Dcal-Mcal*phi;
        if(norm(Uk-Uk_prev)<tolerance)
        end
        Uk_prev = Uk;
        if(l>10)
            warning('l = %d -> Slow Uk convergence, consider increasing tolerance\n', l)
        end
        break;
    end
    pk(:,k) = Pk;                       % Not necessary, but computed for analysis
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A1*xk(:,k)+B*uk(:,k);
    break;
end


%% (e) Plot
figure
subplot(2,2,1)
hold on
grid on
plot(0:k_sim,data(1).x(1,:),'b',LineWidth=2)
plot(0:k_sim,data(2).x(1,:),'r',LineWidth=2)
plot(0:k_sim,data(3).x(1,:),'g--',LineWidth=2)
plot(0:k_sim,xk(1,1:length(xk)-1),LineWidth=2)
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','Quasi LPV')
xlabel('k')
ylabel('x1')
subplot(2,2,2)
hold on
grid on
plot(0:k_sim,data(1).x(2,:),'b',LineWidth=2)
plot(0:k_sim,data(2).x(2,:),'r',LineWidth=2)
plot(0:k_sim,data(3).x(2,:),'g--',LineWidth=2)
plot(0:k_sim,xk(2,1:length(xk)-1),LineWidth=2)
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','Quasi LPV')
xlabel('k')
ylabel('x2')
subplot(2,2,3)
hold on
grid on
plot(0:k_sim-1,data(1).u,'b',LineWidth=2)
plot(0:k_sim-1,data(2).u,'r',LineWidth=2)
plot(0:k_sim-1,data(3).u,'g--',LineWidth=2)
plot(0:k_sim,uk,LineWidth=2)
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','Quasi LPV')
xlabel('k')
ylabel('u')
subplot(2,2,4)
hold on
grid on
plot(data(1).x(1,:),data(1).x(2,:),'b',LineWidth=2)
plot(data(2).x(1,:),data(2).x(2,:),'r',LineWidth=2)
plot(data(3).x(1,:),data(3).x(2,:),'g--',LineWidth=2)
plot(xk(1,1:length(xk)-1),xk(2,1:length(xk)-1),LineWidth=2)
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','Quasi LPV')
xlabel('x1')
ylabel('x2')

figure
hold on
grid on
plot(0:k_sim-1,data(1).t,'b',LineWidth=2)
plot(0:k_sim-1,data(2).t,'r',LineWidth=2)
plot(0:k_sim-1,data(3).t,'g--',LineWidth=2)
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('T')

% figure
% subplot(2,2,1)
% hold on
% grid on
% stairs(0:k_sim,data(1).x(1,:),'b',LineWidth=2)
% stairs(0:k_sim,data(2).x(1,:),'r',LineWidth=2)
% stairs(0:k_sim,data(3).x(1,:),'g--',LineWidth=2)
% legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
% xlabel('k')
% ylabel('x1')
% subplot(2,2,2)
% hold on
% grid on
% stairs(0:k_sim,data(1).x(2,:),'b',LineWidth=2)
% stairs(0:k_sim,data(2).x(2,:),'r',LineWidth=2)
% stairs(0:k_sim,data(3).x(2,:),'g--',LineWidth=2)
% legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
% xlabel('k')
% ylabel('x2')
% subplot(2,2,3)
% hold on
% grid on
% stairs(0:k_sim-1,data(1).u,'b',LineWidth=2)
% stairs(0:k_sim-1,data(2).u,'r',LineWidth=2)
% stairs(0:k_sim-1,data(3).u,'g--',LineWidth=2)
% legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
% xlabel('k')
% ylabel('u')
% subplot(2,2,4)
% hold on
% grid on
% plot(data(1).x(1,:),data(1).x(2,:),'b',LineWidth=2)
% plot(data(2).x(1,:),data(2).x(2,:),'r',LineWidth=2)
% plot(data(3).x(1,:),data(3).x(2,:),'g--',LineWidth=2)
% legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
% xlabel('x1')
% ylabel('x2')
% 
% figure
% hold on
% grid on
% stairs(0:k_sim-1,data(1).t,'b',LineWidth=2)
% stairs(0:k_sim-1,data(2).t,'r',LineWidth=2)
% stairs(0:k_sim-1,data(3).t,'g--',LineWidth=2)
% legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
% xlabel('k')
% ylabel('T')