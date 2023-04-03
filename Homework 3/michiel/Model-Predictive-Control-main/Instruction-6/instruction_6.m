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
objective = objective + x(:,N+1)'*P*x(:,N+1); % Terminal cost 
constraints = [constraints, M_N*x(:,N+1)<=b_N]; % Terminal set constraints
MPC_Linear = optimizer(constraints,objective,[],x(:,1),u); % Define object such that you dont need Yalmip overhead every time instance

k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic
    Uk = MPC_Linear(xk(:,k)); % Get the optimal control sequence
    tk(:,k) = toc;
    uk(:,k) = Uk(1:nu);
    % Use discrete nonlinear model
    xk(:,k+1) = [1,Ts;-Ts*k0/M*exp(-xk(1,k)),1-Ts*hd/M]*xk(:,k)+B0*uk(:,k); % Simulate next step
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
    constraints = [constraints, x(:,i+1)==[1,Ts;-Ts*k0/M*exp(-x(1,i)),1-Ts*hd/M]*x(:,i)+B0*u(:,i)]; % Nonlinear constraints
    constraints = [constraints, xmin<=x(:,i)<=xmax, umin<=u(:,i)<=umax];
end
objective = objective + x(:,N+1)'*P*x(:,N+1); % Terminal cost
constraints = [constraints, M_N*x(:,N+1)<=b_N]; % Terminal set 
MPC_Nonlinear = optimizer(constraints,objective,[],x(:,1),u); % Constraints became nonlinear
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
% Initialization
% Design variables
yalmip('clear')
x = sdpvar(nx,N+1,'full');
u = sdpvar(nu,N,'full');
rho_var = sdpvar(1,N,'full');
objective = 0;
constraints = [];
% Plant model
A = [1, Ts;
    -Ts*(k0/M), 1-Ts*(hd/M)];

B = [0;Ts/M];
for i = 1:N
    objective = objective + x(:,i)'*Q*x(:,i)+u(:,i)'*R*u(:,i);
    constraints = [constraints, x(:,i+1)==[A(1,1),A(1,2);A(2,1)*rho_var(i),A(2,2)]*x(:,i)+B*u(:,i)]; % This constraint covers the prediction model
    constraints = [constraints, xmin<=x(:,i)<=xmax, umin<=u(:,i)<=umax];
end
    objective = objective + x(:,N+1)'*P*x(:,N+1); % Terminal cost 
    constraints = [constraints, M_N*x(:,N+1)<=b_N]; % Terminal set constraints
    MPC_LPV = optimizer(constraints,objective,[],[x(:,1);rho_var'],u); % what more inputs?

x0 = [-1;1];
k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
trials = 5;

% rho = ones(N,1)*exp(-xk(:,k));
rho = ones(N,1)*exp(-x0(1));
for k = 1:k_sim
%     rho = ones(N,1)*exp(-x0(1));
        for l = 1:10
        Uk = MPC_LPV([xk(:,k);rho]);
        if length(Uk) ~= 5
        error('Solver cant find a Uk')
        end
        Xk = predict_model(A,B,rho,N,Uk',xk(:,k)); % Problem is that Xk goes outside of constriant bound of x
        rho = newrho(Xk,nx);
        Uk_hist{l} = Uk;
        uk(:,k) = Uk(1:nu);
        % Caluclate norm(Uk(l=1)-Uk(l=10)) and check iff sufficient 
        end
    xk(:,k+1) = (A.*[0,0;rho(1),0])*xk(:,k)+B*uk(:,k); % Simulate next step (Use rho from optimization)
end
data(4).x = xk;
data(4).u = uk;

%% (e) Plot
figure
subplot(2,2,1)
hold on
grid on
stairs(0:k_sim,data(1).x(1,:),'b')
stairs(0:k_sim,data(2).x(1,:),'r')
stairs(0:k_sim,data(3).x(1,:),'g--')
stairs(0:k_sim,data(4).x(1,:),'m--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','LPV')
xlabel('k')
ylabel('x1')
subplot(2,2,2)
hold on
grid on
stairs(0:k_sim,data(1).x(2,:),'b')
stairs(0:k_sim,data(2).x(2,:),'r')
stairs(0:k_sim,data(3).x(2,:),'g--')
stairs(0:k_sim,data(4).x(2,:),'m--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','LPV')
xlabel('k')
ylabel('x2')
subplot(2,2,3)
hold on
grid on
stairs(0:k_sim-1,data(1).u,'b')
stairs(0:k_sim-1,data(2).u,'r')
stairs(0:k_sim-1,data(3).u,'g--')
stairs(0:k_sim-1,data(4).u,'m--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','LPV')
xlabel('k')
ylabel('u')
subplot(2,2,4)
hold on
grid on
plot(data(1).x(1,:),data(1).x(2,:),'b')
plot(data(2).x(1,:),data(2).x(2,:),'r')
plot(data(3).x(1,:),data(3).x(2,:),'g--')
plot(data(4).x(1,:),data(4).x(2,:),'m--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon','LPV')
xlabel('x1')
ylabel('x2')

%%
% figure
% hold on
% grid on
% stairs(0:k_sim-1,data(1).t,'b')
% stairs(0:k_sim-1,data(2).t,'r')
% stairs(0:k_sim-1,data(3).t,'g--')
% legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
% xlabel('k')
% ylabel('T')