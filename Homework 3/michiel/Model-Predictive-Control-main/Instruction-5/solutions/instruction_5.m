% 5LMB0 Model predictive control
% Instruction 5
clear all
close all

%% 2.1 Admissible and invariant set
% First-order model
Ac = 0;
Bc = 1;
Cc = 1;
Dc = 0;
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);

xmin = -2;
xmax = 2;
umin = -0.5;
umax = 0.5;

% Q,R,P design
Q = 1*eye(nx);
R = 1*eye(nu);
[K,P,~] = dlqr(A,B,Q,R); % In MATLAB documentation, u(k)=-Kx(k)
K_lqr = -K; % In this lecture, we use the form u(k)=+Kx(k)

% (a) terminal constraint set
% state constraints (A_x * x <= b_x)
X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
% input constraints (A_u * u <= b_u)
U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);
% constraints admissible set
CA_set = Polyhedron([-eye(nu);eye(nu)]*K_lqr,[-umin;umax]);
% or
CA_set = Polyhedron(U_set.A*K_lqr,U_set.b);
% input constraint admissible set
IA_set = CA_set&X_set;
% invariant set within input constraint admissible set
model = LTISystem('A', A+B*K_lqr,'Ts',Ts);
INV_set = model.invariantSet('X',X_set);

figure
subplot(2,2,1)
plot(X_set); title('X_{set}');
subplot(2,2,2)
plot(U_set); title('U_{set}');
subplot(2,2,3)
plot(CA_set); title('CA_{set}');
subplot(2,2,4)
plot(INV_set); title('INV_{set}');

%% 2.2 Admissible and invariant set
% Second-order model
Ac = [0 1;0 0];
Bc = [0;-1];
Cc = eye(2);
Dc = [0;0];
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);

xmin = [-10;-10];
xmax = [10;10];
umin = -1;
umax = 1;

% Q,R,P design
Q = 1*eye(nx);
R = 1*eye(nu);
[K,P,~] = dlqr(A,B,Q,R); % In MATLAB documentation, u(k)=-Kx(k)
K_lqr = -K; % In this lecture, we use the form u(k)=+Kx(k)

% state constraints (A_x * x <= b_x)
X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
% input constraints (A_u * u <= b_u)
U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);
% constraints admissible set
CA_set = Polyhedron([-eye(nu);eye(nu)]*K_lqr,[-umin;umax]);
% or
CA_set = Polyhedron(U_set.A*K_lqr,U_set.b);
% input constraint admissible set
IA_set = CA_set&X_set;

model = LTISystem('A', A+B*K_lqr,'Ts',Ts);

% (a) invariant set within state constraints set
INV_set = model.invariantSet('X',X_set);
% plot
figure;
subplot(2,3,1)
plot(INV_set); title('Inv. Set without C.A. set');
% plot trajectories
subplot(2,3,2)
plot(INV_set,'Color','red');
title('Trajectories, unsaturated u')
hold on
extr = INV_set.V';
clear u_list
for i = 1:size(extr,2)
    x = extr(:,i);
    for k = 0:200
        u = K_lqr*x(:,end);
        x = [x,A*x(:,end)+B*u];
    end
    plot(x(1,:),x(2,:),'.-','LineWidth',1);
end
% plot trajectories
subplot(2,3,3)
plot(INV_set,'Color','red');
title('Trajectories, saturated u')
hold on
extr = INV_set.V';
clear u_list
for i = 1:size(extr,2)
    x = extr(:,i);
    for k = 0:200
        u = min(umax,max(umin,K_lqr*(x(:,end))));
        x = [x,A*x(:,end)+B*u];
    end
    plot(x(1,:),x(2,:),'.-','LineWidth',1);
end

%% (b) invariant set within constraint admissible set
INV_set = model.invariantSet('X',IA_set);
% plot
subplot(2,3,4)
plot(INV_set); title('Inv. Set with C.A. set');
% plot trajectories
subplot(2,3,5)
plot(INV_set,'Color','red');
title('Trajectories, unsaturated u')
hold on
extr = INV_set.V';
clear u_list
for i = 1:size(extr,2)
    x = extr(:,i);
    for k = 0:200
        u = K_lqr*x(:,end);
        x = [x,A*x(:,end)+B*u];
        u_list(i,k+1) = u;
    end
    plot(x(1,:),x(2,:),'.-','LineWidth',2);
end
% plot trajectories
subplot(2,3,6)
plot(INV_set,'Color','red');
title('Trajectories, saturated u')
hold on
extr = INV_set.V';
clear u_list
for i = 1:size(extr,2)
    x = extr(:,i);
    for k = 0:200
        u = min(umax,max(umin,K_lqr*(x(:,end))));
        x = [x,A*x(:,end)+B*u];
        u_list(i,k+1) = u;
    end
    plot(x(1,:),x(2,:),'.-','LineWidth',2);
end

%% (d) Compute the quadratic lyaponov function
Z = Q+K_lqr'*R*K_lqr;
O = sdpvar(nx,nx);
Y = sdpvar(nu,nx);
L1 = [O,(A*O+B*Y)',O;(A*O+B*Y),O,zeros(size(O));O,zeros(size(O)),Z^-1]>=0;
L2 = O>=1e-9;
constraints = L1 + L2;
diagnostics = optimize(constraints);
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end
P = value(O)^-1;

c1 = 6;
[X,Y] = meshgrid(-c1:0.1:c1,-c1:0.1:c1);
clear LF
for i = 1:size(X,1)
    for j = 1:size(Y,1)
        LF(i,j) = [X(i,j);Y(i,j)]'*P*[X(i,j);Y(i,j)];
    end
end
figure
plot(INV_set,'Color','red');
hold on
contour(X,Y,LF,200)
title('Quadratic Lyaponuv Function level sets with the Invariant Set')

%% (f) feasible set of states
M_N = INV_set.A;
b_N = INV_set.b;

N_list = [1,3,5,8,12,15,19,25,30,40,50,75];
% N_list = [1,5,10,15];
figure; 
hold on;
for nn = 1:length(N_list)
    N = N_list(nn);
    [Phi, Gamma] = ABN2PhiGamma(A,B,N);
    [ W, L, c] = getWLc_TS( A, B, xmax, xmin, umax, umin, Gamma, Phi, M_N, b_N);
    Feasibleset_x0_U = Polyhedron('A',[-W L],'B',c);
    Feasibleset_x0(nn) = projection(Feasibleset_x0_U,1:2);
end
for nn = length(N_list):-1:1
    plot(Feasibleset_x0(nn),'color','yellow','Alpha',1-0.9/nn)
end
title('Feasible set of states for increasing N');
plot(INV_set)

%% (g) Terminal constraint for x(k+1) = Ax(k) + Bu(k)
system = LTISystem('A', A, 'B', B);
system.x.min = xmin;
system.x.max = xmax;
system.u.min = umin;
system.u.max = umax;
InvSet = system.invariantSet();
figure
plot(InvSet)
xlabel('x_1')
ylabel('x_2')
title('Invariant set of x(k+1)=Ax(k)+bu(k)')

%% 2.3 Constrained MPC Simulation
clear all
close all
% Second-order model
Ac = [0 1;0 0];
Bc = [0;-1];
Cc = eye(2);
Dc = [0;0];
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);
xmin = [-10;-10];
xmax = [10;10];
umin = -1;
umax = 1;
x0 = [-6;3];

% Q,R,P,Terminal Set
Q = 1*eye(nx);
R = 1*eye(nu);
[K,P,~] = dlqr(A,B,Q,R); % In MATLAB documentation, u(k)=-Kx(k)
K_lqr = -K; % In this lecture, we use the form u(k)=+Kx(k)

X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);
CA_set = Polyhedron(U_set.A*K_lqr,U_set.b);
IA_set = CA_set&X_set;
model = LTISystem('A', A+B*K_lqr,'Ts',Ts);
INV_set = model.invariantSet('X',IA_set);
M_N = INV_set.A;
b_N = INV_set.b;

% Compact formulation
N = 50;
[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,R,P,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;
[W, L, c] = getWLc_TS(A,B,xmax,xmin,umax,umin,Gamma,Phi,M_N,b_N);

% (a) Feasible set of states - [-W L][x0;Uk]<=c
Feasibleset_x0_U = Polyhedron('A',[-W L],'B',c);
Feasibleset_x0 = projection(Feasibleset_x0_U,1:2);

k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
% opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym')
for k = 1:k_sim
    [Uk,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[]);
    if exitflag ~= 1
        warning('exitflag quadprog =%d\n', exitflag)
        if exitflag == -2
            disp('Optimization problem is infeasible. \n')
        end
    end
    uk(:,k) = Uk(1:nu);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end
figure
subplot(1,2,1)
stairs(0:k_sim,xk')
xlabel('$k$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
title('constrained MPC')
axis([0 k_sim -5 10])
subplot(1,2,2)
stairs(0:k_sim-1,uk)
xlabel('$k$','Interpreter','latex')
ylabel('$u$','Interpreter','latex')
title('constrained MPC')
axis([0 k_sim -1 1])

figure
plot(Feasibleset_x0)
hold on
plot(INV_set)
plot(xk(1,:),xk(2,:),'o')
xlabel('x_1')
ylabel('x_2')
title('Trajectory with Terminal set')

% (b) Feasible set of states without terminal constraint
[ W, L, c] = getWLc( A, B, xmax, xmin, umax, umin, Gamma, Phi);
Feasibleset_x0_U = Polyhedron('A',[-W L],'B',c);
Feasibleset_x0 = projection(Feasibleset_x0_U,1:2);

k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
% opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym')
for k = 1:k_sim
    [Uk,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[]);
    if exitflag ~= 1
        warning('exitflag quadprog =%d\n', exitflag)
        if exitflag == -2
            fprint('Optimization problem is infeasible. \n')
        end
    end
    uk(:,k) = Uk(1:nu);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end
figure
subplot(1,2,1)
stairs(0:k_sim,xk')
xlabel('$k$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
title('constrained MPC')
axis([0 k_sim -10 10])
subplot(1,2,2)
stairs(0:k_sim-1,uk)
xlabel('$k$','Interpreter','latex')
ylabel('$u$','Interpreter','latex')
title('constrained MPC')
axis([0 k_sim -1 1])

figure
plot(Feasibleset_x0)
plot(xk(1,:),xk(2,:),'o')
xlabel('x_1')
ylabel('x_2')
title('Trajectory with Terminal set')
