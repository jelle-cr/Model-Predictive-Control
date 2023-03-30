close all;
clear all;
clc 

load 'stateMatrices.mat';

%A(1,2)=A(1,2)+0.101096;        %shrinks feasible set drasticly, determined
                                %before new A matrix was given
%A(1,2)=A(1,2)/10;

Ts = 2.5*10^(-6);       % Given as 2.5 micro seconds
t = 300;
N = 10;
[n,m] = size(B);
x_ss = [10;1];
u_ss = 0.52;
unew = (x_ss-A*x_ss);
unew = unew(1)/B(1)
u_ss = unew;
x0 = [0.1; 1];         % response 1, Original initial condition for x bar

x0 = [0.5; 2];  %stable
%x0 = [-3; -1];  %stable
x0 = [2; 0.5];  %stable
x0 = [-3.5; 1]; %unstable
%x0 = [-5; 2];  %unstable

%% Constraints          %defined on x bar
xmin = [-10;-1];
xmax = [5;2];
umin = -0.52;
umax = 0.43;

%% Cost matrix definition
Q = 1*eye(n);
R = 0.1;

%% LMI
O = sdpvar(n,n);                            % O is symmetric, therefore O = O'
Y = sdpvar(m,n);
Con1 = [O, (A*O+B*Y)', O, Y';
    (A*O+B*Y), O, zeros(n,n), zeros(n,m);
    O, zeros(n,n), Q^-1, zeros(n,m);
    Y, zeros(m,n), zeros(m,n), R^-1]>=0;
Con2 = O>=1e-9;
constraints = Con1 + Con2;
diagnostics = optimize(constraints);
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end
P = value(O)^-1;
K = value(Y)*value(O)^-1;

%[K,P,~] = dlqr(A,B,Q,R);       
%K = -K;
%% Terminal set                                 % From instruction 5
X_set = Polyhedron([-eye(n);eye(n)],[-xmin;xmax]);
CA_set = Polyhedron([-eye(m);eye(m)]*K,[-umin;umax]);
IA_set = CA_set&X_set;

model = LTISystem('A', A+B*K,'Ts',Ts);
INV_set = model.invariantSet('X',IA_set);
M_N = INV_set.A;
b_N = INV_set.b;

%% Required matrices
[Ccal, Dcal, Ecal, Mcal] = caligraphicMatricesTerminalSet(umin,umax,xmin,xmax,N,n,m,M_N,b_N);
[phi, gamma] = predictionModel(A,B,N,n,m);

L = Mcal*gamma + Ecal;
W = -Dcal-Mcal*phi;

omega = blkdiag(kron(eye(N-1),Q), P);       % Kronecker product
psi = kron(eye(N),R);

G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*phi;

Feasibleset_x0_U = Polyhedron('A',[-W L],'B',Ccal);
Feasibleset_x0 = projection(Feasibleset_x0_U,1:2);

%% MPC
k_sim = t;
xbar = [x0 zeros(n,k_sim+1)];     %state at each time index k
xbar2 = [x0 zeros(n,k_sim+1)];
xk = [x0+x_ss zeros(n,k_sim+1)];     %state at each time index k
ubar = zeros(m,k_sim+1);          %input at each time index k
uk = zeros(m,k_sim+1);          %input at each time index k
                                %From Summary_con_MPC slide:
                                %quadprog(G,Fx(k),L,Ccal+Wx(k)
opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym');
for k = 1:k_sim+1
    [Uk,~,exitflag] = quadprog(G,F*xbar(:,k),L,Ccal+W*xbar(:,k),[],[],[],[],[],opt);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
            break;
        end
    end
    ubar(:,k) = Uk(1:m);
    uk(:,k) = ubar(:,k) + u_ss;
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);

    xbar(:,k+1) = xk(:,k+1) - x_ss;            % Proper calculation method
    %xbar(:,k+1) = A*xbar(:,k)+B*ubar(:,k);      % Improper calculation method
end

%% Plotting
close all;
font = 20;
%terminalSetPlot(font,X_set,Feasibleset_x0,INV_set);
trajectoryPlot(font,k_sim,xk,uk,xmin+x_ss,xmax+x_ss,umin+u_ss,umax+u_ss,'original');
trajectoryPlot(font,k_sim,xbar,ubar,xmin,xmax,umin,umax,'shifted');

%% Extra plots
load stateResponse1.mat;
load stateResponse2.mat;
load stateResponse3.mat;
load stateResponse4.mat;

terminalSetPlotWithTrajectory(font,X_set,Feasibleset_x0,INV_set,x1,x2,x3,x4,xbar);

% N_list = [5:5:30];
% figure; 
% hold on;
% for nn = 1:length(N_list)
%     N = N_list(nn);
%     [Ccal, Dcal, Ecal, Mcal] = caligraphicMatricesTerminalSet(umin,umax,xmin,xmax,N,n,m,M_N,b_N);
%     [phi, gamma] = predictionModel(A,B,N,n,m);
%     L = Mcal*gamma + Ecal;
%     W = -Dcal-Mcal*phi;
%     Feasibleset_x0_U = Polyhedron('A',[-W L],'B',Ccal);
%     Feasibleset_x0(nn) = projection(Feasibleset_x0_U,1:2);
% end
% plot(X_set,'color','green','Alpha',0.5);
% for nn = length(N_list):-1:1
%     plot(Feasibleset_x0(nn),'color','red','Alpha',1-0.9/nn)
% end
% title('Feasible set of states for increasing N');
% plot(INV_set,'color','yellow');
