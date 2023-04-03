% Non-linear Model predictive control Homework instruction 6
% Jelle Cruijsen & Michiel Wind
yalmip('clear')
clear all, close all, clc

% rmpath("C:\Program Files\Mosek\10.0\toolbox\r2017a")
% addpath("C:\Program Files\Mosek\10.0\toolbox\r2017a")


% Simulation variables
t = 300; 
Ts = 0.4;
N = 5;

% Linearzed Model definition
% [A,B,C,D,sys] = modelselect('aircraft','discrete',Ts); % change
% [n,m] = size(B); 
% p = height(C);

% Initial conditions
x0 = [1;1];

% Performance tuning
Q = 10*eye(p);
R = 0.1*eye(m);

% Calculate terminal set
%% LMI (Linearized A and B?)
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

%% Construct the constraints
yalmip('clear')
nu = 1;
nx = 2;
A = -0.5*eye(2);
B = 1;
x = sdpvar(repmat(2,1,N+1),repmat(1,1,N+1));
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
cost = 0;
constraints = [];
for i = 1:N
cost = cost + x{i}'*Q*x{i} + u{i}'*R*u{i};
constraints = [constraints, x{i+1} == A*x{i} + B*u{i}];
constraints = [constraints, -1 <= u{i}<= 1, -5 <= x{i+1} <= 5];
end
controller = optimizer(constraints, cost,[],x{1},[u{:}]);
%%
% Change
% [phi, gamma] = predictionModel(A,B,N,n,m);
% omega = kron(eye(N),Q);       % Kronecker product
% psi = kron(eye(N),R);
% C_bar = kron(eye(N),C);


% constrained MPC Quadprog
constraint = 'constrained';    
[Ccal, Dcal, Ecal, Mcal,Ebar] = caligraphicMatricesExtended(umin,umax,xmin,xmax,ymin,ymax,dumin,dumax,N,p,n,m);
% [D,M,E,S,c] = DMESc_constraints(constr,N,B);
% Mcal = M;
uk = [u0 zeros(m,t)];
xk = [x0 A*x0+B*u0 zeros(n,t-1)];
yk = [y0 C*xk(:,2) zeros(p,t-1)];
H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
L = Mcal*gamma + Ecal+Ebar*T;
W = -Dcal-Mcal*phi;
i = 0;
for k = 2:t              %time starts at k = 0, xk is known for k = 1
    i = i+2;                      %used for indexing due to reference being twice as long;
    Rk = ref((i+1):(i+p*N));
    v = [uk(:,k-1); zeros(2*(N-1),1)];
    f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
    [Uk,fval,exitflag] = quadprog(H,f,L,Ccal+W*xk(:,k)+Ebar*v,[],[],[],[],[],[]);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
            break;  %optimization failed, break loop then plot results
        end
    end
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k+1) = C*xk(:,k+1);
end
