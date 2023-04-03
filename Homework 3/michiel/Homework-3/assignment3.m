% Assignment 3 Michiel Wind & Jelle Cruijsen
% clc;close all;clear all
% rmpath("C:\Program Files\Mosek\10.0\toolbox\r2017a")
% addpath("C:\Program Files\Mosek\10.0\toolbox\r2017a")

% Simulation parameters
t = 300;
N = 3;
k_sim = 18;

trials = 5; % Finding rho trials

% Model definition
A_sub = [1, 31.4159;
         0, 0.999];
B_sub = [0;0.01];

A = blkdiag(A_sub,A_sub);
B = blkdiag(B_sub,B_sub);
B_rho = [0;-0.005;0;0.005];
[nx,nu] = size(B);
x0 = [0;0.15;0;-0.15];

xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);

% Performance
Q = diag([0.1,10,0.1,10]);

% Constraints      
xmin = [-10;-10;-10;-10];
xmax = [10;10;10;10];
umin = [-5;-5];
umax = [5;5];

%% Simulation
yalmip('clear')
x = sdpvar(nx,N,'full');
u = sdpvar(nu,N-1,'full');
rho_var = sdpvar(1,N-1,'full');
objective = 0;
constraints = [];

for i = 1:N-1 % 
    objective = objective + x(:,i)'*Q*x(:,i); 
    constraints = [constraints, x(:,i+1)==A*x(:,i)+B*u(:,i)+B_rho*rho_var(i)];
    constraints = [constraints, xmin<=x(:,i)<=xmax, umin<=u(:,i)<=umax];
end

MPC_LPV = optimizer(constraints,objective,[],[x(:,1);rho_var'],u);
rho = ones(N-1,1)*sin(x0(1)-x0(3));
tolerance = 1e-16;
Uk_hist{1} = zeros(size(nu));
for k = 1:k_sim
        for l = 1:100
            Uk = MPC_LPV([xk(:,k);rho]);
            Uk = reshape(Uk,[nu*(N-1),1]);
            [phi,gamma,gamma_bar] = predict_model(A,B,B_rho,rho,N,Uk',xk(:,k));
%             phi
%             gamma
%             gamma_bar
%             break;
            Xk = phi*xk(:,k) + gamma*Uk + gamma_bar*rho;
            rho = newrho(Xk,nx);
            Uk_hist{l+1} = Uk;
            uk(:,k) = Uk(1:nu);
            if norm(Uk_hist{l+1}-Uk_hist{l}) < tolerance
            fprintf('Tolerance bound met on iteration %d \n',l)
            break
            end
        end
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k)+B_rho*rho(1);
end

%% Plotting
t = k_sim;
font = 18;
figure;
subplot(1,2,2);
plot(0:t,xk(:,1:length(xk))',LineWidth=2);
grid on;
ax = gca;
ax.FontSize=font-6;
xlim([0 t]);
xlabel('$k$',FontSize=font,Interpreter='latex');
ylabel('$x$',FontSize=font,Interpreter='latex');
legend({'$x_1$','$x_2$','$x_3$','$x_4$'},FontSize=font,Interpreter="latex");

subplot(1,2,1);
plot(0:t-1,uk(:,1:length(uk))',LineWidth=2);
grid on;
ax = gca;
ax.FontSize=font-6;
xlim([0 t]);
xlabel('$k$',FontSize=font,Interpreter='latex');
ylabel('$u$',FontSize=font,Interpreter='latex');
legend({'$u_1$','$u_2$'},FontSize=font,Interpreter="latex");
