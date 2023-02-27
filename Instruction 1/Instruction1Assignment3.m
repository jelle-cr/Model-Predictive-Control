close all;
clear all;
clc;

load carModelSecondOrder.mat
N = 5;
x0 = [3;1];

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

[nx,nu] = size(B);

Q = 1*eye(nx);
R = 1*eye(nu);
[K,P,~] = dlqr(A,B,Q,R); % In MATLAB documentation, u(k)=-Kx(k)
K_lqr = -K; % In this lecture, we use the form u(k)=+Kx(k)

%%LQR simulation
k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_lqr*xk(:,k);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end
figure
subplot(1,3,1)
plot(0:k_sim,xk)
title('No saturation');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');

% (e) LQR Stability
eig(A+B*K_lqr) % check if inside the unit circle

% (f) Input saturation
k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_lqr*xk(:,k);
    uk(:,k) = min(max(uk(:,k),-1*ones(nu,1)),1*ones(nu,1)); % saturation -[1,1]
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

subplot(1,3,2)
plot(0:k_sim,xk)
title('Saturation [-1,1]');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');

k_sim = 500;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_lqr*xk(:,k);
    uk(:,k) = min(max(uk(:,k),-0.1*ones(nu,1)),0.1*ones(nu,1)); % saturation -[0.1,0.1]
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

subplot(1,3,3)
plot(0:k_sim,xk)
title('Saturation [-0.1,0.1]');
legend({'$x_1$','$x_2$'},Interpreter='latex',FontSize=14);
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');

%%Prediction model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = double.empty;
gamma = zeros(N*nx,N*nu);
for i = 1:N
    phi = [phi; A^i];
    for j = 1:i
        gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = A^(i-j)*B;
    end
end

%%Cost function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = blkdiag(kron(eye(N-1),Q), P);
psi = kron(eye(N),R);

G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*phi;

K_mpc = -[eye(nu) zeros(nu,N-1)]*inv(G)*F;

k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_mpc*xk(:,k);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
plot(0:k_sim,xk)
title('Unconstrained MPC');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');


