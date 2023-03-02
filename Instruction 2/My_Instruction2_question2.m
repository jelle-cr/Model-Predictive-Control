clear all;
close all;
clc;

N = 100;

%% Choose order of model
load modelSecondOrder.mat   %modelFirstOrder.mat or modelSecondOrder.mat
%x0 = 1;                    %first order model
x0 = [3;1];                 %second order model
[n,m] = size(B);

%% Question 2.2a
%[CcalTest, DcalTest, EcalTest, McalTest] = caligraphicMatricesTest(umin,umax,xmin,xmax,N,n,m);
[Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m);

%% Question 2.2b
[phi, gamma] = predictionModel(A,B,N,n,m);
L = Mcal*gamma + Ecal;
W = -Dcal-Mcal*phi;

%% Cost function
Q = 1*eye(n);
R = 1*eye(m);
[K,P,~] = dlqr(A,B,Q,R);                    % In MATLAB documentation, u(k)=-Kx(k)

omega = blkdiag(kron(eye(N-1),Q), P);       %Kronecker product
psi = kron(eye(N),R);

G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*phi;

%% Question 2.2c
k_sim = 500;
xk = [x0 zeros(n,k_sim)];       %state at each time index k
uk = zeros(m,k_sim);            %input at each time index k
Uk = zeros(m*N,k_sim);          %predicted inputs at each time index k for N future time steps
                                %From Summary_con_MPC slide:
                                %quadprog(G,Fx(k),L,Ccal+Wx(k)
opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym');
for k = 1:k_sim
    [U,~,exitflag] = quadprog(G,F*xk(:,k),L,Ccal+W*xk(:,k),[],[],[],[],[],opt);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
        end
    end
    Uk(:,k) = U;
    uk(:,k) = U(1:m);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
subplot(1,2,1);
plot(0:k_sim,xk');
xlabel('$k$','Interpreter','latex',Fontsize=14);
ylabel('$x$','Interpreter','latex',Fontsize=14);
legend('x_1','x_2',Fontsize=14);
xlim([0 k_sim]);

subplot(1,2,2);
stairs(0:k_sim-1,uk);
xlabel('$k$','Interpreter','latex',Fontsize=14);
ylabel('$u$','Interpreter','latex',Fontsize=14);
xlim([0 k_sim]);
ylim([umin-0.1 umax+0.1]);
sgtitle('Constrained MPC')

%% Question 2.2d
figure
hold on
for k = 1:k_sim
    stairs(k-1:k+N-2,Uk(:,k))           %predicted input Uk
end
stairs(0:k_sim-1,uk,'k',LineWidth=2);   %actual input uk
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 k_sim]);
ylim([umin-0.1 umax+0.1]);
title('Open loop predicted trajecotry of input')