close all;
clear all;
clc 

load 'stateMatrices.mat';

N = 10;
[n,m] = size(B);
x_ss = [10;1];
u_ss = 0.52;
x0 = [0.1; 1];          %original initial condition for x bar
%x0 = [5; -1];

shift = 'normal';       %switch between normal or shifted system       
shift = 'shifted';

%% Constraints          %defined on x bar
xmin = [-10;-1];
xmax = [5;2];
umin = -0.52;
umax = 0.43;

if(strcmp(shift,'shifted')~=1)
    xmin = xmin+x_ss;
    xmax = xmax+x_ss;
    umin = umin+u_ss;
    umax = umax+u_ss;
    x0 = x0+x_ss;
end

%% Necessary matrices
[Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m);
[phi, gamma] = predictionModel(A,B,N,n,m);

L = Mcal*gamma + Ecal;
W = -Dcal-Mcal*phi;

%% Cost function
Q = 1*eye(n);
R = 0.1;
[K,P,~] = dlqr(A,B,Q,R);                    % We need to find P through lyap
K = -K;

omega = blkdiag(kron(eye(N-1),Q), P);       % Kronecker product
psi = kron(eye(N),R);

G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*phi;

%% MPC
k_sim = 300;
xk = [x0 zeros(n,k_sim+1)];     %state at each time index k
uk = zeros(m,k_sim+1);          %input at each time index k
                                %From Summary_con_MPC slide:
                                %quadprog(G,Fx(k),L,Ccal+Wx(k)
opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym');
for k = 1:k_sim+1
    [Uk,~,exitflag] = quadprog(G,F*xk(:,k),L,Ccal+W*xk(:,k),[],[],[],[],[],opt);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
            break;
        end
    end
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

%% Plotting
font = 14;
figure
subplot(1,2,1);
plot(0:k_sim,xk(:,1:length(xk)-1)',LineWidth=2);
xlim([0 k_sim]);
xlabel('$k$',Interpreter="latex",Fontsize=font);
if(strcmp('shifted',shift))
    ylabel('$\bar{x}$',Interpreter="latex",Fontsize=font);
    legend({'$\bar{x}_1$ (voltage)','$\bar{x}_2$ (current)'},Interpreter="latex",Fontsize=font);
else
    ylabel('$x$',Interpreter="latex",Fontsize=font);
    legend({'$x_1$ (voltage)','$x_2$ (current)'},Interpreter="latex",Fontsize=font);
end

subplot(1,2,2);
stairs(0:k_sim,uk,LineWidth=2);
xlim([0 k_sim]);
ylim([umin-0.1 umax+0.1]);
xlabel('$k$',Interpreter="latex",Fontsize=font);
if(strcmp('shifted',shift))
    ylabel('$\bar{u}$  (duty cycle)',Interpreter="latex",Fontsize=font);
else
    ylabel('$u$ (duty cycle)',Interpreter="latex",Fontsize=font);
end
sgtitle('Constrained MPC for the buck converter',Interpreter="latex",Fontsize=font+4)

