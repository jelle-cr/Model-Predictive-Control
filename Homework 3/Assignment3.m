close all;
clear all;
clc 

load stateMatrices.mat

Ts = 0.1;
N = 2;
t = 20;
[n,m] = size(B1);

x0 = [0; 0.15; 0; -0.15];

umin = [-5; -5];
umax = [5; 5];
xmin = [-10; -10; -10; -10];
xmax = [10; 10; 10; 10];
%xmin = 10000*xmin;
%xmax = 10000*xmax;
%umin = 10000*umin;
%umax = 10000*umax;

%% Prediction model
[phi, gamma, gamma_bar] = predictionModel(A,B1,B2,N-1,n,m);     % Cost function goes to N-1

%% Cost function
Q = diag([0.1,10,0.1,10]);
omega = kron(eye(N-1),Q);       % Kronecker product
G = 2*gamma'*omega*gamma;
F = 2*gamma'*omega*phi;
F_bar = 2*gamma'*omega*gamma_bar;

%% Constraints                  % Need alteration
[Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N-1,n,m);        % Also to N-1
L = Mcal*gamma + Ecal;
W = -Dcal-Mcal*phi;
W_bar = -Mcal*gamma_bar;

% Feasible set
%Feasibleset_x0_U = Polyhedron('A',[-W L],'B',Ccal);
%Feasibleset_x0 = projection(Feasibleset_x0_U,1:2);

%% MPC
k_sim = t;
xk = [x0 zeros(n,k_sim+1)];     %state at each time index k
uk = zeros(m,k_sim);            %input at each time index k
pk = zeros(N-1,k_sim);
Pk = sin(x0(1)-x0(3))*ones(N-1,1);
Xk = zeros(n*(N-1),1);
Uk_prev = zeros(m*(N-1),1);
tolerance = 0.0005;
Loop = 100;

opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym');
for k = 1:k_sim+1
    for l = 1:Loop
        [Uk,~,exitflag] = quadprog(G,F*xk(:,k)+F_bar*Pk,L,Ccal+W*xk(:,k)+W_bar*Pk,[],[],[],[],[],opt);
        if exitflag ~= 1
            warning('exitflag quadprog = %d\n', exitflag)
            if exitflag == -2
                sprintf('Optimization problem is infeasible.')
                break;
            end
        end
        Xk=phi*xk(:,k)+gamma*Uk+gamma_bar*Pk;
        Pk = rhoFunction(Xk,n);
        if(norm(Uk-Uk_prev)<tolerance)
            break;
        end
        Uk_prev = Uk;
        if(l>10)
            warning('l = %d -> Slow Uk convergence, consider increasing tolerance\n', l)
        end
    end
    pk(:,k) = Pk;                       % Not necessary, but computed for analysis
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A*xk(:,k)+B1*uk(:,k)+B2*Pk(1);
end

%% Plotting                 % Choose between stairs and plot
font = 18;
figure;
subplot(1,2,2);
plot(0:t,xk(:,1:length(xk)-1)',LineWidth=2);        
grid on;
ax = gca;
ax.FontSize=font-6;
xlim([0 t]);
xlabel('$k$',FontSize=font,Interpreter='latex');
ylabel('$x$',FontSize=font,Interpreter='latex');
legend({'$x_1$','$x_2$','$x_3$','$x_4$'},FontSize=font,Interpreter="latex");

subplot(1,2,1);
plot(0:t,uk(:,1:length(uk))',LineWidth=2);
grid on;
ax = gca;
ax.FontSize=font-6;
xlim([0 t]);
xlabel('$k$',FontSize=font,Interpreter='latex');
ylabel('$u$',FontSize=font,Interpreter='latex');
legend({'$u_1$','$u_2$'},FontSize=font,Interpreter="latex");

sgtitle(sprintf('Simulation for $N$ = %d',N),FontSize=font,Interpreter="latex")

%figure;
%plot(Feasibleset_x0);          % Doesn't work properly due to us having 4 states
