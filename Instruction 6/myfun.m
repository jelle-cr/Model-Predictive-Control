function OBJ = myfun(u,x0)
%------------------------------
% Global variables can be used to avoid duplicate definition
M = 1;
k0 = 0.33;
hd = 1.1;
Ts = 0.4;
A0 = [1,Ts;-Ts*k0/M,1-Ts*hd/M];
B0 = [0;Ts/M];
[nx,nu] = size(B0);
N = 5;
Q = 1*eye(nx);
R = 0.1*eye(nu);
[~,P,~] = dlqr(A0,B0,Q,R);
%------------------------------
x = zeros(nx,N+1);
x(:,1) = x0;
OBJ = 0;
for i = 1:N
    OBJ = OBJ + x(:,i)'*Q*x(:,i)+u(:,i)'*R*u(:,i);
    x(:,i+1) = [1,Ts;-Ts*k0/M*exp(-x(1,i)),1-Ts*hd/M]*x(:,i)+B0*u(:,i);
end
OBJ = OBJ + x(:,N+1)'*P*x(:,N+1);
end