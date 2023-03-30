function [CON,CONe] = mycon(u,x0,M_N,b_N)
%------------------------------
% Global variables can be used to avoid duplicate definition
M = 1;
k0 = 0.33;
hd = 1.1;
Ts = 0.4;
% A0 = [1,Ts;-Ts*k0/M,1-Ts*hd/M];
B0 = [0;Ts/M];
[nx,nu] = size(B0);
N = 5;
xmin = [-2.65;-2];
xmax = [+2.65;+2];
umin = -4.5;
umax = +4.5;
%------------------------------
CON = [];
CONe = [];
x = zeros(nx,N+1);
x(:,1) = x0;
for i = 1:N
    x(:,i+1) = [1,Ts;-Ts*k0/M*exp(-x(1,i)),1-Ts*hd/M]*x(:,i)+B0*u(:,i);
    CON = [CON; -x(:,i+1)+xmin; x(:,i+1)-xmax; -u(:,i)+umin; u(:,i)-umax];
end
CON = [CON; M_N*x(:,N+1)-b_N];
end