close all;
clear all;
clc;

load carModelSecondOrder.mat
N = 5;

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

P = 1*ones(size(A));%rand(size(A));
Q = 2*ones(size(A));%rand(size(A));
R = 3*ones(width(B),width(B));%rand(height(B),height(B));

[nx,nu] = size(B);

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
%omega = kron(eye(N-1),Q);
%omega = [omega zeros(height(omega),width(P)); zeros(height(P),width(omega)) P];
omega = blkdiag(kron(eye(N-1),Q), P);
psi = kron(eye(N),R);

G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*phi;

Kmpc = -[eye(nu) zeros(nu,N-1)]*inv(G)*F;


