function [phi,gamma,gamma_bar] = predict_model(A,B,B_rho,rho,N,Uk,x0)
phi = cell(N-1,1);
gamma = cell(N-1,N-1);
gamma_bar = cell(N-1,N-1);

% Construct gamma
zero_vec = zeros(size(A*B));
for j = 1:N-1
    for i = 1:N-1
        if i == j
        gamma{i,j} = B;
        elseif i > j
        gamma{i,j} = A^((i-1)-(j-1)) * B;
        else
        gamma{i,j} = zero_vec;
        end
    end
end

% Construct gamma_bar
zero_vec = zeros(size(A*B_rho));
for j = 1:N-1
    for i = 1:N-1
        if i == j
        gamma_bar{i,j} = B_rho;
        elseif i > j
        gamma_bar{i,j} = A^((i-1)-(j-1)) * B_rho;
        else
        gamma_bar{i,j} = zero_vec;
        end
    end
end

for i = 1:N-1
phi{i} = A^i;
end
phi = cell2mat(phi);
gamma = cell2mat(gamma);
gamma_bar = cell2mat(gamma_bar);
end
