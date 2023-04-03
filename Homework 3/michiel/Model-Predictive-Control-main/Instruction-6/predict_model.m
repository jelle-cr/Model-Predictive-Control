function Xk = predict_model(A,B,rho,N,Uk,x0)
phi = cell(N,1);
gamma = cell(N,N);

zero_vec = zeros(size(A*B));
gamma_hist = ones(size(A*B));
for i = 1:N-1
gamma_hist = gamma_hist.*(A.*[0,0;rho(i),0])*B;
gamma_vec{i} = gamma_hist;
end


for j = 1:N
    for i = 1:N
        if i == j
        gamma{i,j} = B;
        elseif i > j
        gamma{i,j} = gamma_vec{i-1};
        else
        gamma{i,j} = zero_vec;
        end
    end
end

phi_hist = ones(size(A));
for i = 1:N
phi_hist = phi_hist * (A.*[0,0;rho(i),0]);
phi{i} = phi_hist;
end
phi = cell2mat(phi);
gamma = cell2mat(gamma);

Xk = phi * x0 + gamma*Uk;
end
