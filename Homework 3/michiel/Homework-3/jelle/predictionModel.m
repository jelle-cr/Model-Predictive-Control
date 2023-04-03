function [phi, gamma, gamma_bar] = predictionModel(A,B1,B2,N,n,m)
    phi = double.empty;
    gamma = zeros(N*n,N*m);
    gamma_bar = zeros(N*n,N);
    for i = 1:N
        phi = [phi; A^i];
        for j = 1:i
            gamma((i-1)*n+1:n*i,(j-1)*m+1:j*m) = A^(i-j)*B1;
            gamma_bar((i-1)*n+1:n*i,(j-1)+1:j) = A^(i-j)*B2;
        end
    end
end