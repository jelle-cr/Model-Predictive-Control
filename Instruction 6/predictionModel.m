function [phi, gamma] = predictionModel(A1,A2,A3,A4,A5,B,N,n,m)
    phi = [A1;A2*A1;A3*A2*A1;A4*A3*A2*A1;A5*A4*A3*A2*A1];
    zer = [0;0];
    gamma = [B zer zer zer zer; A2*B B zer zer zer; A3*A2*B A2*B B zer zer; A4*A3*A2*B A3*A2*B A2*B B zer; A5*A4*A3*A2*B A4*A3*A2*B A3*A2*B A2*B B];
    %phi = double.empty;
    %gamma = zeros(N*n,N*m);
%     for i = 1:N
%         phi = [phi; A^i];
%         for j = 1:i
%             gamma((i-1)*n+1:n*i,(j-1)*m+1:j*m) = A^(i-j)*B;
%         end
%     end
end