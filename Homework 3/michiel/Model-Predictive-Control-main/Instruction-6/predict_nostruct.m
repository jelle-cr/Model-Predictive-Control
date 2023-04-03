function Xk = predict_nostruct(A,B,rho,N,Uk,x0)
Xk = zeros(N,1);

Xk(1) =  A.*[0,0;rho(1),0]*x0+B*Uk(1);
Xk(2) =  A.*[0,0;rho(1),0]*A.*[0,0;rho(2),0]*x0+A.*[0,0;rho(2),0]*B*Uk(2)+A.*[0,0;rho(1),0]*B*Uk(1);
for k = 2:N
   for i = 1:k
    Xk(k) = Xk(k) +
   end
end

end

