function Pk = rhoFunction(Xk,n)
    x1 = Xk(1:n:length(Xk),1);
    x2 = Xk(2:n:length(Xk),1);
    
    for i = 1:length(Xk)/n
        Pk(i,1) = exp(-x1(i));
    end
end