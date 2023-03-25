function Pk = rhoFunction(Xk,n)
    x1 = Xk(1:4:length(Xk),1);
    x2 = Xk(2:4:length(Xk),1);
    x3 = Xk(3:4:length(Xk),1);
    x4 = Xk(4:4:length(Xk),1);
    
    for i = 1:length(Xk)/n
        Pk(i,1) = sin(x1(i)-x3(i));
    end
end