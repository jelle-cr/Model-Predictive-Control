function [Ccal, Dcal, Ecal, Mcal] = caligraphicMatricesTerminalSet(umin,umax,xmin,xmax,N,n,m,M_N,b_N)
    Ccal = b_N;
    b_i = [-umin; umax; -xmin; xmax];
    for i = 0:(N-1)
        Ccal = [b_i; Ccal];
    end
    
    Mcal = M_N;
    M_i = [zeros(m,n); zeros(m,n); -eye(n); eye(n)];
    for i = 1:(N-1)
        Mcal = blkdiag(M_i, Mcal);
    end
    Mcal = [zeros(height(M_i),width(Mcal)); Mcal];
    
    Dcal = [M_i; zeros((height(M_i)*(N-1)+height(M_N)),width(M_i))];
    
    Ecal = double.empty;
    E_i = [-eye(m); eye(m); zeros(n,m); zeros(n,m)];
    for i = 0:(N-1)
        Ecal = blkdiag(Ecal, E_i);
    end
    Ecal = [Ecal; zeros(height(M_N),width(Ecal))];
end