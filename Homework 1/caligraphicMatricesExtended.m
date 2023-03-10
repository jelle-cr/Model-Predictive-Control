function [Ccal, Dcal, Ecal, Mcal,Ebar] = caligraphicMatricesExtended(umin,umax,xmin,xmax,ymin,ymax,dumin,dumax,N,p,n,m)
    Ccal = double.empty;
    for i = 0:(N-1)
        b = [-umin; umax; -xmin; xmax; -ymin; ymax; -dumin; dumax];
        Ccal = [Ccal; b];
    end
    Ccal = [Ccal; -xmin; xmax; -ymin; ymax];
    
    Mcal = double.empty;
    M_i = [zeros(2*m,n); -eye(n); eye(n); 0 1 0 -7.74; 0 -1 0 7.74; zeros(2*m,n)];
    for i = 1:(N-1)
        Mcal = blkdiag(Mcal, M_i);
    end
    Mcal = blkdiag(Mcal, [-eye(n); eye(n); 0 1 0 -7.74; 0 -1 0 7.74]);
    Mcal = [zeros(height(M_i),width(Mcal)); Mcal];
    
    Dcal = [M_i; zeros((height(M_i)*(N-1)+n*2+2),width(M_i))];  %N-1 times the rows of M_i + the rows of M_N
    
    Ecal = double.empty;
    E_i = [-eye(m); eye(m); zeros(2*n,m); zeros(6,m)];
    for i = 0:(N-1)
        Ecal = blkdiag(Ecal, E_i);
    end
    Ecal = [Ecal; zeros(n*2+2,width(Ecal))];

    Ebar = double.empty;
    E_i = [zeros(2*n,m); zeros(6,m); -eye(m); eye(m)];
    for i = 0:(N-1)
        Ebar = blkdiag(Ebar, E_i);
    end
    Ebar = [Ebar; zeros(height(xmax)+height(xmin)+height(ymax)+height(ymin),width(Ebar))];
end