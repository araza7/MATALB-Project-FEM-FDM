function P_analytical = solution_1D(f, c1, c2, L1, L2, x)
    w = 2 * pi * f;    
    k1 = w / c1;       
    k2 = w / c2;       
    
    A=complex(zeros(4)); 
    A(1,:)=[1 exp(-1i*k1*L1) 0 0]; 
    A(2,:)=[exp(-1i*k1*L1)  exp(0)  -exp(0)  -exp(-1i*k2*(L2))]; 
    A(3,:) = [(-1i*k1)*exp(-1i*k1*L1), -(-1i*k1)*exp(0), -(-1i*k2)*exp(0), (-1i*k2)*exp(-1i*k2*L2)];
    A(4,:)=[0 0 (-1/1i/w)*(-1i*k2)*exp(-1i*k2*L2)*c2-exp(-1i*k2*L2)  -(-1/1i/w)*(-1i*k2)*c2-1]; 
    B=[1;0;0;0]; 
    X=A\B; 
    P_analytical=zeros(numel(x),1); 
    for ii=1:numel(x) 
        if(x(ii)<=L1 & x(ii)>=0) 
            P_analytical(ii)=X(1)*exp(-1i*k1*x(ii))+X(2)*exp(-1i*k1*abs(x(ii)-L1)); 
        elseif(x(ii)<=L1+L2) 
            P_analytical(ii)=X(3)*exp(-1i*k2*abs(x(ii)-L1))+X(4)*exp(-1i*k2*abs(x(ii)-L1-L2)); 
        else 
            P_analytical(ii)=NaN; 
        end 
    end 
end


