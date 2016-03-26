function [Z, expVokBool] = getZ(M, MI, phi) 
    N = size(MI,1);
    b = zeros(N,1);
    b(N) = 1;
    A = speye(N) - M;
    Z = A\b;
    Z = full(Z);
    j = 0;
     while(1)
        j = j+1 ;     
        Zprev = Z;
        Z = getNextZ(Z, M, MI, b,phi);
        if mod(j,10) == 0
           %j
           norm(log(Z) - log(Zprev));
           if norm(log(Z) - log(Zprev)) < 0.0001
             break;
           end
           %norm(log(Z) - log(Zprev))
           if (j > 300)             
            break ;
           end
        end
     end  
     % Check feasible
     residual = norm(log(Z) - log(Zprev));
     minele = min(Z(:));
     expVokBool = 1;
     if minele == 0 || minele < OptimizeConstant.NUM_ERROR
       expVokBool = 0;
       fprintf('min zero');
     end 
     if residual > 10 || (~isreal(Z))
       expVokBool = 0;
     end
end

function Znext = getNextZ(Z, M, MI, b, phi)
    n = size(M,1);
    e = ones(n,1);
    Zsp =  spdiags(Z,0,n,n);
    U = MI * Zsp + MI * realmin;
    [~,~,sMI] = find(MI * realmin);
    [~,~,sU] = find(U);
    [i,j,sPhi] = find(phi);
    s = (sU - sMI) .^ sPhi;
    X = sparse(i,j,s,n,n);    
    Z =  (M .* X) * e;
    Znext = Z + b;   
end
