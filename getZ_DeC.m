function [Z, expVokBool] = getZ_DeC(M,B)
    global mu;
    
    % get sacles
    MI = sparse(M); 
    MI(find(M)) = 1;
    a = mu;
    k = 1 ./ mu;
    phi = sparse((k * a') .* MI);    
    A = speye(size(M)) - M;
    Z = A\B;
    Z = full(Z);
    j = 0;
     while(1)
        j = j+1 ;     
        Zprev = Z;
        Z = getNextZ(Z, M, MI, B, phi);
        if mod(j,10) == 0
            %j
            %norm(log(Z) - log(Zprev))
           if norm(log(Z) - log(Zprev)) < 0.0001
             break;
           end
           if (j > 200)             
            break ;
           end
        end
     end
    
     % Check feasible
     minele = min(Z(:));
     expVokBool = 1;
     if minele == 0 || minele < OptimizeConstant.NUM_ERROR
       expVokBool = 0;
     end 
     
     Zabs = abs(Z); % MAYBE SET TO VERY SMALL POSITIVE VALUE? 
     D = Z - Zprev;
     resNorm = norm(D(:));
     if resNorm > OptimizeConstant.RESIDUAL
       expVokBool = 0;
     end    
end

function Znext = getNextZ(Z, M, MI, B, phi)
    ndests = size(B,2);
    for d = 1:ndests
       U = sparse(bsxfun(@times,Z(:,d)',MI));
       X = MI;
       X(find(MI)) =  U(find(MI)) .^ phi(find(MI));
       Z(:,d) = diag(M * X');
    end
    Znext = Z + B;     
end