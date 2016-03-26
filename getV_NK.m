function [Z, nI, expVokBool] = getV_NK(U, M, MI, phi) 
    global mu;
    N = size(MI,1);
    b = zeros(N,1);
    b(N) = 1;
    Ux = spdiags(1./mu,0,N,N)* U;  
    A = speye(N) - M;
    Z = A\b;
    V = log(Z);
    V_CI = V;
    Vprev = 0;
    j = 0;
    while(1)
        j = j+1 ;
        u = norm(V_CI - Vprev);
        Vprev = V;
        if u >OptimizeConstant.TAU_SWITCH_NFXP
            [V,V_CI] = getNextV_CI(V, Ux, MI, b,phi);
        else
            [V,V_CI] = getNextV_NK(V, Ux, MI, b,phi);
        end             
        if norm(V_CI - Vprev) < OptimizeConstant.TAU_CONVERGE_NFXP
            break;
        end
        if (j > 1300)             
            break ;
        end
    end 
     nI = j;
     residual  = norm(V_CI - Vprev);
     % Check feasible
     expVokBool = 1;
     if max(V) > 100 
       expVokBool = 0;
       fprintf('Large value functions');
     end 
     if residual > 0.01 || (~isreal(V))
       expVokBool = 0;
     end
     V = V .* mu;
end

function [Vnext,Vnext_CI] = getNextV_CI(V, Ux, MI, b, phi) % Newton Karotovic iteration
    N = size(Ux,1);
    e = ones(N,1);
    I = speye(N);
    Vsp =  spdiags(V,0,N,N);
    X = Ux + phi * Vsp;
    X(find(MI)) = exp(X(find(MI)));
    Vnext = log(X*e + b);
    Vnext_CI = Vnext;
end

function [Vnext,Vnext_CI] = getNextV_NK(V, Ux, MI, b, phi) % Newton Karotovic iteration
    N = size(Ux,1);
    e = ones(N,1);
    I = speye(N);
    Vsp =  spdiags(V,0,N,N);
    X = Ux + phi * Vsp;
    X(find(MI)) = exp(X(find(MI)));
    Y =  spdiags(exp(-V),0,N,N) * (X .* phi);
    Vnext_CI = log(X*e + b);
    Vnext_NK = V - (I - Y)\(V - Vnext_CI);
    if norm(Vnext_NK - V)>norm(Vnext_CI - V)
        Vnext = Vnext_NK;
    else
        Vnext = Vnext_CI;
    end
    
end