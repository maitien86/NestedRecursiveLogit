%   Get analytical Hessian matrix
%   Link size is included
%%
function [LL, grad, Hessian, Hs] = getAnalyticalHessian_Znested()
    
    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global isLinkSizeInclusive;
    global lastIndexNetworkState;
    global mu;
    % ----------------------------------------------------
    % If Link size is included
    % MU IS NORMALIZED TO ONE
    [lastIndexNetworkState, maxDest] = size(incidenceFull);
    
    Mfull = getM(Op.x, isLinkSizeInclusive);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(Op.x, isLinkSizeInclusive);
    UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    
    %% Set LL value
    LL = 0;
    grad = zeros(1, Op.n);
    
    %% Initialize
    M = MregularNetwork;
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    U = UregularNetwork;
    U(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    U(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    for i = 1:Op.n
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end

    %% b  and B matrix:
    N = size(M,1);
    b = sparse(zeros(N,1));
    b(N) = 1;
    B = sparse(zeros(N, maxDest - lastIndexNetworkState));
    B(N,:) = ones(1,maxDest - lastIndexNetworkState);
    for i = 1: maxDest - lastIndexNetworkState
        B(1:lastIndexNetworkState,i) = Mfull(:, i+lastIndexNetworkState);
    end
    
    %% Compute Z by iterative method: 
    % B = B(:,1:4);
    [Z, expVokBool]   = getZ(M, B);
    if (expVokBool == 0)
            LL = OptimizeConstant.LL_ERROR_VALUE;
            grad = ones(Op.n,1);
            disp('The parameters not fesible')
            return; 
    end
    %% Compute V
    V = log(Z);
    V = (bsxfun(@times,mu,V));

    %% Get gradient  Z
    gradExpV = objArray(Op.n);
    MI = sparse(M); 
    MI(find(M)) = 1;
    a = mu;
    k = 1 ./ mu;
    kf = k(1:lastIndexNetworkState,1);
    phi = sparse((k * a') .* MI);    
    for i = 1:Op.n
        gradZi = zeros(size(B));
        gradt = Atts(i).Value(:,lastIndexNetworkState+1 : maxDest) .* Mfull(:,lastIndexNetworkState+1 : maxDest);
        gradt = bsxfun(@times,kf,gradt);
        gradt(lastIndexNetworkState+1,:) = sparse(zeros(1,maxDest - lastIndexNetworkState));        
        gradt = sparse(gradt);
        gradM = M .* (AttLc(i).Value);
        gradM = bsxfun(@times,k,gradM);
        for d = 1: size(B,2)
            % Compute P, N
            Zd = sparse(bsxfun(@times,Z(:,d)',MI));
            X = MI;
            X(find(MI)) =  Zd(find(MI)) .^ (phi(find(MI))-1);
            N = gradM .* X;
            P = phi .* (M .* X);
            gradZi(:,d) = (speye(size(M)) - P)\(N * Z(:,d) + gradt(:,d));  
        end
        gradExpV(i).value = gradZi;  
    end
    
    %% Compute gradient of Vs
    gradV = objArray(Op.n);
    for i = 1:Op.n
        X = gradExpV(i).value ./ Z;
        gradV(i).value = sparse(bsxfun(@times,mu,X));
    end

    %% Evaluate hessian
    hessianM = objMatrix(Op.n, Op.n);
    hessianB = objMatrix(Op.n, Op.n);
    hessianZ = objMatrix(Op.n, Op.n);
    Hessian = zeros(Op.n);
    H = zeros(Op.n);
    Hs = objArray(nbobs);
       
    %% Get second order derivative of M and B
    for i = 1: Op.n
        for j = 1: Op.n
            u = M .* AttLc(i).Value .* AttLc(j).Value;
            hessianM(i,j).value = bsxfun(@times,k .* k, u);
            u = Atts(i).Value(:,lastIndexNetworkState+1 : maxDest) .*  Atts(j).Value(:,lastIndexNetworkState+1 : maxDest) .* Mfull(:,lastIndexNetworkState+1 : maxDest);
            u = bsxfun(@times,kf .* kf,u);
            u(lastIndexNetworkState+1,:) = sparse(zeros(1,maxDest - lastIndexNetworkState));        
            hessianB(i,j).value = sparse(u);            
            
        end
    end
    %% Get Hessian
    for i = 1: Op.n
        for j = 1: Op.n
            Hes = zeros(size(B));     
            gradMi = M .* (AttLc(i).Value);
            gradMi = bsxfun(@times,k,gradMi);
            gradMj = M .* (AttLc(j).Value);
            gradMj = bsxfun(@times,k,gradMj);
            for d = 1:size(B,2)
                Zd = sparse(bsxfun(@times,Z(:,d)',MI));
                X = MI;
                X(find(MI)) =  Zd(find(MI)) .^ (phi(find(MI))-1);  
                Gd = hessianM(i,j).value .* X;
                gd = (gradExpV(i).value(:,d) .* gradExpV(j).value(:,d)) ./ Z(:,d);
                Qd = phi .* (phi - 1) .* (M .* X);
                Pd = phi .* (M .* X);
                Ui = phi .* (gradMi .* X);
                Uj = phi .* (gradMj .* X);
                Hes(:,d) = (speye(size(M)) - Pd)\(Gd * Z(:,d) + Ui * gradExpV(j).value(:,d) + Uj * gradExpV(i).value(:,d) + Qd * gd +  hessianB(i,j).value(:,d));  
            end
            hessianZ(i,j).value = Hes;
        end
    end
    %% Compute second derivetive of V(.)
    hessianV = objMatrix(Op.n, Op.n);
    for i = 1:Op.n
        for j = 1: Op.n
            X = hessianZ(i,j).value ./ Z;
            Ui =  gradExpV(i).value ./ Z;
            Uj =  gradExpV(j).value ./ Z;
            hessianV(i,j).value = bsxfun(@times,mu, X - Ui .* Uj); 
        end
    end
    %% Compute the LL, grad and hessian    
    gradVd = zeros(size(Z,1),Op.n);
    hesVd = objMatrix(Op.n, Op.n);
    for n = 1:nbobs    
        dest = Obs(n, 1);
        orig = Obs(n, 2);       
        Vd = V(:,dest - lastIndexNetworkState);    
        lnPn = 0;
        for i = 1: Op.n
             Gradient(n,i) = 0;
             gradVd(:,i) =  gradV(i).value (:,dest - lastIndexNetworkState);
             for j = 1:Op.n
                 hesVd(i,j).value = hessianV(i,j).value(:,dest - lastIndexNetworkState);
             end
        end
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);       
        sumInstH = zeros(Op.n,Op.n);
        path = Obs(n,:);
        lpath = size(find(path),2);
        % Compute regular attributes
        for i = 2:lpath - 1
            sumInstU = sumInstU + (Ufull(path(i),path(i+1)) + Vd(min(path(i+1),lastIndexNetworkState + 1)) - Vd(path(i)))/mu(path(i));
            for j = 1:Op.n
                sumInstX(j) = sumInstX(j) + (Atts(j).Value(path(i),path(i+1)) + gradVd(min(path(i+1),lastIndexNetworkState + 1),j) - gradVd(path(i),j))/mu(path(i));
            end
            for x = 1: Op.n
                for y = 1: Op.n
                    sumInstH(x,y) = sumInstH(x,y) + (hesVd(x,y).value(min(path(i+1),lastIndexNetworkState + 1)) - hesVd(x,y).value(path(i)))/ mu(path(i));
                end
            end
        end
        Hs(n).value = sumInstH;
        Hessian = Hessian + (sumInstH - Hessian)/n;
        Gradient(n,:) = Gradient(n,:) + sumInstX;
        lnPn = lnPn + sumInstU ;  
        LL =  LL + (lnPn - LL)/n;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end

function [gradExpV] = getGradExpV(M, Att, op, expV)
    I = speye(size(M));  
    A = I - M; 
    for i = 1:op.n
        u = M .* (Att(i).Value); 
        v = u * expV; 
        gradExpV(i) = Matrix2D(A\v); 
    end
end
