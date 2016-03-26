% Compute the loglikelohood value and its gradient.
% Based on V -  Relaxing scales mu.
%%
function [LL, grad] = getLL_nested()

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
    global Scale;
    global isFixedMu;
    global LSatt;
    global LinkSize;
    global SampleObs;

    %% If Mu is fixed
    if isFixedMu == true
       [LL, grad] = getLL_nested_fixedMu();
       return;
    end
    
    
    %% Setting sample
    if isempty(SampleObs)
        sample = 1:nbobs;
    else
        sample = SampleObs;
        nbobs = size(find(SampleObs),2);
    end
        
    %% Get M, U
    [lastIndexNetworkState, ~] = size(incidenceFull);     
    LL = 0;
    grad = zeros(1, Op.n);
    mu = getMu(Op.x);
    N = lastIndexNetworkState + 1;
    
    %% Compute gradient of Mu
    gradMu = (zeros(N,Op.n));
    for i = Op.m+1 : Op.n
         gradMu(:,i) = mu .* Scale(:,i - Op.m);
    end 
    %%
    %% Compute MI;    
    M = incidenceFull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));      
    MI = sparse(M); 
    MI(find(M)) = 1;
    
    dPhi = objArray(Op.n - Op.m);
    a = mu;
    k = 1 ./ mu;
    GrdMuCol = objArray(Op.n);
    kM  = spdiags(k,0,N,N); 
    aM =  spdiags(a,0,N,N); 
    for i = 1: Op.n
        GrdMuCol(i).value = spdiags(gradMu(:,i),0,N,N); 
    end
    phi = (kM * MI) .* (MI * aM); 
    e = ones(size(M,1),1); 
    gradPhi = objArray(Op.n);
    for i = 1: Op.n
        if i <= Op.m
           gradPhi(i).value = sparse(MI * 0); 
        else
           gradPhi(i).value =  (MI * spdiags(gradMu(:,i),0,N,N)) .* (kM * MI) - (MI * aM) .*  (spdiags(k .* k .* gradMu(:,i),0,N,N) * MI);
           gradPhi(i).value = sparse(gradPhi(i).value);
        end
    end 
    
    if isLinkSizeInclusive == true
        sizeOfParams = Op.m - 1;
    else
        sizeOfParams = Op.m;
    end 
    AttLc = objArray(Op.n);
    for i = 1 : sizeOfParams
        AttLc(i).value =  sparse(Atts(i).value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    for i = Op.m+1: Op.n
        AttLc(i).value =  sparse(AttLc(1).value * 0);
    end
    gradVd = sparse(zeros(N ,Op.n));
    if isLinkSizeInclusive == false
        prevOD = 0;
        ii = 1;
    else
        prevOD = [0,0];
        ii = 2;
    end
    %% Compute LL and gradient defined over observations
    for n = 1:nbobs
%         fprintf('%d - ',n);
%         if mod(n,30) == 0
%             fprintf('\n');
%         end
        dest = Obs(sample(n), 1);
        if ~isequal(prevOD,Obs(sample(n),1:ii)) ;
            prevOD = Obs(sample(n),1:ii);
            %% set Link Size attributes        
            if isLinkSizeInclusive == true
                LinkSize = LSatt(sample(n)).value;
                Atts(Op.m).value = LinkSize;
                AttLc(Op.m).value =  (LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
                AttLc(Op.m).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
                AttLc(Op.m).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
            end
            %% Compute M matrix
            Mfull = getM(Op.x, kM, isLinkSizeInclusive); % matrix with exp utility for given beta
            M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Mfull(:,dest);
            M(:,lastIndexNetworkState+1) = addColumn;
            M(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            M = sparse(M);
            Ufull = getU(Op.x, isLinkSizeInclusive); % matrix with exp utility for given beta
            U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Ufull(:,dest);
            U(:,lastIndexNetworkState+1) = addColumn;
            U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            U = sparse(U);

            %% get Z 
            MI = sparse(M); 
            MI(find(M)) = 1;
            phi =  (kM * MI) .* (MI * aM);
            phi = sparse(phi);
            [Z, expVokBool] = getZ_NK(M, MI, phi);
            %[Z, expVokBool] = getV_NK(U, M, MI, phi);
            if (expVokBool == 0)
                LL = OptimizeConstant.LL_ERROR_VALUE;
                grad = ones(Op.n,1);
                disp('The parameters not fesible')
                return; 
            end
            %% Compute V
           V = log(Z);
           V = mu .* V; %(bsxfun(@times,mu,V)); 
            %% Compute gradient of V - respect to attributes parameters   
            gradV = objArray(Op.n);    
            for i = 1:Op.n
                % Compute gradient of M
                h = gradMu(:,i) .* log(Z);
                U1 = kM * AttLc(i).value;
                U2 = (kM .* kM .*  GrdMuCol(i).value) * U;
                gradM = sparse(M .* (U1 - U2));
                Zd = MI * spdiags(Z,0,N,N);
                X = MI;
                X(find(MI)) =  Zd(find(MI)) .^ (phi(find(MI)));
                X = spdiags(1 ./ Z,0,N,N) * X;
                H = M .* X;
                U1 = H * spdiags(log(Z),0,N,N);
                U2 = (aM * U1) .* gradPhi(i).value;
                U3 = U1 * GrdMuCol(i).value ;
                S =   aM * (gradM .* X) + U2 - U3;
                gradV(i).value = sparse((speye(size(M)) - H)\( S * e + h));   

            end
        end
        %% Compute gradient of Z    
        lnPn = 0;
        for i = 1: Op.n
             Gradient(n,i) = 0;
             gradVd(:,i) =  gradV(i).value;
        end
        path = Obs(sample(n),:);
        lpath = size(find(path),2);
        sumInstX = zeros(1,Op.n);  
        sumInstU = 0;
        for i = 2:lpath - 1
            midx = min(path(i+1),lastIndexNetworkState + 1);
            sumInstU = sumInstU + (Ufull(path(i),path(i+1)) + V(midx) - V(path(i)))/mu(path(i));
            for j = 1:Op.n
                p1 =   (AttLc(j).value(path(i),midx) + gradVd(midx,j) - gradVd(path(i),j))/mu(path(i));
                p2 =   (Ufull(path(i),path(i+1)) + V(midx) - V(path(i)))/(mu(path(i))^2) * gradMu(path(i),j) ;
                sumInstX(j) = sumInstX(j) + p1 - p2;
            end
        end
        Gradient(n,:) = Gradient(n,:) + sumInstX;
        lnPn = lnPn + sumInstU ;  
        LL =  LL + (lnPn - LL)/n;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
    end
    fprintf('\n');
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end

%%
