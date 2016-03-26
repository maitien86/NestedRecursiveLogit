% Compute the loglikelohood value and its gradient.
% Based on V -  Relaxing scales mu.
%%
function [LL, grad] = getPLL_nested(sample)

    global incidenceFull; 
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global isLinkSizeInclusive;
    global lastIndexNetworkState;
    global mu;
    global LSatt;
    global LinkSize;

    %% If Mu is fixed
    %% Get M, U
    lastIndexNetworkState = size(incidenceFull,1);     
    LL = 0;
    grad = zeros(1, Op.n);
    mu = getMu(Op.x);
    N = lastIndexNetworkState + 1;
    %% Compute phi(a|k)=mu_a / mu_k;    
    Mfull = getM(Op.x, isLinkSizeInclusive); % matrix with exp utility for given beta
    M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));      
    %% 
    if isLinkSizeInclusive == true
        sizeOfParams = Op.m - 1;
    else
        sizeOfParams = Op.m;
    end 
    AttLc = objArray(Op.n);
    for i = 1 : sizeOfParams
        AttLc(i).value =  (Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    for i = Op.m+1: Op.n
        AttLc(i).value =  sparse(zeros(size(AttLc(1).value)));
    end    
    %% Loop over observations
    fprintf(' n = ');
    nobs = size(find(sample),2);
    for t = 1:nobs
        fprintf('%d - ',t);
        if mod(t,20) == 0
            fprintf('\n');
        end
        dest = Obs(sample(t), 1);       
        %% set Link Size attributes
        LinkSize = LSatt(sample(t)).value;
        if isLinkSizeInclusive == true  
            Atts(Op.m).Value = LinkSize;
            AttLc(Op.m).value =  (LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
            AttLc(Op.m).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
            AttLc(Op.m).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
        end
        %% Compute M matrix
        Mfull = getM(Op.x, isLinkSizeInclusive); % matrix with exp utility for given beta
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
       [Z, expVokBool] = getZ(M);
       if (expVokBool == 0)
            LL = OptimizeConstant.LL_ERROR_VALUE;
            grad = ones(Op.n,1);
            disp('The parameters not fesible')
            return; 
       end
       %% Compute V
       V = log(Z);
       V = (bsxfun(@times,mu,V)); 
        %% Compute gradient of Z    
        lnPn = 0;
        sumInstU = 0;
        path = Obs(sample(t),:);
        lpath = size(find(path),2);
        % Compute regular attributes
        for i = 2:lpath - 1
            midx = min(path(i+1),lastIndexNetworkState + 1);
            sumInstU = sumInstU + (Ufull(path(i),path(i+1)) + V(midx) - V(path(i)))/mu(path(i));
        end
        lnPn = lnPn + sumInstU ;  
        LL =  LL + lnPn;
    end
    fprintf('\n');
end

%%
