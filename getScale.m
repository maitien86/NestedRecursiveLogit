%% Compute ET and LS for the Scales
function [] = getScale()
    global incidenceFull;  
    global Scale;
    global EstimatedTime;
    global Nb_out;
    global Uturn;
    %% Get ET
    maxstates = size(EstimatedTime,1);
    ET = zeros(size(EstimatedTime,2),1);
    I = find(EstimatedTime);
    [nbnonzero] = size(I,1);
    for i = 1:nbnonzero
        [k,a] = ind2sub(size(EstimatedTime), I(i));
        ET(a,1) = EstimatedTime(k,a);     
    end
    ET = ET(1:maxstates+1);
    ET = ET / max(ET);
    ONE = ones(maxstates+1,1); 
    %% Compute flow
    lastIndexNetworkState = size(incidenceFull,1);
    beta = [-1,0,0,0]';
    beta = [-2.5,-1,-0.4,-4]';
    [ok, flow] = getFlow(beta);
    if ok == false
        fprintf(' Beta is wrong, try other parameters \n');
        flow = ones(lastIndexNetworkState,1);
    end
    flow = flow(1: maxstates + 1);
    flow = flow / max(flow);
    flow = flow;
    
    %% Travel time average
    MI = sparse(EstimatedTime); 
    MI(find(EstimatedTime)) = 1;
    AT = sum(EstimatedTime,2) ./ sum(MI,2);
    AT(maxstates+1,1) = 0;
    AT = AT / max(AT);
    %% UTurn average
    %UT = sum(Uturn,2) ./ sum(MI,2);
    %UT(maxstates+1,1) = 0;
    %% Number of outgoing links
    Nb = sum(incidenceFull,2);
    Nb = [Nb;0];
    N = length(Nb);
    Nb_out =  spdiags(Nb,0,N,N);
    %% Set Scales
    Scale = zeros(maxstates + 1,1);
    Scale(:,1) = ET;
    Scale(:,2) = flow;
    Scale(:,3) = Nb;
    %Scale(:,2) = ONE;
end