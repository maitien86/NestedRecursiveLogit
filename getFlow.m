%   Compute Link Size attribute from data
%   
%%
function [ok, Flow] = getFlow(beta)
    global incidenceFull; 
    global Obs;     % Observation
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;    
    % ----------------------------------------------------
    mu = 1; % MU IS NORMALIZED TO ONE
    Obs_temp = Obs;
    Obs_temp(find(Obs)) =  Obs(find(Obs)) + 1;
    % Save data
    IF = (incidenceFull);
    ET = (EstimatedTime);
    TA = (TurnAngles);
    LT = (LeftTurn);
    UT = (Uturn);
    
    incidenceFull = transferMarix(incidenceFull);
    EstimatedTime = transferMarix(EstimatedTime);
    TurnAngles = transferMarix(TurnAngles);
    LeftTurn = transferMarix(LeftTurn);
    Uturn = transferMarix(Uturn);    
    lastIndexNetworkState = size(incidenceFull,1);
    dest = unique(Obs_temp(:,1));
    orig = unique(Obs_temp(:,2));
    incidenceFull(1,orig') = 1;
    EstimatedTime(1,orig') = 0.05;
    incidenceFull(dest,lastIndexNetworkState) = 1;
    EstimatedTime(dest,lastIndexNetworkState) = 0.05;
    
    M = getStandardM(beta,false);
    [expV, expVokBool] = getExpV(M);
    if (expVokBool == 0)
       ok = false;
       Flow = [];
       disp('ExpV is not fesible')
       return; 
    end  
    P = getP(expV, M);
    G = sparse(zeros(size(expV)));
    G(1) = 1;
    I = speye(size(P));
    F = (I-P')\G;                        
    if (min(F) < 0)                
        ToZero = find(F <= 0);
        for i=1:size(ToZero,1)
            F(ToZero(i)) = 1e-9;
        end
    end
    ips = F;
    ips = ips(2:lastIndexNetworkState - 1 ); 
    Flow = ips;
    % Recover real size
    incidenceFull = IF;
    EstimatedTime = ET;
    TurnAngles = TA;
    LeftTurn = LT;
    Uturn = UT;
    ok = true;
end