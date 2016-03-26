%   Get MUtility
%%
function Mfull = getM(x,kM,isLS)   

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;    
    global LinkSize;
    global isFixedUturn;  
    global Nb_out;
    u1 = x(1) * EstimatedTime;
    u2 = x(2) * TurnAngles;
    u3 = x(3) * LeftTurn;
    if isFixedUturn == false
        u4 = x(4) * Uturn;
        if isLS == true
            u5 = x(5) * LinkSize;
        else
            u5 = 0 * LeftTurn;
        end
    else 
        u4 = -20 * Uturn;
        if isLS == true
            u5 = x(4) * LinkSize;
        else
            u5 = 0 * LeftTurn;
        end
    end
    u = sparse(u1 + u2 + u3 + u4 + u5);
    sizeU = size(u,1);
    kM = kM(1:sizeU,1:sizeU); 
    u = kM * u;
    
    expM = u ;
    expM(find(incidenceFull)) = exp(u(find(incidenceFull)));
    Mfull = incidenceFull .* expM;

end