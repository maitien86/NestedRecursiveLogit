%   Get MUtility
%%
function Mfull = getStandardM(x, isLS)   

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;    
    global LSatt;
    global isFixedUturn;  
    u1 = x(1) * EstimatedTime;
    u2 = x(2) * TurnAngles;
    u3 = x(3) * LeftTurn;
    if isFixedUturn == false
        u4 = x(4) * Uturn;
        if isLS == true
            u5 = x(5) * LSatt;
        else
            u5 = 0;
        end
    else 
        u4 = -20 * Uturn;
        if isLS == true
            u5 = x(4) * LSatt;
        else
            u5 = 0;
        end
    end
    u = (u1 + u2 + u3 + u4 + u5);
    expM = u ;
    expM(find(incidenceFull)) = exp(u(find(incidenceFull)));
    Mfull = incidenceFull .* expM;
end