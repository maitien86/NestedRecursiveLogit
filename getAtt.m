%   Get attribute 
function [] = getAtt()

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;
    global Op;
    global Atts;
    Atts  = objArray(Op.n);   
    Incidence = incidenceFull;
    Atts(1).value = (Incidence .* EstimatedTime);
    Atts(2).value = (Incidence .* TurnAngles);
    Atts(3).value = (Incidence .* LeftTurn);
    Atts(4).value = (Incidence .* Uturn);
end
