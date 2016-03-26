%   Get mu
%%
function Mu = getMu(x)   
    global Scale;
    global Op;
    t = Op.n - Op.m;
    y = Op.x(Op.m+1: Op.n);
    Mu = exp(Scale(:,1:t) * y); 
end