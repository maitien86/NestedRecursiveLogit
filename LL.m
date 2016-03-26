function [f g] = LL(x)
    global Op;
    Op.x = x;
    [f g] = getLL_nested();
    Op.nFev  = Op.nFev + 1;
end