
% get 1/mu * Grad(V0)
%%
function [gradV0] = getGradV0(M, Att, op, expV, origin)
    gradExpV = getGradExpV(M, Att, op, expV);
    gradV0 = zeros(1,op.n);
    for i = 1:op.n
        gradV0(i) =  gradExpV(i).Value(origin)/expV(origin);
    end    
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













































































