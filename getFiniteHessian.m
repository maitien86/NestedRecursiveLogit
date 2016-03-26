%   Get finite Hessian matrix
%   Link size is included
%%
function Hessian = getFiniteHessian()
    global Op;
    x = Op.x;
    step = 1e-7;
    H = eye(Op.n) * step;
    Hessian = zeros(Op.n);
    [~, grad1] = getLL_nested();
    for i= 1 : Op.n
        Op.x = x + H(:,i);
        [~, grad2] = getLL_nested();
        Hessian(:,i) = (grad2 - grad1)/step;
    end
    Op.x = x;
end