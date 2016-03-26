%   Get analytical Hessian matrix
%   Link size is included
%%
function Hessian = getHessian()
   [~,~,Hessian,~] = getAnalyticalHessian_nested();    
end

