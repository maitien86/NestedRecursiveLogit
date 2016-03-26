global Op;
Op.x = [-1;-1;-1;-2];
[~, grad1, Hessian1, ~] = getAnalyticalHessian_Vnested();
Hessian1 = -Hessian1;
eps = 1e-7; 
h = eps * [0 0 0 1]';
Op.x = Op.x +h;
[~, grad2] = getLL_nested()
Hesi = (grad2 - grad1)/eps;
Hessian1
Hesi