global Op;
step = 0.000001;
x = Op.x;
[val1, gr1] = LL(Op.x);
gr1
H  = eye(Op.n);
G = [];
for i=1:Op.n
    h = step * H(:,i);
    Op.x = x +h;
    [val2, gr2] = LL(Op.x);
    G(i) = (val2 - val1)/step;
end
G'