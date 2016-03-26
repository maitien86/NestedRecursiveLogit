% %% Test nested logit model :
% A =[-2,-3,-4,-4,-3.5,-3];
% P = zeros(6,1);
% %  Nested logit parameter ...
% m1 = 0.80;
% m2 = 0.50;
% 
% expA = exp(A);
% g1 = [1,2,3];
% g2 = [4,5,6];
% V(1) = m1 * log(sum(exp(A(g1)/m1)));
% V(2) = m2 * log(sum(exp(A(g2)/m2)));
% expV = exp(V);
% for i=1:3
%     P(i) = expV(1)/sum(expV) * exp(A(i)/m1)/ sum(sum(exp(A(g1)/m1)));
% end
% for i=4:6
%     P(i) = expV(2)/sum(expV) * exp(A(i)/m2)/ sum(sum(exp(A(g2)/m2)));
% end
% %sum(P(1:3))
% P(2)/P(3)
% P
% sum(P)




%% Test cross nested logit model :
% A =[-2,-3,-4,-4,-3.5,-3];
% P = zeros(6,1);
% %  Nested logit parameter ...
% m1 = 0.8;
% m2 = 0.5;
% a31 = 0.99999;
% a41 = 0.40;
% a32 = 1 - a31;
% a42 = 1 - a41;
% expA = exp(A);
% g1 = [1,2,3,4];
% g2 = [3,4,5,6];
% 
% B = A;
% B(3) = B(3) + log(a31);
% B(4) = B(4) + log(a41);
% C = A;
% C(3) = C(3) + log(a32);
% C(4) = C(4) + log(a42);
% 
% V(1) = m1 * log(sum(exp(B(g1)/m1)));
% V(2) = m2 * log(sum(exp(C(g2)/m2)));
% expV = exp(V);
% for i=1:2
%     P(i) = expV(1)/sum(expV) * exp(B(i)/m1)/ sum((exp(B(g1)/m1)));
% end
% for i=5:6
%     P(i) = expV(2)/sum(expV) * exp(C(i)/m2)/ sum((exp(C(g2)/m2)));
% end
% P(3)  = expV(1)/sum(expV) * exp(B(3)/m1)/ sum((exp(B(g1)/m1))) +  expV(2)/sum(expV) * exp(C(3)/m2)/ sum((exp(C(g2)/m2)));
% P(4)  = expV(1)/sum(expV) * exp(B(4)/m1)/ sum((exp(B(g1)/m1))) +  expV(2)/sum(expV) * exp(C(4)/m2)/ sum((exp(C(g2)/m2)));
% %sum(P(1:3))
% P
% sum(P)
% 

%% Test Link-nested logit model
N = zeros(9,6);
Path =[-2,-3,-4,-4,-3.5,-3];
link = -[1,1,2,2,2,1,1.5,1,1];
nest = [1,1,1,2,3,4,5,5,5,6,7,8,9,9];
alts = [1,2,3,1,2,3,4,5,6,4,5,6,3,4];
for i = 1:size(nest,2)
    N(nest(i),alts(i)) = link(nest(i))/Path(alts(i));
end
mu = 0.6;
V = zeros(1,9); 
% compute value function
for i = 1:9
    nestP = find(N(i,:));
    for j = 1: size(nestP,2)
        r = nestP(j);
        V(i) = V(i) + exp( (Path(r) + log(N(i,r)))/mu);
    end
    V(i) = mu * log (V(i));
end
Pro_nest = exp(V)/ sum(exp(V));
Prob = zeros(1,6);
for a = 1:6    
    for n = 1:9
        if N(n,a)> 0
           Prob(a) = Prob(a) + Pro_nest(n) * exp( (Path(a) + log(N(n,a)))/mu)/ exp(V(n)/mu);  
        end
    end
    
end 
Prob'

















