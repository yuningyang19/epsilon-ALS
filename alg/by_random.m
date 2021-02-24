function v =  by_random(M,~)
% uniformly draw a normalized vector y and return v = M*y;
    y = randn(size(M,2),1); y = y/norm(y);
    v = M*y;
end