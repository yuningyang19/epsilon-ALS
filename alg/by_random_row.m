function v =  by_random_row(M,~)
% uniformly pick a row and then normalize it to give y. v = M*y;
    idx = randi(size(M,1));
    y = M(idx,:)'; y = y/norm(y);
    v = M*y;
end