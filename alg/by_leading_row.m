function v =  by_leading_row(M,~)
% return the i-th column of V by M*y, where y is the row of M with largest
% magnitude (y is also normalized).
    [l,idx] = max( sqrt(sum(M.^2,2)) );
    y = M(idx,:)'/l;
    v = M*y;
end