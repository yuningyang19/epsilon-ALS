function v =  by_svds(M,options)
% return the i-th column of V by the leading singular value and the left
% singular vector.
    [xi,l,~] = svds(M,1,'L',options);
     v = xi*l;             % the i-th col. is Mi*yi;
end