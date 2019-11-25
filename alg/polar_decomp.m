function [U] = polar_decomp(V)
    [P,Lambda,Q] = svd(V,'econ');
    U = P*Q';
    
end