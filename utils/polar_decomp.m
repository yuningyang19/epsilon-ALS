function [U] = polar_decomp(V) 
    [P,~,Q] = svd(V,'econ');
    U = P*Q';
    
end