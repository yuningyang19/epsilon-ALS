function [ x ] = rank1approx( B )
% B is a d-th order tensor, x is a tuple of vectors that approximates the
% leading singular vectors of B with a theoretical lower bound.

d = ndims(B); sz = size(B);
x = {};
if d == 1
    x{1} = B/norm(B(:)); 
    if size(B,1) < size(B,2)
        x{1}= x{1}';
    end
elseif d == 2
    if size(B,1) == 1
        x{1} = B/norm(B(:));
    elseif size(B,2) == 1
        x{1} = B'/norm(B(:));
    else
        [x{1},l,x{2}] = svds(B,1,'L');
    end
else
%     mB_col_sz = prod(sz)/(sz(end)*sz(end-1));
%     mB_row_sz = sz(end)*sz(end-1);
%     mB = reshape(B,mB_col_sz,mB_row_sz);
    mB = tens2mat(B,1:d-2,d-1:d);
    [~,l,x_mminus1_m] = svds(mB,1,'L');
    Xx_mminus1_m = reshape(x_mminus1_m,sz(end-1),sz(end)); %step 4
    tuple_x_mminus1_m = rank1approx(Xx_mminus1_m);
    x_mminus_1 = tuple_x_mminus1_m{1}; x_m = tuple_x_mminus1_m{2};
    
    tensor_X_1_mminus2 = contract(B,tuple_x_mminus1_m,[d-1 d]);     %step 5
    tuple_x_1_mminus2 = rank1approx(tensor_X_1_mminus2);
    
    
    x = tuple_x_1_mminus2; lenx = length(x); 
    x{lenx+1} = x_mminus_1; x{lenx+2} = x_m;    
end
        


end

