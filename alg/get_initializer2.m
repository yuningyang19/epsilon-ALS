function [ U0 ] = get_initializer2( A, t, R )
% generate a good initializer for the epsilon-ALS.
% t the number of columnwise orthonormal factor matrices.


    d = ndims(A);
    sz = size(A);
    if d<t 
        error('the number of ...')
    end
    
    U0 = {};
    B= A;
    for j = d: -1: d-t+1
%     the order of the tensor should be j.
        B1 = tens2mat(B,j,setdiff(1:j,j));
        [Un,~,~]=svd(B1);
        U0{j} = Un(:,1:R);
        B0 = tmprod(B,U0{j}',j);
        B = tmprod(B0,ones(1,R),j); B = squeeze(B); % reduce the order by 1.
    end
    % after the above computation, the order of B0 should be d-t+1.
%     uj = {};
    
     for i = 1: R
         ei = zeros(1,R); ei(i) = 1;
        Bi = tmprod(B0,ei,d-t+1);
%         B{i} = Bi/frob(Bi);
        x = rank1approx(Bi);            % x{j} is the j-th column of the i-th factor.
        for j = 1: d-t
            U0{j}(:,i) = x{j};
        end
    end
        
        
 end

