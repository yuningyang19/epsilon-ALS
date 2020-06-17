function [ U0 ] = get_initializer( A, t, R )
% generate a good initializer for the epsilon-ALS.
% t the number of columnwise orthonormal factor matrices.


    d = ndims(A);
    sz = size(A);
    if d<t 
        error('the number of ...')
    end
    
    U0 = {};
    
    for j = d-t+1: d
%         B1 = reshape(A,sz(j), prod(sz)/sz(j));
        B1 = tens2mat(A,j,setdiff(1:d,j));
        [Un,~,~]=svd(B1,'econ');
        U0{j} = Un(:,1:R);
    end
    
    uj = {};
    B = {};
    for i = 1: R
        for j = d-t+1: d
          tmp= U0{j};
          uj{j-d+t} = tmp(:,i);          
        end  
        Bi = contract(A,uj,[d-t+1:d]);
        B{i} = Bi/frob(Bi);
        x = rank1approx(Bi);            % x{j} is the j-th column of the i-th factor.
        for j = 1: d-t
            U0{j}(:,i) = x{j};
        end
    end
        
        
 end

