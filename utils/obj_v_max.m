function [ v ] = obj_v_max( A,U ) 
%OBJ_V_MAX Summary of this function goes here
%   Detailed explanation goes here

    v = 0;
    R = size(U{1},2); 
    d = ndims(A);
    for i = 1: R
        u1={};
        for j = 1: d
            u1{j} = U{j}(:,i);
        end
        v = v + contract(A,u1,1:d)^2;
    end

end

