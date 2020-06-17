function [ norm_diff ] = diff_U_Uout( U, Uout ) 
%DIFF_U_UOUT Summary of this function goes here
%   Detailed explanation goes here
norm_diff = sqrt(sum(cellfun(@(u,v)(u(:)-v(:))'* ...
        (u(:)-v(:)),U,Uout)))/sqrt(sum(cellfun(@(u)u(:)'*u(:),U)));

end

