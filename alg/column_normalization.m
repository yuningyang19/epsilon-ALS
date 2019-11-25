function A = column_normalization(A)

    col_norm_A  = sum(abs(A).^2).^(1/2);   % norm of each column of Vn -> ||Vn_i||   
    A = A*diag( 1./col_norm_A  );        % line 5 of the paper, every column of t_Vn should be normalized.
end