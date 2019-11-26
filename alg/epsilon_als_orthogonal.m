function [U,output] = epsilon_als_orthogonal(T,U0,options)
% epsilon ALS for solving tensor CPD where the last t latent factor
% matrices are assumed to be columnwisely orthonormal.
%   [U,output] = epsilon_als_orthogonal(T,U0) computes the factor matrices
%   U{1}, ..., U{d}.
%     The algorithm is initialized with the factor matrices U0{n}.
%   The structure output returns additional information:
%

 
%      output.iterations - The number of iterations.
%      output.relstep    - successive relative difference.
%
%    options:
%      options.e1 and optoins.e2
%                                  - The perturbation parameters. Should be
%                                     nonnegative.
%      options.MaxIter = 2000      - The maximum number of iterations.
%      options.Order = 1:d        - The order in which to update the factor
%                                   matrices. Must be a permutation of 1:d.
%      options.TolX = 1e-4        - The tolerance for output.relstep.

%   Author: Yuning Yang (yyang@gxu.edu.cn)
% Y. Yang, The Epsilon-Alternating Least Squares for Orthogonal Low-Rank
%                   Tensor Approximation and Its Global Convergence,
%                   http://arxiv.org/abs/1911.10921

d = ndims(T);
if d < 3, error('the order should be >= 3.'); end
if ~isfield(options,'t'),  error('there should be at least one factor matrix being columnwise orthogonal.');    end
t = options.t;

if ~isfield(options,'e1'), options.e1 = 10^(-8); end
if ~isfield(options,'e2'), options.e2 = 10^(-8); end
e1 = options.e1;
e2 = options.e2;
% R = size(U0{1},2);

% Check the initial factor matrices U0.
U = U0(:).';
R = size(U{1},2);
if any(cellfun('size',U,2) ~= R)
    error('all factor matrices should have the same column size.');
end
if any(cellfun('size',U,1) ~= size(T))
    error('the row size of the factor matrices should be equal to the size of the data tensor.');
end
if ~isfield(options,'MaxIter'), options.MaxIter = 2000; end
if ~isfield(options,'Order'), options.Order = 1:d; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-5; end
if ~isfield(options,'TolX'), options.TolX = 1e-5; end


M = arrayfun(@(n)tens2mat(T,n).',1:d,'UniformOutput',false);


first = options.Order(1);
last = options.Order(end);
D = kr(U([d:-1:first+1 first-1:-1:1]));

D = M{first}-D*U{first}.';

output.info = false;
output.iterations = 0;
output.relstep = [];

Omega = 0;
while ~output.info
    
    % Save current iterate.
    U1 = U;
    for n = 1: d-t          % Update non-orthogonal factor matrices.
        
        krUn = kr(U([d:-1:n+1 n-1:-1:1])); % KR product of the other factor matrices.
        Tn = M{n}';        % the mode-n unfolding of the tensor T, Mn: n_1\times produ_{j\neq 1}n_j
        Vn = Tn*krUn;      % line 4 of the paper, update the V_j's simultaneously by using KR prod..
        
        if n == 1  % update Omega
            squared_sigma_n = dot(U{n},Vn);      % the <U{n}_i,Vn_i> is exactly sigma_i, if n=1.
            norm_sigma_n = norm(squared_sigma_n); % ||sigma||
            Omega = diag(  squared_sigma_n/norm_sigma_n   ); % Omega = diag( omega_i   ) = diag( sigma_i/||sigma|| )
        end
        
        tilde_Vn = Vn*Omega + e1*U{n};    % line 5 of the paper. Omega = diag(omega_i^k)
        
        col_norm_Vn = sum(abs(tilde_Vn).^2).^(1/2);   % norm of each column of tilde_Vn -> ||tilde_Vn_i||=1
        
        U{n} = tilde_Vn*diag( 1./col_norm_Vn  );        % line 5 of the paper, every column of t_Vn should be normalized.
        
        %      U{n} = tilde_Vn./repmat(col_norm_Vn,size(U{n},1),1);   % the same operation. Which is faster?
        
    end
    for n = d-t+1: d %   Factors related to polar decomp.
        krUn = kr(U([d:-1:n+1 n-1:-1:1])); % KR product of the other factor matrices.
        Tn = M{n}';        % the mode-n unfolding of the tensor T, Mn: n_1\times produ_{j\neq 1}n_j
        Vn = Tn*krUn;      % line 10 of the paper, update the V_j's simultaneously by using KR prod..
        tilde_Vn = Vn*Omega + e2*U{n};    % line 5 of the paper. Omega = diag(omega_i^k)
        U{n} = polar_decomp(tilde_Vn);
    end
    
    output.iterations = output.iterations+1;
    
    output.relstep(end+1) = sqrt(sum(cellfun(@(u,v)(u(:)-v(:))'* ...
        (u(:)-v(:)),U,U1)))/sqrt(sum(cellfun(@(u)u(:)'*u(:),U)));
    if output.relstep(end) <= options.TolX, output.info = 1; end
    if output.iterations >= options.MaxIter, output.info = 2; end
    
end





