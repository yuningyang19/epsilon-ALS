% rng('default')clc; clear;addpath algaddpath utilsaddpath Tensorlabn=10;dim = [	20  20 20   ]; sz = dim;r = 10;nl = 0.1;t=1;MaxIter = 2000;%%%%%% settingomega = abs(randn(r,1));Omega = diag(omega);A01 = rand(dim(1),r)-1/2; A01 = column_normalization(A01); A1= A01*Omega;B1 = rand(dim(2),r)-1/2; B1 = column_normalization(B1);C1 = rand(dim(3),r)-1/2;[P,L,Q] = svd(C1,'econ');C1 = P; Utrue = {A01,B1,C1 };dim = [size(A1,1) size(B1,1) size(C1,1)  ]; U1={A1,B1,C1 };T1 = cpdgen(U1);N = randn(dim);T = T1/frob(T1) + nl* N/frob(N);%%%% end of the setting of the problem%%%% beginning of different epsilon.options.TolX = 1e-4;options.MaxIter = MaxIter;options.t = t;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% good initializer% can also use get_initializer2, get_initializer3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tic[ U0 ] = get_initializer( T, t, r );t_init = toc;g_v_good_init = obj_v_max(T,U0);options.e1 = 0.000; options.e2 = 0.000;tic[U,out] = epsilon_als_orthogonal(T,U0,options);t01=toc;iter01 = out.iterations; U00{1} = column_matching(A01,U{1});U00{2} = column_matching(B1,U{2});U00{3} = column_matching(C1,U{3}); norm_diff01   = diff_U_Uout( Utrue, U00  );gval01 =  obj_v_max( T,U );%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%options.e1 = 0.00000001; options.e2 = 0.00000001;tic[U,out] = epsilon_als_orthogonal(T,U0,options);t11=toc;iter11 = out.iterations;U1{1} = column_matching(A01,U{1});U1{2} = column_matching(B1,U{2});U1{3} = column_matching(C1,U{3}); norm_diff11   = diff_U_Uout( Utrue, U1  );%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%options.e1 = 0.000001; options.e2 = 0.000001;tic[U,out] = epsilon_als_orthogonal(T,U0,options);t21=toc;iter21 = out.iterations;U2{1} = column_matching(A01,U{1});U2{2} = column_matching(B1,U{2});U2{3} = column_matching(C1,U{3}); norm_diff21   = diff_U_Uout( Utrue, U2  );%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%options.e1 = 0.0001; options.e2 = 0.0001;tic[U,out] = epsilon_als_orthogonal(T,U0,options);t31=toc;iter31 = out.iterations;U{1} = column_matching(A01,U{1});U{2} = column_matching(B1,U{2});U{3} = column_matching(C1,U{3}); norm_diff31   = diff_U_Uout( Utrue, U  );fprintf('iter: e0 %d, e10^-8 %d, e10^-6 %d, e10^-4 %d \n', iter01,iter11,iter21,iter31 )fprintf('rel-err: e0 %.2f, e10^-8 %.2f, e10^-6 %.2f, e10^-4 %.2f \n', norm_diff01,norm_diff11,norm_diff21,norm_diff31)fprintf('time: e0 %.2f, e10^-8 %.2f, e10^-6 %.2f, e10^-4 %.2f \n', t01,t11,t21,t31)