function [ Uout ] = column_matching( U,Uhat ) 
 % the columns of Uhat might be misplaced, and we need to adjust the place
 % of these columns according to U.
 R = size(U,2);
 
 decision_mat = U'*Uhat;   % e.g., row 1 of decisoin_mat denotes the 
 % closedness of the first column of Uhat to every column of U.
decision_mat = abs(decision_mat); % because ambiguity of the sign, we need only 
% the absolute values of the decision matrix.

if verLessThan('matlab','9.7') % version less than Matlab 2019 a simplified H algorithm is used.
[~,place] = max(decision_mat'); % max takes the maximum value of each column. Thus 
% we need to transpose the decision matrix.
% place records the true place that the columns of Uhat should appear.
% for example, if place= [2 1 3], then it maps [1 2 3] -> [2 1 3];
else 
   [place,~] = matchpairs(decision_mat',-1,'max');  % version above Matlab 2019, use the built-in function.
   place = place(:,1);
end

% put the largest entry of each row of decision_mat to 1, and other to 0,
% then it becoms a permutation matrix. Then Uout = Uhat*decision_mat.
Uout = Uhat(:,place);   % this is exactly the permutation.

negative = diag(Uout'*U); negative = negative <0;       % adjust the sign
Uout(:,negative) = -Uout(:,negative);
end

