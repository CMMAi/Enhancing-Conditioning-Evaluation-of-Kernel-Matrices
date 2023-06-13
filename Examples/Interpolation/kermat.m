function mat=kermat(X, Y)
% creates kernel matrix
% for two point sets X and Y  
global RBFscale
% [nx dx]=size(X);
% [ny dy]=size(Y);
% if dx~=dy
%     error('Unequal space dimension for kermat arguments');
% end
DM = distsqh(X,Y);
% [row1, col1] = find(isnan(DM));

mat=frbf(DM/RBFscale^2,0);
% [row2, col2] = find(isnan(mat));

