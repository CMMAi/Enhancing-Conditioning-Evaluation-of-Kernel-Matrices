function [IPoints, Cpoints, EPoints, ex, ey]=getPoints01(n,nc,nt,a,b)
%%%%%%%%%%%%%
%% Amir: hand-made function haltonseq to select the halton points
p2 = haltonseq(50*n,2); q2=p2(1000:end,:);
pts=q2*b*2*2-2*b; % go beyond the [0,1]^2
indx = find(pts(:,1)<b & pts(:,1)>a & ...
            pts(:,2)<b & pts(:,2)>a ); % select in domain [0,1]^2
pts2 = pts(indx,:);  
xic=pts2(1:n,1);  yic=pts2(1:n,2); 
%%plot(xic,yic,'bo') % plot to check.
%%%%%%%%%%%%% 
%%%%Finally, collecting the points
IPoints = [xic yic];
Cpoints = [xic(1:nc) yic(1:nc)]; % selecting "nc" center points
%test mesh grid
n1D=sqrt(nt); % the 1 D number of points
h1D=(b-a)/(n1D-1); %
% Evaluation grid RS: changed Sep 30 
[ex, ey]=meshgrid(a+h1D/2:h1D:b-h1D/2,a+h1D/2:h1D:b-h1D/2);
EPoints=[ex(:) ey(:)];
