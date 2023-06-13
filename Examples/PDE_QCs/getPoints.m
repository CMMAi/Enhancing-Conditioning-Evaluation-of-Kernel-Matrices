function [bdyp,innerp,center,testp,ka,kb,kc,ki]=getPoints(D,N,plottext)
% uniformely interior points
%Inner and Boundary points
Allpoints0=xlsread('points441.xlsx');
if N==441
    Allpoints=xlsread('points441.xlsx');
else
    gridpoints=linspace(-0.127,0.127,sqrt(N));
    [x,y]=meshgrid(gridpoints);
    Allpoints=[x(:),y(:)];
end

indx=find(Allpoints(:,1)==0.127 | Allpoints(:,1)==-0.127 | ...
          Allpoints(:,2)==0.127 | Allpoints(:,2)==-0.127);
bdyp = Allpoints(indx,:);
innerp = Allpoints(setdiff([1:size(Allpoints,1)],indx),:);
kb=size(bdyp,1); ki=size(innerp,1); ka=kb+ki; kc=kb+ki;
%% Ghost points
t=2*pi*(1:kb)/kb;  
xc=D*cos(t);yc=D*sin(t);
p2=haltonset(2,'Skip',1500); q2=net(p2,5000);
pts2=q2*2*D-D;
in2=inpolygon(pts2(:,1),pts2(:,2),xc,yc);
pts2=pts2(in2,:);
%% center points
center = [pts2(1:kc,1) pts2(1:kc,2)];
%% test points
indxt=find(Allpoints0(:,2)==0 & Allpoints0(:,1)>=0);
testp = Allpoints0(indxt,:);
testp(:,2) = ones(11,1)*0.012;
% testp(1,2) = 0;


if plottext~='plot'
    return
end
% plotting. The rbf data must be specified outside... 
figure()
plot(innerp(:,1),innerp(:,2),'.r',bdyp(:,1),bdyp(:,2),'.g'...
    ,center(:,1),center(:,2),'xk',testp(:,1),testp(:,2),'bo',...
    'MarkerSize',8);axis equal;hold off
legend('interior pts','boundary pts','center pts','evaluation pts')
    set(gca,'FontSize',18);hold off

end