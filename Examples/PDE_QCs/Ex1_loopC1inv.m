%Edited
clear
close all
% static loading
% clamped B.C.
% solving the main equation with Gauss Jordan
% MQ RBF coupled with ghost and traditinal center points
% regular points
%% operators NMQ RBF
%%%% scale shape: r/c
rbf    = @(c,r)        sqrt(1+(r./c).^2);
dxrbf  = @(c,r,dx)     dx./(c.^2*(r.^2/c.^2 + 1).^(1/2));
dyrbf  = @(c,r,dy)     dy./(c.^2*(r.^2/c.^2 + 1).^(1/2));
dxxrbf = @(c,r,dx)     (c.^2*(r.^2./c.^2 + 1) - dx.^2)...
                        ./(c.^4*(r.^2./c.^2 + 1).^(3/2));
dyyrbf = @(c,r,dy)     (c.^2.*(r.^2./c.^2 + 1) - dy.^2)...
                        ./(c.^4*(r.^2./c.^2 + 1).^(3/2));
dxyrbf = @(c,r,dx,dy) -(dx.*dy)./(c.^4*(r.^2/c.^2 + 1).^(3/2));

%% Defining the parametrs
H=0.012;  a=0.254;
kappa=5.0/6.0;
rho=4180;
K1=12.2e10; K2=2.4e10;
c11 = 23.43e10;  c22 = 23.43e10;  c12 = 5.741e10;
c44 = 7.019e10;  c55 = 7.019e10;  c66 = (c11 - c12) / 2.0; 
Me = c66;
R = Me*0.4; %%% change this either "Me*0" or "Me*0.4"
R5=R; R6=R; R1=0.0; R2=0.0;
d1=c11*H^3/12; d2=c22*H^3/12; d3=(c12+c66)*H^3/12; d4=kappa*H*c44;
d5=kappa*H*R5; d6=c66*H^3/12; d7=kappa*H*c55; d8=kappa*H*R6;
%% import MLPG results
[MLPG_psi1,MLPG_psi2,MLPG_u30,MLPG_W3] = comparison(R,Me);
MLPGresultsList = [MLPG_psi1,MLPG_psi2,MLPG_u30,MLPG_W3];
%% Points
D = 0.4; % the ghost circle radius
N = 15^2; % number of total points
[bdyp,innerp,center,testp,ka,kb,kc,ki]=getPoints(D,N,'plot'); %ghost 
% center=[innerp;bdyp]; % Kansa 
%% DifferenceMatrix
dx = DifferenceMatrix(innerp(:,1),center(:,1));
dy = DifferenceMatrix(innerp(:,2),center(:,2));
%% Force matrix
q=-1e7; %loading
Z1=zeros(ka,1);
Z2=zeros(ka,1);
Z3_bdy=zeros(kb,1);
Z3_inn=q*ones(ki,1);
Z4=zeros(ka,1);
force=[Z1;Z2;Z3_bdy;Z3_inn;Z4];
%% Distance matrix & Shape parameter
DM1=DistanceMatrix(bdyp,center);
DM2=DistanceMatrix(innerp,center);
DM3=DistanceMatrix(testp,center);
%% loop on the shape parameter
Clist  = linspace(0.0001,0.5,100);
Lc       = length(Clist);
rng()
v = randn(4*N,1); %random points for "sigmaN" usage
for i = 1:Lc
    c = Clist(i);
    %%%derivatives
    RBF  = rbf(c,DM1);%RBF for the BC
    RBFi = rbf(c,DM2); %RBF for interior
    %operators for the interior 
    S1   = dxrbf(c,DM2,dx);         S2   = dyrbf(c,DM2,dy);
    S11  = dxxrbf(c,DM2,dx);        S22  = dyyrbf(c,DM2,dy);
    S12  = dxyrbf(c,DM2,dx,dy);
    %%% A matrix: 1st governing eq.
    A11 =  d1*S11+d6*S22-d7*RBFi;   A12 =  d3*S12;
    A13 = -d7*S1;                   A14 = -d8*S1;
    %2nd
    A21 =  d3*S12;                  A22 =  d6*S11+d2*S22-d4*RBFi;
    A23 = -d4*S2;                   A24 = -d5*S2;
    %3rd
    A31 = d7*S1;                    A32 = d4*S2;
    A33 = d7*S11+d4*S22;            A34 = d8*S11+d5*S22;
    %4th
    A41 = (R1+R6)*S1;               A42 = (R2+R5)*S2;
    A43 = R6*S11+R5*S22;            A44 = K1*S11+K2*S22;
    %zero matrix
    Z   = zeros(kb,ka);
    
    % construct the kernel matrix
    A   = [RBF Z   Z   Z     %1st BC
           A11 A12 A13 A14   %1st govening equation
        
           Z   RBF Z   Z     %2nd BC
           A21 A22 A23 A24   %2nd govening equation
        
           Z   Z   RBF Z     %3rd BC
           A31 A32 A33 A34   %3rd govening equation
        
           Z   Z   Z   RBF   %4th BC
           A41 A42 A43 A44]; %4th govening equation
    
    %%%get \sigma_n from Matlab SVD
%     tic
%     [U,S,V]      = svd(A);
%     SD           = diag(S);
%     smin      = SD(end);
%     ti(i) = toc; % record the time for \sigma_n calcuation
    
    %%%get \sigma_n from the porposed
    tic
    [qq,rr]   = qr(A);
    [v1,smin] = sigmaN(A,v,qq,rr);
    ti(i) = toc; % record the time for \sigma_n calcuation
    
    %%solvers
    x   = gauss_jordan(A,force);
    
    %%calculate the ECN 
    kappaEffe(i) = norm(force)/(norm(x)*smin);
    
    %Result
    GRBFg_psi1 =  rbf(c,DM3)*x(1:ka,1);
    GRBFg_psi2 =  rbf(c,DM3)*x(ka+1:2*ka,1);
    GRBFg_u30  =  1000*rbf(c,DM3)*x(2*ka+1:3*ka,1);
    GRBFg_w3   = -1000*rbf(c,DM3)*x(3*ka+1:4*ka,1);
    
    %make a list
    resultlist(:,i) = [GRBFg_psi1;GRBFg_psi2;GRBFg_u30;GRBFg_w3];
    
    %max abs. error difference between RBF and MLPG
    errPsi1 = norm((GRBFg_psi1-MLPG_psi1),inf);
    errPsi2 = norm((GRBFg_psi2-MLPG_psi2),inf);
    errU30  = norm((GRBFg_u30-MLPG_u30),inf);
    errW3   = norm((GRBFg_w3-MLPG_W3),inf);
    %make a list dim: (4variables X 100scales)
    errlist(:,i) = [errPsi1;errPsi2;errU30;errW3];
end
%% display the avarage time per scale for getting \sigma_n
fprintf('CPU time/c  = %3.2f\n', mean(ti))
fprintf('-----------------------\n')
%% plots
%%%%%plot error vs ECN
var = ["\Psi_1","\Psi_2","u_{30} x 10^{-3}","-W_3 x 10^{-3}"];
if length(Clist) > 1 % no plot if only one "c" implemented.
for ii=1:length(var)
    figure()
    yyaxis left
    semilogy(Clist,errlist(ii,:),'b-','LineWidth',2);
    ax1=gca; ax1.YColor='b';
    yyaxis right
    semilogy(Clist,kappaEffe,'r-','LineWidth',2);hold on
%     semilogy(scalist,loocv,'g-','LineWidth',2)
    ax2=gca; ax2.YColor='r';
    yyaxis left
    xlabel('shape parameter')
    ylabel('Maxerr')
    yyaxis right
    ylabel('Eff. cond. no.')
    legend('error','ECN')
    title(sprintf('Results for %s',var(ii)))
    set(gca,'FontSize',18);hold off
end
end
%%
%%%%%plot RBF solution vs MLPG
%plot the results on "11" evaluation points.
%% results
[maxkappaEffe,n1] = max(kappaEffe); 
% "n1" is chosen scale parameter index 
%because it leads to largest ECN
% plot the results corresponds to "n1"

for ii=1:length(var)
    figure()
    plot(testp(:,1)./a,MLPGresultsList(:,ii),'ko','MarkerSize',8,'LineWidth',2);hold on
    plot(testp(:,1)./a,resultlist((ii-1)*11+1:11*ii,n1),'r-','LineWidth',2);
    xlabel('x/a')
    ylabel(var(ii))
    legend('MLPG','GRBF')
%     title(sprintf('Results for %s',var(ii)))
    set(gca,'FontSize',18);hold off
end
