% This program uses "sigmaN" to calculate the smallest 
%singular value and compares with SVD in terms of CPU time
% this does not go through solutions. For solution see "doTable01.m"
% see instrunction in "doTable01.m"
%%    
clear; 
% close all
global RBFscale
global RBFpar  
global RBFtype 
RBFtype='g' 
RBFpar=2
%% inputs
CASE =1; %Test functions
n = round(logspace(log10(500), log10(10000), 50));
nt = 900; % number of test point
RBFscale = 0.1; % choose a scale parameter
%% get the domain boundary
[~,a,b] =Force(CASE,1,1); 
%%% initialization
timem = zeros(length(n),1);
timep = zeros(length(n),1);
%%%%%%%%%%%%%%%%%%%%%
%% LOOP
disp(sprintf('f%d:',CASE))
%output display
for j=1:length(n)
    disp(sprintf('n=%d',n(j)));
    %get the points
    [coll,cntr,test] = getPoints01(n(j),n(j),nt,a,b); % get points
    %distance matrix
    v = randn(n(j),1); %random points for "sigmaN" usage
    A = kermat(coll,cntr);
    
    %%% get \sigma_n from the proposed method
    tic
    [qq,rr]   = qr(A);
    [v1,smin] = sigmaN(A,v,qq,rr); 
    timep(j)  = toc;
    
    %%% get \sigma_n from MATLAB
    tic
    [U,S,V] = svd(A);SD = diag(S);Smin=SD(end);
    timem(j)= toc;
end

%% plot
plot(n,timem,'r-*','LineWidth',2);hold on
plot(n,timep,'b-o','LineWidth',2)
xlabel('n')
ylabel('CPU(s)')
legend('Matlab alg.','proposed alg.')
set(gca,'FontSize',18);
fprintf('=============================================\n')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%