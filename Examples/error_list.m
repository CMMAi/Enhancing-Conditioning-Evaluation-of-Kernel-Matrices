function error_list(type, par, n, nt, CASE)

%%    
global RBFscale
global RBFpar  
global RBFtype 
RBFtype=type
RBFpar=par
%% inputs
c(:,1) = linspace(0.01,1,50);
%% get the domain boundary
[~,a,b] =Force(CASE,1,1); % 
% figure(1),plot(coll(:,1),coll(:,2),'r.',cntr(:,1),cntr(:,2),'ko')
%initialization
error=zeros(length(c),length(n)); 
ECNm=error;ECNp=error;
sp=error; sm=error;
%%%%%%%%%%%%%%%%%%%%%
%% LOOP
disp(sprintf('f%d:',CASE))
for j=1:length(n)
    disp(sprintf('n: %d',n(j)))
    %get the points 
    [coll,cntr,test] = getPoints01(n(j),n(j),nt,a,b); % get points
    %RHS
    [rhs,~,~]=Force(CASE,coll(:,1),coll(:,2)); 
    %exact
    [Exact,~,~]=Force(CASE,test(:,1),test(:,2)); 
    rng(1);
    v = randn(n(j),1); %random points for "sigmaN" usage
    for i=1:length(c)
        RBFscale = c(i);
        %%%get interpolation and evaluation matrices
        A = kermat(coll,cntr);
        Aeval = kermat(test,cntr);
        %%%solve the linear system using SVD
        [U,S,V]  =svd(A); SD = diag(S);
        %%% get the coefficients of "\Alpha=rhs*A^-1"
        cc=U'*rhs; y=S\cc; Alpha=V*y;
        %%%calculate the error
        error(i,j) = norm(Exact-Aeval*Alpha,inf);
        
        %%%calculate the ECN from proposed algorithm
        %%%get the \sigma_n from the porposed meghtod
        [qq,rr] = qr(A);
        [v1,sp(i,j)] = sigmaN(A,v,qq,rr); % sigma_n
        ECNp(i,j) = 1/sp(i,j)*norm(rhs,Inf)/norm(Alpha,Inf); % ecn
        
        %%%calculate the ECN from MATLAB SVD
        %%%get the \sigma_n from Matlab SVD
        sm(i,j) = SD(end); %\sigma_n
        ECNm(i,j) = 1/sm(i,j)*norm(rhs,Inf)/norm(Alpha,Inf); % ecn
        
        %%% check the progress!
        if rem(i,25)==0
            disp(sprintf('scale: %d out of %d are solved',i,length(c)))
        end
        
    end
end
%%
%output display
disp('ready for LaTex')
disp('n      &c(mat.)  &err(mat.)  &&c(prop.)&err(prop.) ');
for j=1:length(n)
    %get maximum ECN resulted by Matlab SVD
    [~,nm] = max(ECNm(:,j)); 
    %get maximum ECN resulted by proposed alg.
    [~,np] = max(ECNp(:,j)); 
    %disp
    disp(sprintf('%d    &%1.2f     &%2.1e     &%1.2f     &%2.1e \\\\'...
                 ,n(j), c(nm), error(nm,j), c(np), error(np,j) ) )
end
%% plot
for j=1:length(n)
    figure()
    yyaxis left
    semilogy(c,error(:,j),'k-','LineWidth',2);
    ax1=gca; ax1.YColor='k';
    yyaxis right
    semilogy(c,ECNm(:,j),'r-x','LineWidth',1,'MarkerSize',6);hold on
    semilogy(c,ECNp(:,j),'m-.o','LineWidth',1,'MarkerSize',8);hold off
    ax2=gca; ax2.YColor='k';
    yyaxis left
    title(sprintf('n= %d',n(j)))
    xlabel('C')
    ylabel('Maxerr')
    yyaxis right
    ylabel('ECN')
    legend('error','mECN','pECN')
    set(gca,'FontSize',18);
    str=sprintf('fig_err_ecn_T%s_n%d_F%d',...
        RBFtype, n(j), CASE);
    saveas(gcf,str,'fig')
%     saveas(gcf,str,'epsc')
    %%% plot the \sigma_n
    figure()
    semilogy(c,sp(:,j),'r-o','LineWidth',2,'MarkerSize',6);hold on
    semilogy(c,sm(:,j),'m-.x','LineWidth',1,'MarkerSize',15);hold off
    % semilogy(c,tol,'k-','LineWidth',2);
    title(sprintf('n= %d',n(j)))
    xlabel('C')
    ylabel('\sigma_{n}')
    legend('proposed alg.','MATLAB alg.')
    set(gca,'FontSize',16);
    str=sprintf('fig_Sn_T%s_n%d_F%d',...
        RBFtype, n(j), CASE);
    ylim([10^-20 10^0])
    yticks([10^-20 10^-15 10^-10 10^-5 10^0])
    saveas(gcf,str,'fig')
%     saveas(gcf,str,'epsc')
  
    
end
fprintf('=============================================\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
