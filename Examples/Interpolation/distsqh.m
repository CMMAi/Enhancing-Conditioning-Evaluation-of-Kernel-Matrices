function dst=distsqh(p, q)
% calculates a np*nq matrix of halved squares of point distances 
% with two point sets p and q of np and nq points each.
% [np pdim]=size(p);
% [nq qdim]=size(q);
% if pdim~=qdim 
%      error('point sets of unequal dimension')
% end
% if pdim==1
%     dst=(repmat(p,1,nq)-repmat(q',np,1)).^2/2;
%     return
% end
% dst=p*q';
% cp=sum((p.*p)')/2; % squared norms of p points, as row, halved
% cq=sum((q.*q)')/2; % squared norms of q points, as row, halved
% dst=repmat(cp',1,nq)+repmat(cq,np,1)-dst;



%%
[M,s] = size(p); N = length(q); 
DM = zeros(M,N);
for i=1:s
    [dr,cc]=ndgrid(p(:,i),q(:,i));
    DM = DM +(dr-cc).^2;
end
% DM=sqrt(DM);
dst = DM/2;

%%
% dst = pdist2(p,q).^2/2;
