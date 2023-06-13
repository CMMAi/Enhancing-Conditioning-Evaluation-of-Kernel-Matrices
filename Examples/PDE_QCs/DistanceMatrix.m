function DM = DistanceMatrix(dsites,ctrs)
s = 2;
k = 3; N = (2^k+1)^s;
neval = 10; M = neval^s;
[M,s] = size(dsites); [N,s] = size(ctrs);
DM = zeros(M,N);
for d=1:s
[dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);