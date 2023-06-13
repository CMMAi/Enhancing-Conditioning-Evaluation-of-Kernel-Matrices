
function [v,s] = sigmaN(A,v,Q, R)

% <Purpose>
%   Computing the smallest singular value of an upper triangular matrix.
%
% <Syntax>
%   [v,s] = sigular_n(A);
%
% <Input parameters>
%   1.     A
%
% <Output Parameters>
%   1.     v  -- the smallest singular vector.
%   2.     s  -- the smallest singular value.

n = size(A,2);

% [~,R] = qr(A);

%v = randn(n,1);

tmp = 1/norm(v);v = tmp*v;
u = zeros(1,n);
max_iter = 18;
eps1 = eps;%10^-14; %
for k = 1:n
    if (abs(R(k,k)) < eps1)
        R(k,k) = eps1;
    end
end

for k = 1:max_iter
    %--------------------------- forward substitution
    u(1) = v(1)/R(1,1);
    for kk = 2:n
        u(kk) = (v(kk) - u(1:kk-1)*R(1:kk-1,kk))/R(kk,kk);
    end
    tmp = 1/norm(u);
    u = tmp*u;
    %--------------------------- backward substitution
    v(n) = u(n)/R(n,n);
    for kk = (n-1):-1:1
        v(kk) = (u(kk) - R(kk,kk+1:n)*v(kk+1:n))/R(kk,kk);
    end
    s = 1/norm(v);
    v = s*v;
    %---------------------------
    if (k > 2)
        if ( abs(s - s_old)/s < eps1 )
           break;
        end
        if ( s < eps1 )
           break;
        end
    end
    s_old = s;
end
%-----------------------------------------------------------------------
% End of function
%-----------------------------------------------------------------------
