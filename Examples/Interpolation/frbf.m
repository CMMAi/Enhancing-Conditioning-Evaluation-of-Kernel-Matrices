function y =  frbf (s, k)   
% k-th derivative of standard rbf kernel in f form, i.e.
% as a function of s=r^2/2.
global RBFtype;   % 'g'=gaussian, 'mq' = multiquadric
global RBFpar;
switch lower(RBFtype)
    case ('g')   % 'g'=gaussian,
        par=RBFpar;
        if mod(k,2)==0
            y=exp(-s);
        else
            y=-exp(-s);
        end
    case ('mq') % 'mq' = multiquadric, inverse or not...
        fac=1;
        ord=k;
        par=RBFpar;
        while ord>0
            ord=ord-1;
            fac=fac*par;
            par=par-2;
        end
        y=fac*(1+2*s).^(par/2);    
    case ('p')  % powers 
        fac=1;
        ord=k;
        par=RBFpar;
        while ord>0
            ord=ord-1;
            fac=fac*par;
            par=par-2;
        end
        y=fac*(2*s+eps).^(par/2);  
    case ('tp')  % thin-plate 
        fac=1;
        ord=k;
        par=RBFpar;
        su=0;
        while ord>0
            ord=ord-1;
            if ord==k-1
                su=1;
            else
                su=su*par+fac;
            end           
            fac=fac*par;
            par=par-2;
        end
        y=(2*s+eps).^(par/2);
        y=fac*y.*log(2*s+eps)/2 +su*y;    
    case ('ms')  % Matern/Sobolev
        y=(-1)^k*besselk(RBFpar-k,sqrt(2*s+eps)).*(sqrt(2*s+eps)).^(RBFpar-k);
    case ('w') % Wendland functions.
        % we use only those which are pos. def. in dimension at most 3.
        [coeff, expon]=wendcoeff(3+2*k, RBFpar-k);
        r=sqrt(2*s);
        ind=find(r<=1);
        u=zeros(size(ind));
        sp=ones(size(ind));
        sloc=r(ind);
        for i=1:length(coeff)
            u=u+coeff(i)*sp;
            sp=sp.*sloc;   
        end
        u=u.*(1-sloc).^expon;
        y=zeros(size(s));
        y(ind)=(-1)^k*u;
    otherwise
        error('RBF type not implemented')
end