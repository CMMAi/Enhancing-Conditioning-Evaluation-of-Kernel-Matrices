function [rhs,a,b]=Force(CASE,x,y)
% test functions
switch CASE
    case{1}
        f1 = @(x,y) .75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
        f2 = @(x,y) .75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
        f3 = @(x,y) .5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
        f4 = @(x,y) .2*exp(-((9*x-4).^2+(9*y-7).^2));
        F  = @(x,y) f1(x,y) + f2(x,y) + f3(x,y) - f4(x,y);
        a=0; b=1;
    case{2}
        F = @(x,y) 0.0025./((x-1.01).^2+(y-1.01).^2);
        a=0; b=1;
    case{3}
        F = @(x,y) (tanh(9*y-9*x)+1)/9;
        a=0; b=1;
    case{4}
        F = @(x,y) (1.25+cos(5.4*y))./(6*(1+(3*x-1).^2));
        a=0; b=1;       
    case{5}
        F = @(x,y) 1/3*exp(-81/16*((x-0.5).^2+(y-0.5).^2));
        a=0; b=1;
    case{6}
        F = @ (x,y) 1/3*exp(-81/4*((x-0.5).^2+(y-0.5).^2));
        a=0; b=1;
    case{7}
        F =@(x,y) 1/9*(64-(81*((x-0.5).^2+(y-0.5).^2)))-0.5;
        a=0; b=1;

end
rhs=F(x,y);
end