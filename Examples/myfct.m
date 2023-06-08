function fval=myfct(Points, ind)
% the list of unctions for examples
x=Points(:,1);
y=Points(:,2);
switch ind
    case 1
        fval=0.75*exp(-0.25*((9*x-2).^2+(9*y-2).^2)) + ...
            0.75*exp(-1/49*(9*x+1).^2-0.1*(9*y+1).^2)+...
            0.5*exp(-1/4*((9*x-7).^2+(9*y-3).^2))-...
            0.2*exp(-(9*x-4).^2-(9*y-7).^2);
    case 2
        fval=(1/9)*(tanh(9*y-9*x)+1);
    case 3
        fval=(1.25+cos(5.4*y))./(6*(1+(3*x-1).^2));
    case 4
        fval=(1/3)*exp((-81/16)*((x-1/2).^2+(y-1/2).^2));
    case 5
        fval=(1/3)*exp((-81/4)*((x-1/2).^2+(y-1/2).^2));
    case 6
        fval=(1/9)*(64-81*((x-1/2).^2+(y-1/2).^2))-1/2;
    case 7
        fval=exp(-1000*((x-1/2).^6+(y-1/2).^6));
    case 8
        x=-3+6*x;
        y=-3+6*y;
        fval=3*(1-x).^2.*exp(-x.^2-(y+1).^2)-10*(x/5-x.^3-y.^5).*exp(-x.^2-y.^2)-(1/3)*exp(-(x+1).^2-y.^2);
    case 9
        x=-2+4*x;
        y=-2+4*y;
        fval=sin(3*x).*cos(3*y);
    case 10
        fval=x.*exp(-x.^2-y.^2);
    case 11
        fval=(1/9)*(64-81*(abs(x-1/2)+abs(y-1/2)))-1/2;
    case 12
        fval=0.0025./((x-1.01).^2+(y-1.01).^2);
    case 13
        fval=0.05./sqrt(((x-1.01).^2+(y-1.01).^2));
    otherwise
        error('Case not implemented')
end