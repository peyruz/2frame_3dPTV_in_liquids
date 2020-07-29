function [c1Ref,c2Ref,dRF]=selfCalibOpt7(c1Est, c2Est,sz01,sz02,xi1r, xi2r,PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y,n,e3)

xini=zeros(1,7); % initial guess
    
g = @(x)selfCalibOptObjFun7(x,c1Est,c2Est,sz01,sz02,xi1r, xi2r,PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y,e3,n);
options = optimset('Display','iter','MaxFunEvals',60000,'MaxIter',60000,'TolFun',1e-8,'TolX',1e-8);
[x, ~] = fminsearch(g,xini,options);


c1Ref=c1Est+[x(1); x(2); x(3)+x(7)];
c2Ref=c2Est+[x(4); x(5); x(6)+x(7)];
dRF=x(7);

end

