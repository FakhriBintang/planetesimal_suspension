function [xq,cxq,cmq]  =  equilibrium(T0,c,P,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,PhDg)

perCm = (perCm-cphs0)./(cphs1-cphs0);
perCx = (perCx-cphs0)./(cphs1-cphs0);
perT  = (perT -Tphs0) /(Tphs1-Tphs0);

T   = max(0,min(1,(T0 - P*clap -Tphs0)./(Tphs1-Tphs0)));

a = 15;
b = 0.04;
ind1 = T>=perT+b;
ind2 = T< perT-b;
ind3 = T>=perT-b & T<perT+b;

cx1 = max(0,min(1,          perCx .*erfc(PhDg(1).*(T-perT)./(1-perT))));
cx2 = max(0,min(1, perCx+(1-perCx).*erfc(PhDg(2).*(T-   0)./   perT) ));

dcdT = -2*perCx*PhDg(1)/sqrt(pi)./(1-perT);
cx1(T<perT) = perCx + dcdT.*(T(T<perT)-perT);

cxq = zeros(size(T));
if perT > 0
    cxq(ind1) =  cx1(ind1);
    cxq(ind2) =  cx2(ind2);
    cxq(ind3) = (cx1(ind3).^-a+cx2(ind3).^-a).^-(1/a);
else
    cxq = cx1;
end

cm1 = max(0,min(1,          perCm .*erf(PhDg(3).*(1   -T)./(1-perT))./erf(PhDg(3))));
cm2 = max(0,min(1, perCm+(1-perCm).*erf(PhDg(4).*(perT-T)./(  perT))./erf(PhDg(4))));

cmq = zeros(size(T));
if perT > 0
    cmq(ind1) =  cm1(ind1);
    cmq(ind2) =  cm2(ind2);
    cmq(ind3) = (cm1(ind3).^a+cm2(ind3).^a).^(1/a);
else
    cmq = cm1;
end

cxq = max(cphs0,min(c,cphs0 + cxq.*(cphs1-cphs0)));
cmq = min(cphs1,max(c,cphs0 + cmq.*(cphs1-cphs0)));

xq  = max(0,min(1, (c-cmq)./(cxq-cmq) ));


