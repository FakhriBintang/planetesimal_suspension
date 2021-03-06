function [xq,cxq,cmq]  =  equilibrium(T0,c,P,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,PhDg,TINY)

T   = max(0,min(1,(T0 - P*clap -Tphs0)./(Tphs1-Tphs0))) ;
    perCm = (perCm-cphs0)./(cphs1-cphs0);
    perCx = (perCx-cphs0)./(cphs1-cphs0);
    perT  = (perT -Tphs0)./(Tphs1-Tphs0);
    
    cx1 = max(TINY,min(1-TINY,          perCx .*erfc(PhDg(1).*(T-perT)./(1-perT))));
    cx2 = max(TINY,min(1-TINY, perCx+(1-perCx).*erfc(PhDg(2).*(T     )./   perT) ));
    
    cxq = zeros(size(T));
    cxq(T>=perT) = cx1(T>=perT);
    cxq(T< perT) = cx2(T< perT);
    
    cm1 = max(TINY,min(1-TINY,          perCm .*erf(PhDg(3).*(1   -T)./(1-perT))./erf(PhDg(3))));
    cm2 = max(TINY,min(1-TINY, perCm+(1-perCm).*erf(PhDg(4).*(perT-T)./(  perT))./erf(PhDg(4))));
    
    cmq = zeros(size(T));
    cmq(T>=perT) = cm1(T>=perT);
    cmq(T< perT) = cm2(T< perT);
    
    cxq = min(c,cphs0 + cxq.*(cphs1-cphs0));
    cmq = max(c,cphs0 + cmq.*(cphs1-cphs0));
    
    xq  = max(TINY,min(1-TINY, (c-cmq)./(cxq-cmq) ));

if size(xq,1)>1
    xq([1 end],:) = xq([2 end-1],:);
    xq(:,[1 end]) = xq(:,[2 end-1]);
    cxq([1 end],:) = cxq([2 end-1],:);
    cxq(:,[1 end]) = cxq(:,[2 end-1]);
    cmq([1 end],:) = cmq([2 end-1],:);
    cmq(:,[1 end]) = cmq(:,[2 end-1]);
end

