function [xq,cxq,cmq]  =  equilibrium_noPer(T0,c,P,Tphs0,Tphs1,cphs0,cphs1,clap,PhDg,TINY)

T   = max(0,min(1,(T0 - P*clap -Tphs0)./(Tphs1-Tphs0))) ;

% cx1 = max(TINY,min(1-TINY,          perCx .*erfc((2+PhDg).*(T-perT)./(1-perT)))); 
% cx2 = max(TINY,min(1-TINY, perCx+(1-perCx).*erfc((1+PhDg).*(T     )./   perT) ));
cxq = max(TINY,min(1-TINY,          erfc((2+PhDg).*T))); 


% cxq = zeros(size(T));
% cxq(T>=perT) = cx1(T>=perT);
% cxq(T< perT) = cx2(T< perT);

cmq = max(TINY,min(1-TINY,          perCm .*erf((1.0+PhDg/10).*(1   -T))./erf((1.0+PhDg/10))));
% cm1 = max(TINY,min(1-TINY,          perCm .*erf((1.0+PhDg/10).*(1   -T)./(1-perT))./erf((1.0+PhDg/10))));
% 
% cm2 = max(TINY,min(1-TINY, perCm+(1-perCm).*erf((0.9+PhDg/10).*(perT-T)./(  perT))./erf((0.9+PhDg/10))));

% cmq = zeros(size(T));
% cmq(T>=perT) = cm1(T>=perT);
% cmq(T< perT) = cm2(T< perT);

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

