function [advn] = advection(a,u,w,dx,dz,scheme,type)

wp  = w(2:end  ,2:end-1);
wm  = w(1:end-1,2:end-1);
up  = u(2:end-1,2:end  );
um  = u(2:end-1,1:end-1);

vx  = (up+um)./2;
vz  = (wp+wm)./2;

vxp = max(vx,0); vxm = min(vx,0);
vzp = max(vz,0); vzm = min(vz,0);

Div_v = diff(w(:,2:end-1),1,1)./dz + diff(u(2:end-1,:),1,2)./dx;

agh                    = zeros(size(a)+2);
agh(2:end-1,2:end-1)   = a;

agh([1 2 end-1 end],:) = agh([3 3 end-2 end-2],:);
agh(:,[1 2 end-1 end]) = agh(:,[3 3 end-2 end-2]);

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);


switch scheme
    
    case 'centred FD' % not recommended
        
        advn = vx.*(aip-aim)./h + vz.*(ajp-ajm)./h;
        
        
    case 'flxdiv'
        
        advn = ((ajp+acc)./2.*wp - (ajm+acc)./2.*wm)./dz ...
             + ((aip+acc)./2.*up - (aim+acc)./2.*um)./dx;
        
        
    case 'fromm'
        
        advn   =     up .*(-(aipp-aip)./dx./8 + (aip + acc)./dx./2 + (acc-aim )./dx./8) ...
               - abs(up).*(-(aipp-aip)./dx./8 + (aip - acc)./dx./4 - (acc-aim )./dx./8) ...
               -     um .*(-(aip -acc)./dx./8 + (acc + aim)./dx./2 + (aim-aimm)./dx./8) ...
               + abs(um).*(-(aip -acc)./dx./8 + (acc - aim)./dx./4 - (aim-aimm)./dx./8) ...
               +     wp .*(-(ajpp-ajp)./dz./8 + (ajp + acc)./dz./2 + (acc-ajm )./dz./8) ...
               - abs(wp).*(-(ajpp-ajp)./dz./8 + (ajp - acc)./dz./4 - (acc-ajm )./dz./8) ...
               -     wm .*(-(ajp -acc)./dz./8 + (acc + ajm)./dz./2 + (ajm-ajmm)./dz./8) ...
               + abs(wm).*(-(ajp -acc)./dz./8 + (acc - ajm)./dz./4 - (ajm-ajmm)./dz./8);
        
        
    case 'first upwind'
        
        axp   = (aip-acc)./dx;
        axm   = (acc-aim)./dx;
        azp   = (ajp-acc)./dz;
        azm   = (acc-ajm)./dz;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        advn  = daxdt + dazdt;
        
        
    case 'second upwind'
        
        axp   = (-3*acc+4*aip-aipp)/2/dx;
        axm   = ( 3*acc-4*aim+aimm)/2/dx;
        azp   = (-3*acc+4*ajp-ajpp)/2/dz;
        azm   = ( 3*acc-4*ajm+ajmm)/2/dz;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        advn  = vxp.*axm + vxm.*axp + vzp.*azm + vzm.*azp + acc.*Div_v;
        
        
    case 'third upwind'
        
        axp   = (-2*aim-3*acc+6*aip-aipp)/6/dx;
        axm   = ( 2*aip+3*acc-6*aim+aimm)/6/dx;
        azp   = (-2*ajm-3*acc+6*ajp-ajpp)/6/dz;
        azm   = ( 2*ajp+3*acc-6*ajm+ajmm)/6/dz;
        
        advn  = vxp.*axm + vxm.*axp + vzp.*azm + vzm.*azp + acc.*Div_v;
        
%         axp   = (-2*aim-3*acc+6*aip-aipp)/6/dx;
%         axm   = ( 2*aip+3*acc-6*aim+aimm)/6/dx;
%         azp   = (-2*ajm-3*acc+6*ajp-ajpp)/6/dz;
%         azm   = ( 2*ajp+3*acc-6*ajm+ajmm)/6/dz;
%         
%         daxdt = vxp.*axm + vxm.*axp;
%         dazdt = vzp.*azm + vzm.*azp;
        
%         advn  = daxdt + dazdt;
end

if strcmp(type,'adv')
    advn = advn - acc.*Div_v;
end

advn = advn([1 1:end end],[1 1:end end]);                         % periodic top/bot boundaries

end
