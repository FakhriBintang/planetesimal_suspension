function [dff] = advect_simple(mode,a,w,u,h,rp)
% a = quantity
% k = diffusivity
% h = grid spacing
% rp = radial distance (if spherical)
switch mode
    case 'cartesian'
        qz    = - (k(1:end-1,:)+k(2:end,:))./2 .* ddz(a,h);
        qx    = - (k(:,1:end-1)+k(:,2:end))./2 .* ddx(a,h);
        dff   = - ddz(qz(:,2:end-1),h)  ...                                      % heat diffusion
                - ddx(qx(2:end-1,:),h);
    case 'spherical'
        % ignore horizontal component
        rw    = (rp(1:end-1)+rp(2:end))./2;
        qz    = - (a(1:end-1)+a(2:end))./2 .* (rw.^2) .* w;
        dff   = - ddz(qz,h)./(rp(2:end-1).^2);
end

end