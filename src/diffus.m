function [dff] = diffus(mode,a,k,h,rp)
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
        qz    = - (k(1:end-1,2)+k(2:end,2))./2 .* (rw.^2) .* ddz(a(:,2),h);
        dff   = - ddz(qz,h)./(rp(2:end-1).^2);
end

end