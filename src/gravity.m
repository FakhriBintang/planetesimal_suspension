function [gz, gzP] = gravity(rho,zW)
    rw = flip(zW)'; % convert depth to radial distance
    % %% trapezoidal method
    % rhoW = (rho(1:end-1,:)+ rho(2:end,:))./2;
    % m_trapz = flip(4*pi()/3.*cumtrapz(flip(rw).^3,flip(rhoW)));

    %% shell method
    m_shell = rho(2:end-1,:)*4*pi()/3.*(rw(1:end-1).^3-rw(2:end).^3);
    m_sh = [flip(cumsum(flip(m_shell))); 0 0 0]; 

    %% calculate graavity
    G = 6.672e-11;
    gmin = 0.01;
    % gz_t = G.*m_trapz./(rw.^2); gz_t(isnan(gz_t)) = 0;
    gz = G.*m_sh./(rw.^2);   gz(isnan(gz)) = 0; gz = gz+gmin;
    % gzU = (gz(1:end-1,1:end-1)+gz(2:end,1:end-1)+gz(1:end-1,2:end)+gz(2:end,2:end))./4;
    % gzU = [gzU(1,:); gzU; 0 0];
    gzP = (gz(1:end-1,:)+gz(2:end,:))./2;
    gzP = [gzP(1,:); gzP; gzP(end,:)];
end