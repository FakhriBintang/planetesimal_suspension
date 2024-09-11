function [dndt,H] = rad_decay(npar,tau,E_decay)
dndt = -npar/tau;
H    = abs(dndt *E_decay); % calculates radiactive heating in W/Kg
end