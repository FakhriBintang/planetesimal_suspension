function [dndt,H] = rad_decay(npar,tau,E_decay)
dndt = -npar/tau;
H    = abs(dndt *E_decay);
end