%% numerical radiogenic heating
% current version does not implement partitioning of radioactive isotopes

%% 
thalf       = 7.17e5*NUM.yr; % half life(s)
lambda      = log(2)/thalf;

%% numerically solve radioactive decay

% get values from previous solution

% numerically solve decay
dNdt = -lambda.*NAlo;

NAl  = NAlo + (NUM.theta.*dNdt + (1-NUM.theta).*dNdto) .* NUM.dt;

% calculate radiogenic heat production
MAT.Hr = -dNdt.*eAl;

