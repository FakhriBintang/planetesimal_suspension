%% test effect of melt frraction 
rhoxSi = 3300;
rhomSi = 2900;
rhoxFe = 6800;
rhomFe = 5900;

% initial silicate melt fraction
phimSi = linspace(0,100,101);

% get segregation velocities
vxSi_seg = 2.*(rhoxSi-rhoTot).*(d^2)*g./9./Eta; % solid Si
vxFe_seg = 2.*(rhoxFe-rhoTot).*(d^2)*g./9./Eta; % solif Fe
vmFe_seg = 2.*(rhomFe-rhoTot).*(d^2)*g./9./Eta; % molten Fe

