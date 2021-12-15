% planetesimal: solve thermo-chemical equations

% print solver header
fprintf(1,'  ---  solve thermo-chemical equations \n');
tic;

%% update phase fractions
WFe                 = SOL.segW + SOL.W;
[advn]              = advection(SOL.phi,SOL.U,WFe,NUM.h,NUM.h,NUM.AdvnScheme,'Flx');
if NUM.step>0
SOL.phi             = phio - advn * NUM.dt;
end
SOL.phi(SOL.phi<0)  = 0;
%% Solve energy equation explicitly

% get advection rate
advn_T              = advection(SOL.T,SOL.U,SOL.W,NUM.h,NUM.h,NUM.AdvnScheme,'adv');
% get diffusion rate
diff_T              = (diff(SOL.T(:,2:end-1),2,1)./NUM.h^2  ...
                    +  diff(SOL.T(2:end-1,:),2,2)./NUM.h^2) ...
                    .* MAT.kT(2:end-1,2:end-1);

% heat capacity density
RhoCp               = MAT.Rhot.*MAT.Cp;

% get total rate of change
SOL.dTdt            = - advn_T(2:end-1,2:end-1) ...
                      + (diff_T + MAT.Hr(2:end-1,2:end-1) ...
                      + SOL.Hs(2:end-1,2:end-1) ...
                      + SOL.Ha(2:end-1,2:end-1))./RhoCp(2:end-1,2:end-1);
SOL.dTdt(isinf(SOL.dTdt)) = 0;


% update temperature solution
if NUM.step>0
SOL.T(2:end-1,2:end-1) = To(2:end-1,2:end-1) + (  NUM.theta .*SOL.dTdt   ...
    + (1-NUM.theta).*    dTdto) .* NUM.dt;
end

% apply top boundary conditions
switch SOL.BCTTop
    case 'isothermal'; SOL.T(1,:) =     To(1,:);
    case 'insulating'; SOL.T(1,:) = SOL.T (2,:);
end

% apply bottom boundary conditions
switch SOL.BCTBot
    case 'isothermal'; SOL.T(end,:) =     To(end  ,:);
    case 'insulating'; SOL.T(end,:) = SOL.T (end-1,:);
end

% apply side boundary conditions
switch SOL.BCTSides
    case 'isothermal'; SOL.T(:,[1 end]) =     To(:,[1 end  ]);
    case 'insulating'; SOL.T(:,[1 end]) = SOL.T (:,[2 end-1]);
end