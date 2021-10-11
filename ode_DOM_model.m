function [Nt] = ode_DOM_model(t, N, C, E, S, eta, beta, r_mort, r_max, K, xB)

% (bacteria: rows, DOM : columns)
[numB, numD] = size(C);

B    = N(1      : numB);
D    = N(numB+1 : numB+numD)';
AgeB = N(numB+numD+2 : end-numD);
AgeD = N(end-numD+1: end)';

% Prevent zero biomass
if any(B==0)
    warning(sprintf('\nBacterial biomass got zero! Corrected to 1e-200.'))
    B(B==0) = 1e-300;
end

% bacterial- and DOM-specific consumption  
consumption = B*((r_max*D)./(K + D)).*C;
consumption_1  = sum(consumption,1);
consumption_2  = sum(consumption,2);

growth      = eta*consumption_2;
excretion   = ((1-eta)*beta*consumption_2')*E;
mortality   = r_mort(1).*B + r_mort(2).*B.*B;
lysis       = mortality'*E;

% ODE's
dB = growth - mortality;
dD = excretion + lysis - consumption_1 + S;
dC = (1-eta)*(1-beta)*sum(consumption_1);

% conditional breakpoints:
% (t>500) && any(dD>0) at t = 870
%if (t>500) && any(dD>0) 
% if (t>300) 
%     t
%     keyboard()
% end

%% Age 

f_excreted = excretion./D;
f_excreted(isinf(f_excreted)) = 0;
f_excreted(isnan(f_excreted)) = 0;
f_lysed = lysis./D;
f_lysed(isinf(f_lysed)) = 0;
f_lysed(isnan(f_lysed)) = 0;
f_supplied = S./D;
f_supplied(isinf(f_supplied)) = 0;
f_supplied(isnan(f_supplied)) = 0;

% excreted age
Cnorm = consumption./repmat(consumption_2, 1, size(C,2));
Cnorm(isnan(Cnorm)) = 0;
A_cons = Cnorm*AgeD';
Enorm  = E./repmat(sum(E),size(E,1),1);
A_excreted = A_cons'*Enorm;

% lysed age
dAgeB    = 1 - (growth./B).*(AgeB-A_cons);
A_lysed  = AgeB'*Enorm;

% age ODE
Age0  = 0; % age of supplied compounds (~DIC age at surface)
dAgeD = 1 - f_excreted.*(AgeD-A_excreted) - f_lysed.*(AgeD-A_lysed) - f_supplied.*(AgeD-Age0);

% Bacterial inflow // Technical variable to prevent biomass decline
%dB = dB +  1e-50 * (B < 1e-100);
% dB = dB + 1e-20;
dB = dB + xB;

Nt= [dB; dD'; dC; dAgeB; dAgeD'];


end
