%% Wrap function for ODE of DOM model
%
% Optional input arguments:
% - string 'plot', followed by string-code for different types of plots
%          'tBD', for time development of B and D
%          'Ogawa', for time development of D compared to Ogawa experiment
%          'mass_conservation', for time development total carbon of system
%
% Output arguments:
% t:     time vector 
% Bt:    time series of bacterial biomass 
% Dt:    time series of compound concentration 
% Ct:    time series of inorganic carbon (CO2) concentration
% At_D : time series of DOM age
% At_B : time series of bacterial biomass age
%
% Number of rows = number of time points.
%
function [t, Bt, Dt, Ct, At_D, At_B] = wrap_ode_DOM_model(tspan, C, E, S, beta, eta, r_mort, r_max, K, B0, D0, varargin)
%% Check and assign input

numB = size(C,1);
numD = size(C,2);

% r_mort: specify as vector [linear mortality, quadratic mortaliy]
assert(isequal(size(r_mort), [1 2]),...
    'Wrong dimensions of r_mort. Give as: [r_mort1 r_mort2].')
assert(isequal(size(S), [1 numD]) | isequal(S, 0),...
    'Wrong dimensions of S.')

% Relative error tolerance for ODE
reltol  = 1e-3; % Default value
if any(strcmp(varargin, 'reltol'))
    ind  = find(strcmp(varargin, 'reltol'));
    reltol = varargin{ind+1};
end

% Absolute error tolerance for ODE
abstol  = 1e-6; % Default value
if any(strcmp(varargin, 'abstol'))
    ind  = find(strcmp(varargin, 'abstol'));
    abstol = varargin{ind+1};
end

% Maximum step size for ODE
maxstep  = 50; %0.1*abs(tspan(end)-tspan(1)); % = Default value
if any(strcmp(varargin, 'maxstep'))
    ind  = find(strcmp(varargin, 'maxstep'));
    maxstep = varargin{ind+1};
end

% Initial age of bacterial biomass and compounds
A0 = zeros(numB+numD,1);
if any(strcmp(varargin, 'Age0'))
    ind  = find(strcmp(varargin, 'Age0'));
    A0 = varargin{ind+1};
end

% Technical variable: inflow of B (to prevent biomass decline to e-300)
% set to zero to get results as before June 2018
xB = 1e-300;
% xB = 0
% xB = 1e-100
if any(strcmp(varargin, 'xB'))
    ind  = find(strcmp(varargin, 'xB'));
    xB = varargin{ind+1};
end

%% Call ODE solver

% Default starting conditions for A and fluxes
C0 = 0;
y0 = [B0(:); D0(:); C0; A0(:)]; 

% ODE solver (use adapted solver which returns non-interpolated time steps, 
% otherwise age calculations can be imprecise)
options = odeset('NonNegative', 1:size(y0,1), 'RelTol', reltol, ...
    'AbsTol', abstol, 'MaxStep', maxstep);
ind = find(strcmp(varargin, 'tvalues'));
if ~isempty(ind) && length(varargin)>=ind+1 && ...
        isa(varargin{ind+1}, 'double') 
    % force the specified time steps (non-interpolated)
    tvalues = varargin{ind+1};
    [t, Nt] = ode45_(@ode_DOM_model, tvalues, y0, options,...
        C, E, S, eta, beta, r_mort, r_max, K, xB);
else % let solver choose timesteps (non-interpolated)
    [t, Nt] = ode45_(@ode_DOM_model, tspan, y0, options,...
        C, E, S, eta, beta, r_mort, r_max, K, xB);
end
   
Bt = Nt(:, 1      : numB);
Dt = Nt(:, numB+1 : numB+numD);
Ct = Nt(:, numB+numD+1);
At_B = Nt(:, numB+numD+2 : numB+numD+1+numB);
At_D = Nt(:, end-numD+1 : end);

% mass conservation?
mass_t = sum(Bt,2)+sum(Dt,2)+Ct-t*sum(S,2);
if all(abs(mass_t(1)-mass_t)<1e-5)
    % all good
elseif all(mass_t(1)-mass_t<1e-4)
    warning('Mass conservation hampered (only <1e-4))')
    figure('color', 'white', 'position', [636,99,720,665])
            subplot(3,1,1)
            plot(t, sum(Bt,2)+sum(Dt,2)+Ct)
            xlabel('time (days)'), ylabel('Total carbon')
            axis tight
            subplot(3,1,2)
            plot(t, t*sum(S,2))
            xlabel('time (days)'), ylabel('Cumulative supply')
            axis tight
            subplot(3,1,3)
            plot(t, sum(Bt,2)+sum(Dt,2)+Ct-t*sum(S,2))
            xlabel('time (days)'), ylabel('Total carbon minus supply')
            axis tight
else
    error('Not mass conserving.')
    figure('color', 'white', 'position', [636,99,720,665])
            subplot(3,1,1)
            plot(t, sum(Bt,2)+sum(Dt,2)+Ct)
            xlabel('time (days)'), ylabel('Total carbon')
            axis tight
            subplot(3,1,2)
            plot(t, t*sum(S,2))
            xlabel('time (days)'), ylabel('Cumulative supply')
            axis tight
            subplot(3,1,3)
            plot(t, sum(Bt,2)+sum(Dt,2)+Ct-t*sum(S,2))
            xlabel('time (days)'), ylabel('Total carbon minus supply')
            axis tight
end

%% Warn if bacteria die out and then come back (zombie bacteria)

die_out = NaN(1,size(Bt,2));
revival = zeros(1,size(Bt,2));
for x = 1:size(Bt,2)
    idx = find(Bt(:,x)==0, 1, 'first');
    if ~isempty(idx)
        die_out(x) = idx;
        revival(x) = any(Bt(die_out(x):end,x)>0);
    end
end
if any(revival)
    warning('Bacteria die out and then come back!')
end


%% Plotting

if any(strcmp(varargin, 'plot'))
    ind  = find(strcmp(varargin, 'plot'));
    switch varargin{ind+1}
        case {'tBD', 'BDt', 'on'}
            figure('color', 'white', 'position', [660,205,800,335])
            [ax] = plotyy(t, sum(Dt,2), t, sum(Bt,2));
            axis(ax, 'tight'), xlabel('time (days)')
            ax(1).YLabel.String = 'D mmolC m^{-3}';
            ax(2).YLabel.String = 'B mmolC m^{-3}';
        case {'age', 'all'}
            figure('color', 'white', 'position', [660,205,800,335])
            subplot(2,1,1)
            [ax] = plotyy(t, sum(Dt,2), t, sum(Bt,2));
            axis(ax, 'tight'), xlabel('time (days)')
            ax(1).YLabel.String = 'D mmolC m^{-3}';
            ax(2).YLabel.String = 'B mmolC m^{-3}';
            subplot(2,1,2)
            meanAge = sum(At_D.*Dt./repmat(sum(Dt,2), 1, numD),2)/365;
            plot(t, meanAge); 
            set(gca, 'XLim', [min(t) max(t)], 'YLim', [0 inf])
            xlabel('time (days)'), ylabel('mean age [y]')
        case {'tBD_Ogawa', 'Ogawa'}
            figure('color', 'white', 'position', [660,205,800,335])
            plot(t, sum(Dt,2), 'color', mycolors('green'), 'linewidth', 2);
            axis('tight'), xlabel('time (days)'), ylabel('mmol C m⁻³')
            gluct = [0,2,4,7,365];
            glucD = [208,30,18,16,11];
            hold on
            plot(gluct, glucD, 'r+', 'linewidth', 2)
            legend('Model DOM', 'DCta DOC')
        case {'mass_conservation', 'mass_conserving'}
            figure('color', 'white', 'position', [636,99,720,665])
            subplot(3,1,1)
            plot(t, sum(Bt,2)+sum(Dt,2)+Ct)
            xlabel('time (days)'), ylabel('Total carbon')
            axis tight
            subplot(3,1,2)
            plot(t, t*sum(S,2))
            xlabel('time (days)'), ylabel('Cumulative supply')
            axis tight
            subplot(3,1,3)
            plot(t, sum(Bt,2)+sum(Dt,2)+Ct-t*sum(S,2))
            xlabel('time (days)'), ylabel('Total carbon minus supply')
            axis tight   
        case {'fluxes'}
            figure('position', [675,541,1168,423])
            subplot(2,4,1)
            hold on
            plot(t, sum(Dt,2), 'linewidth', 2), axis tight, xlabel('days'),title('Dt')
            plot(t, Dt2, 'k:')
            plot(t, Dst2, 'b--')
            legend('Dt', 'sum(fluxes)')
            title('DOM')
            subplot(2,4,2)
            plot(t, fluxes.cons, 'linewidth', 2), axis tight, xlabel('days'), title('consumption')
            subplot(2,4,3)
            plot(t, fluxes.grow, 'linewidth', 2), axis tight, xlabel('days'), title('growth')
            subplot(2,4,4)
            plot(t, fluxes.excr, 'linewidth', 2), axis tight, xlabel('days'), title('excretion')
            subplot(2,4,5)
            plot(t, fluxes.resp, 'linewidth', 2), axis tight, xlabel('days'), title('respiration')
            subplot(2,4,6)
            plot(t, fluxes.supp, 'linewidth', 2), axis tight, xlabel('days'), title('supply')
            subplot(2,4,7)
            plot(t, fluxes.mort, 'linewidth', 2), axis tight, xlabel('days'), title('mortality')
        case {'individual'}
            figure()
            subplot(2,1,1)
            plot(t, Bt)
            set(gca, 'YScale', 'log')
            axis tight
            ylabel('Biomass [mmolC/m³]')
            subplot(2,1,2)
            plot(t, Dt)
            axis tight
            xlabel('Time [d]')
            ylabel('DOC [mmolC/m³]')
            
            
    end
end

end

