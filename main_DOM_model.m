%% Main file for DOM-bacteria-network model

clear variables
close all

%% PARAMETERS AND MODEL SPECIFICATIONS

% Hold the random number generator seed for reproducible results
rng('default')
rng(1)

% Model parameters
numB    = 35; % number of bacterial groups
numD    = 100; % number of DOC compound groups
K       = 10; % half-saturation constant
nexcr   = 30; % number of excreted compounds per bacterial group
nsubs   = 3; % number of compounds taken-up per bacterial group
eta     = 0.2; % growth efficiency of bacteria
beta    = 0.14/(1-eta); % fraction of uptake released back to DOC pool via excretion
r_max   = 1; % maximum growth rate of bacteria
r_mort = [0.02 0]; % mortality rate of bacteria (second entry is for quadratic mortality term - not used here)

% Generate uptake- and release network 
% consumption matrix (bacteria in rows, DOC in columns)
C = get_consumption_matrix(numB, numD, 'universe',...
    nsubs, 'non-normalized');
% excretion matrix (bacteria in rows, DOC in columns)
E = get_excretion_matrix(C, nexcr);

% Supply
Stot        = 0.08; % total default supply [mmolC/m^3/d]
S           = repmat(Stot/numD, 1, numD); % even distribution of total supply across compounds

% Initial condition
sumB0      = 1; % total intital carbon concentration of microbes [mmolC/m^3]
sumD0      = 80; % total intital DOC concentration [mmolC/m^3]
B0         = repmat(sumB0/numB, numB, 1); % all bacterial groups start with the same share of carbon concentration
D0         = repmat(sumD0/numD, numD, 1); % all DOC compound groups start with the same share of carbon concentration

% Solver options
maxstep = 50; % maximum step length
reltol  = 1e-06; % relative error tolerance
abstol  = 1e-20; % absolute error tolerance

%% Minimal working example: Default simulation

% simulate for ten years
tspan = [0 10*365];

% Default simulation, returns:
% - t, the time vector (t x 1)
% - Bt, the carbon concentration of bacterial groups over time [mmolC/m^3] (t x numB)
% - Dt, the carbon concentration of DOC groups over time [mmolC/m^3] (t x numD)
% - It, the carbon concentration of the inorganic pool over time [mmolC/m^3] (t x 1)
% - At_D, the age of carbon in the DOC compound groups over time [days] (t x numD)
% - At_D, the age of carbon in the bacterial biomass over time [days] (t x numB)
[t, Bt, Dt, It, At_D, At_B] = wrap_ode_DOM_model(tspan, C, E, S, ...
    beta, eta, r_mort, r_max, K, B0, D0, 'plot', 'on', 'reltol', reltol, 'maxstep', maxstep);

% Total DOC concentration at the end of the simulation in [mmolC/m^3]
sum(Dt(end,:))

% Total carbon concentration of bacteria at the end of the simulation in [mmolC/m^3]
sum(Bt(end,:))

% Total carbon concentration in the inorganic carbon pool (I) at the end of the simulation in [mmolC/m^3]
It(end)

% Concentration of first DOC compound at after 1 year in [mmolC/m^3]
idx = find(t>365, 1, 'first');
Dt(idx,1)

% Concentration of third bacterial group after 1 year in [mmolC/m^3]
Bt(idx,3)

% Age of second DOC compound group at the end of the simulation in years
% (The DOC is younger than the simulation time due to the supply of DOC)
At_D(end,2)/365

% Age of carbon in biomass of fifth bacterial group at the end of the simulation in years
At_B(end,5)/365

% Age distribution of individual compound groups at the end
figure(); histogram(At_D(end,:)/365)
xlabel('Age at the end of the simulation in years')

% Average concentration-weighted age of DOC in years
meanAge  = sum(At_D.*Dt./repmat(sum(Dt,2), 1, numD),2)/365;
figure(); plot(t/365, meanAge)
xlabel('simulation time [years]'), ylabel('mean age of DOC [years]')

%% Same simulation without DOC supply (--> linear aging)

% no supply
S_none = zeros(1, numD);

% simulate for ten years
tspan = [0 10*365];

% Simulation without supply
[t, Bt, Dt, It, At_D, At_B] = wrap_ode_DOM_model(tspan, C, E, S_none, ...
    beta, eta, r_mort, r_max, K, B0, D0, 'plot', 'on', 'reltol', reltol, 'maxstep', maxstep);

% Total DOC concentration at the end of the simulation in [mmolC/m^3]
sum(Dt(end,:))

% Total carbon concentration of bacteria at the end of the simulation in [mmolC/m^3]
sum(Bt(end,:))

% Total carbon concentration in the inorganic carbon pool (I) at the end of the simulation in [mmolC/m^3]
It(end)

% Concentration of first DOC compound at after 1 year in [mmolC/m^3]
idx = find(t>365, 1, 'first');
Dt(idx,1)

% Concentration of third bacterial group after 1 year in [mmolC/m^3]
Bt(idx,3)

% Age of second DOC compound group at the end of the simulation in years
% (The DOC is younger than the simulation time due to the supply of DOC)
At_D(end,2)/365

% Age distribution of individual compound groups at the end
figure(); histogram(At_D(end,:)/365)
xlabel('Age at the end of the simulation in years')

% Average concentration-weighted age of DOC in years
meanAge  = sum(At_D.*Dt./repmat(sum(Dt,2), 1, numD),2)/365;
figure(); plot(t/365, meanAge)
xlabel('simulation time [years]'), ylabel('mean age of DOC [years]')


%% Working example II: Generate data for Figure 6 ("Toy model set-up")

% Simple toy set-up
numB_simp = 4;
numD_simp = 5;

% Prescribe consumption and excretion matrix
C_simp = [1 0 0 0 1; 0 1 0 0 1; 0 1 1 0 0; 0 0 0 1 1];
E_simp = [0 0.9 0.1 0 0; 0 0 0.9 0.1 0; 0 0 0 0.8 0.2; 0.6 0 0.4 0 0];

% Check matrices
assert(all(all(C_simp+E_simp<=1)), 'microbes should not take up compounds produced by themselves')

% No supply of DOC
S_simp = zeros(1, numD_simp);

% Increase fraction released for illustration purpose (to be able to see
% the re-working of DOC compounds through the network)
beta_simp = 2.2*0.14/(1-eta); 

% leave the other parameters at default level
r_mort_simp = r_mort;
K_simp = K; 

% solver time span
tspan_simp = [0, 600];

% start with equal amounts of all bacteria
B0_simp    = repmat(0.02, 1, numB_simp);

% start with only the first DOC compound present at high concentrations
D0_simp    = [20 zeros(1, numD_simp-1)];

% Simulate toy model
[t_simp, Bt_simp, Dt_simp, Ct_simp, At_simp] = wrap_ode_DOM_model(tspan_simp,...
    C_simp, E_simp, S_simp, beta_simp, eta, r_mort_simp, r_max, K_simp, ...
    B0_simp, D0_simp, 'plot', 'off', ...
    'reltol', reltol, 'abstol', abstol, 'maxstep', maxstep);


%% Plot Figure 6 (Toy model)

figure('color', 'white', 'position', [328,294,365,480])

subplot(3,2,1)
imagesc(C_simp)
caxis([0 1])
ax = gca;
ax.XTick = 1:numD_simp;
ax.YTick = 1:numB_simp;
ax.YTickLabel = arrayfun(@(i) sprintf('B_%d', i), 1:numB_simp, 'Unif', false);
ax.XTickLabel = arrayfun(@(i) sprintf('D_%d', i), 1:numD_simp, 'Unif', false);
textStrings = num2str(C_simp(:),'%1.0f');  
textStrings = strtrim(cellstr(textStrings)); 
[x,y] = meshgrid(1:numD_simp, 1:numB_simp);  
hStrings = text(x(:),y(:),textStrings(:),...    
                'HorizontalAlignment','center');
set(hStrings(C_simp(:) > 0), 'FontWeight', 'bold');
set(hStrings(C_simp(:) == 0), 'Color', [.5 .5 .5]);
title('Uptake matrix {\bf{\itU}}')
axU = gca;
axis equal
a = gca;
b = copyobj(gca, gcf);
set(gca, 'Box', 'off', 'YColor', [1 1 1], 'XColor', [1 1 1], ...
    'XTickLabel', [], 'YTickLabel', [])
set(b, 'Xcolor', 'k', 'YColor', 'k', 'Box', 'off')
b.Title.String = '';
uistack(b,'down')

subplot(3,2,2)
imagesc(E_simp)
caxis([0 1])
colormap([1 1 1]), axis square
ax = gca;
ax.XTick = 1:numD_simp;
ax.YTick = 1:numB_simp;
ax.YTickLabel = arrayfun(@(i) sprintf('B_%d', i), 1:numB_simp, 'Unif', false);
ax.XTickLabel = arrayfun(@(i) sprintf('D_%d', i), 1:numD_simp, 'Unif', false);
textStrings = num2str(E_simp(:),'%1.1f');  
textStrings = strtrim(cellstr(textStrings)); 
for i = 1:length(textStrings)
    if strcmp(textStrings(i), '0.0')
        textStrings(i) = {'0'};
    elseif strcmp(textStrings(i), '1.0')
        textStrings(i) = {'1'};
    end
end
[x,y] = meshgrid(1:numD_simp, 1:numB_simp);  
hStrings = text(x(:),y(:),textStrings(:),...    
                'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim')); 
set(hStrings(E_simp(:) > 0), 'FontWeight', 'bold');
set(hStrings(E_simp(:) == 0), 'Color', [.5 .5 .5]);
title('Release matrix {\bf{\itR}}')
axR = gca;
axis equal
a = gca;
b = copyobj(gca, gcf);
set(gca, 'Box', 'off', 'YColor', [1 1 1], 'XColor', [1 1 1], ...
    'XTickLabel', [], 'YTickLabel', [])
set(b, 'Xcolor', 'k', 'YColor', 'k', 'Box', 'off')
b.Title.String = '';
uistack(b,'down')

%%

subplot(3,2, [3 4])
hold(gca, 'all')
gcols = [1 0 0; 0 0 0; .3 .3 .3; .45 .45 .45; .7 .7 .7; .85 .85 .85];
ps = plot(t_simp, sum(Dt_simp,2), 'color', 'green', 'linewidth', 2);
set(gca, 'ColorOrder', gcols, 'YLim', [-0.5 sum(D0_simp)*1.22])
pi = plot(t_simp, Dt_simp, 'linewidth', 1.7);
uistack(pi(3), 'top')
uistack(pi(5), 'top')
[ld, ldt] = legend([pi; ps], [arrayfun(@(i) sprintf('D_%d', i), 1:numD_simp, 'Unif', false),...
    'total']);
yd = ylabel('DOC concentration');
axD = gca;
axD.XTick = [];
axD.YTick = [];

subplot(3,2, [5 6])
hold(gca, 'all')
set(gca, 'ColorOrder', gray(numB_simp+2), 'YLim', [-0.05 3.1])
ps = plot(t_simp, sum(Bt_simp,2), 'color', 'red', 'linewidth', 2);
pi = plot(t_simp, Bt_simp, 'linewidth', 1.7);
[lb, lbt] = legend([pi;ps], [arrayfun(@(i) sprintf('B_%d', i), 1:numB_simp, 'Unif', false),...
    'total']);
yb = ylabel('Biomass concentration');
xlabel('Time')
axB = gca;
axB.XTick = [];
axB.YTick = [];
set(findall(gcf,'-property','FontName'), 'FontName','Droid Sans', 'Fontsize', 8) 

yb.Position(1) = yd.Position(1);
yd.Position(1) = yd.Position(1);
axB.XLabel.Position(2) = -0.4;
axB.YLabel.Position(1) = -25;
axD.YLabel.Position(1) = -25;
axD.Position(2) = 0.38;

set(ldt(1:6), 'Fontsize', 8, 'FontName', 'Droid Sans')
set(lbt(1:5), 'Fontsize', 8, 'FontName', 'Droid Sans')
set(ld, 'position', [0.7019 0.42 0.2027 0.2173]);
set(lb, 'position', [0.7019 0.1451 0.2027 0.1806]);
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 8)





