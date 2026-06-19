% Compares your latest commit with CO2SYS v3.2.1



%% Run CO2SYS v3.2.1
disp('Running CO2SYS v3.2.1')
rehash path % Circumvent caching of CO2SYS function
system('git checkout v3.2.1'); % Checkout the v3.2.1 tag

%% Set up input conditions
PARvalues = [2250 2100 8.1 400 405];
PARTYPEs = 1:5;
pHSCALEIN_opts = 1:4;
K1K2CONSTANTS_opts = 1:15;
KSO4CONSTANTS_opts = 1:4;
KFCONSTANT_opts = 1;
SALvalue = 33.1;
[P1, P2, P1type, P2type, sal, pHscales, K1K2, ~, KSO4, KF, ...
    BSal, ~, ~] = CO2SYSigen(PARvalues, PARTYPEs, SALvalue, pHSCALEIN_opts, ...
    K1K2CONSTANTS_opts, KSO4CONSTANTS_opts, KFCONSTANT_opts);
tempin = 24;
tempout = 12;
presin = 1;
presout = 1647;
si = 10;
phos = 1;

%% Determine whether to calculate each input row or not
% xrow = 1 + 210; % just do one row, or...
xrow = 1:numel(P1); % ... do all rows (do this for saving output file)
P1 = P1(xrow);
P2 = P2(xrow);
P1type = P1type(xrow);
P2type = P2type(xrow);
sal = sal(xrow);
pHscales = pHscales(xrow);
K1K2 = K1K2(xrow);

tic
[DATA_v3, HEADERS_v3] = ...
    CO2SYS(P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
toc



%% Run current commit CO2SYS
disp('Running CO2SYS current commit/head')
rehash path % Circumvent caching of CO2SYS function
system('git switch -'); % Switch back to your commit

%% Set up input conditions
PARvalues = [2250 2100 8.1 400 405];
PARTYPEs = 1:5;
pHSCALEIN_opts = 1:4;
K1K2CONSTANTS_opts = 1:15;
KSO4CONSTANTS_opts = 1:4;
KFCONSTANT_opts = 1;
SALvalue = 33.1;
[P1, P2, P1type, P2type, sal, pHscales, K1K2, ~, KSO4, KF, ...
    BSal, ~, ~] = CO2SYSigen(PARvalues, PARTYPEs, SALvalue, pHSCALEIN_opts, ...
    K1K2CONSTANTS_opts, KSO4CONSTANTS_opts, KFCONSTANT_opts);
tempin = 24;
tempout = 12;
presin = 1;
presout = 1647;
si = 10;
phos = 1;

%% Determine whether to calculate each input row or not
% xrow = 1 + 210; % just do one row, or...
xrow = 1:numel(P1); % ... do all rows (do this for saving output file)
P1 = P1(xrow);
P2 = P2(xrow);
P1type = P1type(xrow);
P2type = P2type(xrow);
sal = sal(xrow);
pHscales = pHscales(xrow);
K1K2 = K1K2(xrow);

tic
[DATA, HEADERS] = ...
    CO2SYS(P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
toc

%% Compare results
fprintf("\n\nRelative change vs v3.2.1:\n"); ...
fprintf("%20s %20s %20s %20s %20s\n", "Variable", "Mean rel. change", "Min rel. change", "Max rel. change", "# of samples"); ...
for V = 1:length(HEADERS_v3)
    x = DATA_v3(:,V);
    y = DATA(:,V);
    % Only compare non fill-in values
    ix = x ~= -999 & y ~=(-999) & ~isnan(x) & ~isnan(y);
    relerr = abs(y(ix) - x(ix)) ./ abs(x(ix));
    % Only calculate metrics with non-NaN and non-Inf values
    ix2 = ~isinf(relerr) & ~isnan(relerr);
    maxrelerr = max(relerr(ix2));
    minrelerr = min(relerr(ix2));
    meanrelerr = mean(relerr(ix2));
    fprintf("%20s %20.2g %20.2g %20.2g %20i\n", HEADERS_v3{V}, meanrelerr, minrelerr, maxrelerr, sum(ix2))
end
