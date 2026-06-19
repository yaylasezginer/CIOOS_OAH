function henrys_law_const = eh_henryconst(S,T,gas)

% eh_henryconst   Henry's Law constant for a gas in seawater at equilibrium
%             with the atmosphere
%=========================================================================
% 
% USAGE:  henrys_law_const = eh_henryconst(S,T,gas)
%
% DESCRIPTION:
%    Calculates the Henry's Law constant of a gas in seawater at the given 
%    salinity and temperature of the water and a hydrostatic pressure of 0 
%    dbar (surface).  
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%   gas = gas for which to calculate Henry's Law constant, possible values:
%       'N2','O2','Ar','CO2','Ne','He','CH4','Kr','N2O','Xe','CFC-12',
%       'CFC-11','SF6','DMS'
%
% OUTPUT:
%   henrys_law_const = Henry's Law constant for the gas in seawater 
%       [mol kg^-1 atm^-1]
%
% REFERENCE:
%   Appendices D & F, Chemical Oceanography: Element Fluxes in the Sea
%       Emerson, S.R., R.C. Hamme 2022
%       Cambridge University Press
%   and references therein
%
% VERSION 1.0 : 23 March 2022
% AUTHOR: Roberta C. Hamme (University of Victoria) 
% This software is available from http://www.cambridge.org/emerson-hamme
% as part of Chemical Oceanography: Element Fluxes in the Sea (2022) 
% by Steven R. Emerson and Roberta C. Hamme
%=========================================================================

%----------------------
% Check input parameters
%----------------------

% check number of input parameters is correct
if nargin ~= 3
   error('Must pass 3 input parameters')
end %if

% check datatypes of input parameters are correct
validateattributes(S,{'numeric'},{'nonempty'},mfilename,'S:salinity')
validateattributes(T,{'numeric'},{'nonempty'},mfilename,'T:temperature')
validateattributes(gas,{'char'},{'nonempty'},mfilename,'gas')

% check S,T dimensions and verify they have the same shape or are singular
[rs,cs] = size(S);
[rt,ct] = size(T);
if ((rs~=rt) || (cs~=ct)) && (rs+cs>2) && (rt+ct>2)
   error('S & T must have same dimensions or be singular')
end %if

% check that gas is one of the supported values
expectedGases = {'N2','O2','Ar','CO2','Ne','He','CH4','Kr','N2O','Xe','CFC-12','CFC-11','SF6','DMS'};
if sum(strcmp(gas,expectedGases))==0
    error('Expected input parameter gas to match one of these values:\n''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', or ''%s''','N2','O2','Ar','CO2','Ne','He','CH4','Kr','N2O','Xe','CFC-12','CFC-11','SF6','DMS');
end %if

%------
% BEGIN
%------

% calculate scaled temperature
TS = log((298.15-T)./(273.15+T));
% calculate temperature in Kelvin
TK = T+273.15;

% perform calculation, mainly converting from eq conc to Henry's law const
switch gas
    case 'N2'
        molfract = 7.8084e-1;
        gas_equil_conc_umolkg = exp(6.42931 + 2.92704*TS + 4.32531*TS.^2 + 4.69149*TS.^3 + S.*(-7.44129e-3 - 8.02566e-3*TS  -1.46775e-2*TS.^2));
        henrys_law_const = gas_equil_conc_umolkg * 1e-6 ./ ((1-vpress(S,TK)) * molfract);
    case 'O2'
        molfract = 2.0946e-1;
        gas_equil_conc_umolkg = exp(5.80871 + 3.20291*TS + 4.17887*TS.^2 + 5.10006*TS.^3 - 9.86643e-2*TS.^4 + 3.80369*TS.^5 + S.*(-7.01577e-3 - 7.70028e-3*TS - 1.13864e-2*TS.^2 - 9.51519e-3*TS.^3) - 2.75915e-7*S.^2);
        henrys_law_const = gas_equil_conc_umolkg * 1e-6 ./ ((1-vpress(S,TK)) * molfract);
    case 'Ar'
        molfract = 9.34e-3;
        gas_equil_conc_umolkg = exp(2.79150 + 3.17609*TS + 4.13116*TS.^2 + 4.90379*TS.^3 + S.*(-6.96233e-3 - 7.66670e-3*TS - 1.16888e-2*TS.^2));
        henrys_law_const = gas_equil_conc_umolkg * 1e-6 ./ ((1-vpress(S,TK)) * molfract);
    case 'CO2'
        henrys_law_const = exp(-60.2409 + 9345.17./TK + 23.3585*log(TK/100) + S.*(0.023517 - 0.00023656*TK + 0.0047036*(TK/100).^2));
    case 'Ne'
        molfract = 1.818e-5;
        gas_equil_conc_nmolkg = exp(2.18156 + 1.29108*TS + 2.12504*TS.^2 + S.*(-5.94737e-3 - 5.13896e-3*TS));
        henrys_law_const = gas_equil_conc_nmolkg * 1e-9 ./ ((1-vpress(S,TK)) * molfract);
    case 'He'
        molfract = 5.24e-6;
        gas_equil_conc_molkg = exp(-178.1424 + 217.5991*(100./TK) + 140.7506*log(TK/100) - 23.01954*(TK/100) + S.*(-0.038129 + 0.019190*(TK/100) - 0.0026898*(TK/100).^2) -2.55e-6*S.^2);
        henrys_law_const = gas_equil_conc_molkg ./ ((1-vpress(S,TK)) * molfract);
    case 'CH4'
        gas_equil_conc_nmolkg_1atm = exp(-417.5053 + 599.8626*(100./TK) + 380.3636*log(TK/100) - 62.0764*(TK/100)+ S.*(-0.064236 + 0.034980*(TK/100) - 0.0052732*(TK/100).^2));
        henrys_law_const = gas_equil_conc_nmolkg_1atm * 1e-9 ./ (1-vpress(S,TK));
    case 'Kr'
        molfract = 1.14e-6;
        gas_equil_conc_umolkg = exp(-112.6840 + 153.5817*(100./TK) + 74.4690*log(TK/100) - 10.0189*(TK/100) + S.*(-0.011213 - 0.001844*(TK/100) + 0.0011201*(TK/100).^2))/0.0223518;
        henrys_law_const = gas_equil_conc_umolkg * 1e-6 ./ ((1-vpress(S,TK)) * molfract);
    case 'N2O'
        gas_equil_conc_molkg_1atm = exp(-168.2459 + 226.0894*(100./TK) + 93.2817*log(TK/100) - 1.48693*(TK/100).^2+ S.*(-0.060361 + 0.033765*(TK/100) - 0.0051862*(TK/100).^2));
        henrys_law_const = gas_equil_conc_molkg_1atm ./ (1-vpress(S,TK));
    case 'Xe'
        molfract = 8.7e-8;
        gas_equil_conc_molkg = exp(-224.5100 + 292.8234*(100./TK) + 157.6127*log(TK/100) - 22.66895*(TK/100) + S.*(-0.084915 + 0.047996*(TK/100) - 0.0073595*(TK/100).^2) +6.69e-6*S.^2);
        henrys_law_const = gas_equil_conc_molkg ./ ((1-vpress(S,TK)) * molfract);
    case 'CFC-12'
        gas_equil_conc_molkg_1atm = exp(-220.2120 + 301.8695*(100./TK)+114.8533*log(TK/100)-1.39165*(TK/100).^2 + S.*(-0.147718 + 0.093175.*(TK/100)-0.0157340.*(TK/100).^2));
        henrys_law_const = gas_equil_conc_molkg_1atm ./ (1-vpress(S,TK));
    case 'CFC-11'
        gas_equil_conc_molkg_1atm = exp(-232.0411 + 322.5546*(100./TK) + 120.4956*log(TK/100) - 1.39165*(TK/100).^2 + S.*(-0.146531 + 0.093621*(TK/100) - 0.0160693*(TK/100).^2));
        henrys_law_const = gas_equil_conc_molkg_1atm ./ (1-vpress(S,TK));
    case 'SF6'
        gas_equil_conc_molkg_1atm = exp(-82.1639 + 120.152*(100./TK) + 30.6372*log(TK/100) + S.*(0.0293201 - 0.0351974*(TK/100) + 0.00740056*(TK/100).^2));
        henrys_law_const = gas_equil_conc_molkg_1atm ./ (1-vpress(S,TK));
    case 'DMS'
        henrys_law_0sal = exp(3463 ./ TK - 12.20);
        henrys_law_345sal = exp(3547 ./ TK - 12.64);
        henrys_law_const = ((henrys_law_345sal - henrys_law_0sal) / 34.5 .* S + henrys_law_0sal) ./ denscalc(S,T) * 1000;
end %switch
end %main function

function vapor_press_atm = vpress(S,TK)
% Subfunction to calculate water vapour pressure called by primary function if needed

% calculate scaled temperatures
temp_mod = 1-TK./647.096;

% Calculate value of Wagner and Pruss polynomial
Wagner = -7.85951783*temp_mod +1.84408259*temp_mod.^1.5 -11.7866497*temp_mod.^3 +22.6807411*temp_mod.^3.5 -15.9618719*temp_mod.^4 +1.80122502*temp_mod.^7.5;

% Vapor pressure of pure water in atm
vapor_0sal_atm = exp(Wagner * 647.096 ./ TK) .* 217.75;

% Correct vapor pressure for salinity
molality = 31.998 * S ./(1e3-1.005*S);
osmotic_coef = 0.90799 -0.08992*(0.5*molality) +0.18458*(0.5*molality).^2 -0.07395*(0.5*molality).^3 -0.00221*(0.5*molality).^4;
vapor_press_atm = vapor_0sal_atm .* exp(-0.018 * osmotic_coef .* molality);
end %local function

function density = denscalc(S,T)
% Subfunction to calculate density called by primary function if needed

% Convert temperature on ITS-90 scale to IPTS-68 scale for use with density
% equations
T = 1.00024 * T;
  
% Calculate density of pure water
a0 = 999.842594;
a1 =   6.793952e-2;
a2 =  -9.095290e-3;
a3 =   1.001685e-4;
a4 =  -1.120083e-6;
a5 =   6.536332e-9;

dens_0sal = a0 + a1*T + a2*T.^2 + a3*T.^3 + a4*T.^4 + a5*T.^5;

% Correct density for salinity
b0 =  8.24493e-1;
b1 = -4.0899e-3;
b2 =  7.6438e-5;
b3 = -8.2467e-7;
b4 =  5.3875e-9;

c0 = -5.72466e-3;
c1 = +1.0227e-4;
c2 = -1.6546e-6;

d0 = 4.8314e-4;

density = dens_0sal + (b0 + b1*T + b2*T.^2 + b3*T.^3 + b4*T.^4).*S + (c0 + c1*T + c2*T.^2).*S.^1.5 + d0*S.^2;
end %local function