function gas_diffusion_coef = eh_gasdiff(S,T,gas)

% eh_gasdiff   Molecular diffusion coefficient of a gas in seawater
%=========================================================================
% 
% USAGE:  gas_diffusion_coef = eh_gasdiff(S,T,gas)
%
% DESCRIPTION:
%    Calculates the molecular diffusion coefficient of a gas in seawater at 
%    the given salinity and temperature of the water and a hydrostatic 
%    pressure of 0 dbar (surface)
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%   gas = gas for which to calculate diffusion coefficient, possible values:
%       'N2','O2','Ar','CO2','Ne','He','CH4','Kr','Xe','CFC-12','CFC-11','SF6','DMS','Rn'
%
% OUTPUT:
%   gas_diffusion_coef = molecular diffusion coefficient of the gas 
%       [cm^2 s^-1]
%
% REFERENCE:
%   Appendix E, Chemical Oceanography: Element Fluxes in the Sea
%       Emerson, S.R., R.C. Hamme 2022
%       Cambridge University Press
%   and references therein
%
% VERSION 1.0 : 17 March 2022
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
expectedGases = {'N2','O2','Ar','CO2','Ne','He','CH4','Kr','Xe','CFC-12','CFC-11','SF6','DMS','Rn'};
if sum(strcmp(gas,expectedGases))==0
    error('Expected input parameter gas to match one of these values:\n''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', or ''%s''','N2','O2','Ar','CO2','Ne','He','CH4','Kr','Xe','CFC-12','CFC-11','SF6','DMS','Rn');
end %if

%------
% BEGIN
%------

% Set constants for calculation
switch gas
    case 'N2'
        A = 3412e-5;
        Ea = 18.50e3;
    case 'O2'
        A = 4286e-5;
        Ea = 18.70e3;
    case 'Ar'
        A = 2227e-5;
        Ea = 16.68e3;
    case 'CO2'
        A = 5019e-5;
        Ea = 19.51e3;
    case 'Ne'
        A = 1608e-5;
        Ea = 14.84e3;
    case 'He'
        A = 818e-5;
        Ea = 11.70e3;
    case 'CH4'
        A = 3047e-5;
        Ea = 18.36e3;
    case 'Kr'
        A = 6393e-5;
        Ea = 20.20e3;
    case 'Xe'
        A = 9007e-5;
        Ea = 21.61e3;
    case 'CFC-12'
        A = 3600e-5;
        Ea = 20.1e3;
    case 'CFC-11'
        A = 1500e-5;
        Ea = 18.1e3;
    case 'SF6'
        A = 2900e-5;
        Ea = 19.3e3;
    case 'DMS'
        A = 2000e-5;
        Ea = 18.1e3;
    case 'Rn'
        A = 15877e-5;
        Ea = 23.26e3;
end %switch

% Molecular diffusion coefficient of gas in pure water
TK = T+273.15;
gas_diffusion_coef_0sal = A*exp(-Ea./(8.314462*TK));

% Correct molecular diffusion coefficient for salinity
gas_diffusion_coef = gas_diffusion_coef_0sal .* (1-0.049*S/35.5);

