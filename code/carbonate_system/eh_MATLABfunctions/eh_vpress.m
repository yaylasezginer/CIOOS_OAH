function vapor_press = eh_vpress(S,T,units)

% eh_vpress    Vapor pressure of seawater
%=========================================================================
% 
% USAGE:  vapor_press = eh_vpress(S,T,units)
%
% DESCRIPTION:
%    Calculates the vapor pressure of seawater at the given salinity and
%    temperature of the water and a hydrostatic pressure of 0 dbar
%    (surface)
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%
% OPTIONAL INPUT:
%   units = desired output units formatted as a character string
%       'atm' is default if no value supplied
%       possible values: 'atm','Pa','hPa','mbar','Torr','mm-Hg','psi'
%
% OUTPUT:
%   vapor_press = vapor pressure of seawater 
%       [units = 3rd input, default is atm if 3rd input not supplied]
%
% REFERENCE:
%   Guide to Best Practices for Ocean CO2 Measurements
%   Dickson, A.G., C.L. Sabine, J.R. Christian (Eds.) 2007
%   PICES Special Publication 3, 191pp.
%   Chapter 5: Physical and thermodynamic data
%       Based on: Wagner, W., A. Pruss (2002) The IAPWS formulation 1995 
%       for the thermodynamic properties of ordinary water substance for 
%       general and scientific use, J. Phs. Chem. Ref. Data, 31, 387-535.
%       AND Millero, F.J. (1974) Seawater as a multicomponent electrolyte 
%       solution, pp.3-80.  In: The Sea, Vol. 5, E.D. Goldberg Ed.
%
% VERSION 1.1 : 8 March 2022
% AUTHOR: Roberta C. Hamme (University of Victoria) 
% This software is available from http://www.cambridge.org/emerson-hamme
% as part of Chemical Oceanography: Element Fluxes in the Sea (2022) 
% by Steven R. Emerson and Roberta C. Hamme
%=========================================================================

%----------------------
% Check input parameters
%----------------------

% check number of input parameters is correct
if (nargin ~= 2) && (nargin ~= 3)
   error('Must pass 2 or 3 input parameters')
end %if
% set default if no units are supplied
if nargin == 2
    units = 'atm';
end %if

% check datatypes of input parameters are correct
validateattributes(S,{'numeric'},{'nonempty'},mfilename,'S:salinity')
validateattributes(T,{'numeric'},{'nonempty'},mfilename,'T:temperature')
validateattributes(units,{'char'},{'nonempty'},mfilename,'units')

% check S,T dimensions and verify they have the same shape or are singular
[rs,cs] = size(S);
[rt,ct] = size(T);
if ((rs~=rt) || (cs~=ct)) && (rs+cs>2) && (rt+ct>2)
   error('S & T must have same dimensions or be singular')
end %if

% check that units is one of the supported values
expectedUnits = {'atm','Pa','hPa','mbar','Torr','mm-Hg','psi'};
if sum(strcmp(units,expectedUnits))==0
    error('Expected input parameter units to match one of these values:\n''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', or ''%s''','atm','Pa','hPa','mbar','Torr','mm-Hg','psi');
end %if

%------
% BEGIN
%------

% calculate scaled temperatures
TK = T+273.15;
Tmod = 1-TK./647.096;

% Calculate value of Wagner and Pruss polynomial
Wagner = -7.85951783*Tmod +1.84408259*Tmod.^1.5 -11.7866497*Tmod.^3 +22.6807411*Tmod.^3.5 -15.9618719*Tmod.^4 +1.80122502*Tmod.^7.5;

% Vapor pressure of pure water in atm
vapor_0sal_atm = exp(Wagner * 647.096 ./ TK) .* 217.75;

% Correct vapor pressure for salinity
molality = 31.998 * S ./(1e3-1.005*S);
osmotic_coef = 0.90799 -0.08992*(0.5*molality) +0.18458*(0.5*molality).^2 -0.07395*(0.5*molality).^3 -0.00221*(0.5*molality).^4;
vapor_press_atm = vapor_0sal_atm .* exp(-0.018 * osmotic_coef .* molality);

% Convert to desired units
switch units
    case 'atm'
        vapor_press = vapor_press_atm;
    case 'Pa'
        vapor_press = vapor_press_atm * 101325;
    case 'hPa'
        vapor_press = vapor_press_atm * 1013.25;
    case 'mbar'
        vapor_press = vapor_press_atm * 1013.25;
    case 'Torr'
        vapor_press = vapor_press_atm * 760;
    case 'mm-Hg'
        vapor_press = vapor_press_atm * 760;
    case 'psi'
        vapor_press = vapor_press_atm * 14.6959;
end %switch
