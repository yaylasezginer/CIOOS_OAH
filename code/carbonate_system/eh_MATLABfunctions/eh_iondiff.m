function ion_diffusion_coef = eh_iondiff(S,T,ion)

% eh_iondiff   Molecular diffusion coefficient of an ion in seawater
%=========================================================================
% 
% USAGE:  ion_diffusion_coef = eh_iondiff(S,T,ion)
%
% DESCRIPTION:
%    Calculates the molecular diffusion coefficient of an ion in seawater 
%    at the given salinity and temperature of the water and a hydrostatic 
%    pressure of 0 dbar (surface)
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%   ion = ion for which to calculate diffusion coefficient
%       possible values: 'H','K','Na','Ca','Mg','OH','Cl','SO4','HCO3'
%
% OUTPUT:
%   ion_diffusion_coef = molecular diffusion coefficient of the ion 
%       [cm^2 s^-1]
%
% REFERENCE:
%   Appendix E, Chemical Oceanography: Element Fluxes in the Sea
%       Emerson, S.R., R.C. Hamme 2022
%       Cambridge University Press
%   Handbook of Chemistry and Physics (1992/93) 
%       Lide, D.E. (ed.), V.73. Cleveland:CRC Press.
%
% VERSION 1.0 : 8 March 2022
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
validateattributes(ion,{'char'},{'nonempty'},mfilename,'ion')

% check S,T dimensions and verify they have the same shape or are singular
[rs,cs] = size(S);
[rt,ct] = size(T);
if ((rs~=rt) || (cs~=ct)) && (rs+cs>2) && (rt+ct>2)
   error('S & T must have same dimensions or be singular')
end %if

% check that ion is one of the supported values
expectedIons = {'H','K','Na','Ca','Mg','OH','Cl','SO4','HCO3'};
if sum(strcmp(ion,expectedIons))==0
    error('Expected input parameter ion to match one of these values:\n''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', ''%s'', or ''%s''','H','K','Na','Ca','Mg','OH','Cl','SO4','HCO3');
end %if

%------
% BEGIN
%------

% Set constants for calculation
charge = 1;
switch ion
    case 'H'
        lambda = 349.65e-4;
    case 'K'
        lambda = 73.48e-4;
    case 'Na'
        lambda = 50.08e-4;
    case 'Ca'
        lambda = 59.47e-4;
        charge = 2;
    case 'Mg'
        lambda = 53.0e-4;
        charge = 2;
    case 'OH'
        lambda = 198e-4;
    case 'Cl'
        lambda = 76.31e-4;
    case 'SO4'
        lambda = 80.0e-4;
        charge = 2;
    case 'HCO3'
        lambda = 44.5e-4;
end %switch

% Molecular diffusion coefficient of ions in pure water
temp_K = T+273.15;
ion_diffusion_coef_0sal = 8.314462 * temp_K * lambda ./ charge / 96485.3329^2;

% Correct molecular diffusion coefficient for salinity
ion_diffusion_coef = ion_diffusion_coef_0sal .* (1-0.049*S/35.5);

% Convert from units of m^2 s^-1 to cm^2 s^-1
ion_diffusion_coef = ion_diffusion_coef*1e4;
