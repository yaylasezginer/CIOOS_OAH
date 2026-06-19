function dyn_viscosity = eh_dynvisc(S,T)

% eh_dynvisc   Dynamic viscosity of seawater
%=========================================================================
% 
% USAGE:  dyn_viscosity = eh_dynvisc(S,T)
%
% DESCRIPTION:
%    Calculates the dynamic viscosity of seawater at the given salinity and
%    temperature of the water and a hydrostatic pressure of 0 dbar
%    (surface)
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%
% OUTPUT:
%   dyn_viscosity = dynamic viscosity of seawater [Pa s = kg m^-1 s^-1]
%
% REFERENCE:
%   Sharqawy, M. H., J. H. Lienhard, and S. M. Zubair (2010) 
%       The thermophysical properties of seawater:
%       a review of existing correlations and data, 
%       Desalination and Water Treatment, 16, 354-380.
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
if (nargin ~= 2)
   error('Must pass 2 input parameters')
end %if

% check datatypes of input parameters are correct
validateattributes(S,{'numeric'},{'nonempty'},mfilename,'S:salinity')
validateattributes(T,{'numeric'},{'nonempty'},mfilename,'T:temperature')

% check S,T dimensions and verify they have the same shape or are singular
[rs,cs] = size(S);
[rt,ct] = size(T);
if ((rs~=rt) || (cs~=ct)) && (rs+cs>2) && (rt+ct>2)
   error('S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------
  
% Dynamic viscosity of pure water in atm
dyn_viscosity_0sal = 4.2844e-5 + (0.157*(T + 64.993).^2 - 91.296).^-1;

% Correct dynamic viscosity for salinity
A = 1.541 + 1.998e-2*T - 9.52e-5*T.^2;
B = 7.974 - 7.561e-2*T + 4.724e-4*T.^2;
dyn_viscosity = dyn_viscosity_0sal .* (1 + A.*S/1000 + B.*(S/1000).^2);

