function density = eh_dens0(S,T)

% eh_dens0   Density of seawater at atmospheric pressure
%=========================================================================
% 
% USAGE:  density = eh_dens0(S,T)
%
% DESCRIPTION:
%    Calculates the density of seawater at the given salinity and
%    temperature of the water and a hydrostatic pressure of 0 dbar
%    (surface)
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%
% OUTPUT:
%   density = density of seawater [kg m^-3]
%
% REFERENCE:
%   Unesco 1983. Algorithms for computation of fundamental properties of 
%      seawater, 1983. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.
%
%   Millero, F.J. and  Poisson, A.
%     International one-atmosphere equation of state of seawater.
%     Deep-Sea Res. 1981. Vol28A(6) pp625-629.
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
c1 =  1.0227e-4;
c2 = -1.6546e-6;

d0 = 4.8314e-4;

density = dens_0sal + (b0 + b1*T + b2*T.^2 + b3*T.^3 + b4*T.^4).*S + (c0 + c1*T + c2*T.^2).*S.^1.5 + d0*S.^2;

