function dissociation_const = eh_dissconst(S,T,P,reaction)

% eh_dissconst   Equilibrium constant for the specified dissociation 
%             reaction in seawater
%=========================================================================
% 
% USAGE:  gas_equil_conc = eh_dissconst(S,T,P,reaction)
%
% DESCRIPTION:
%    Calculates the equilibrium constants for the dissociation reactions of
%    the carbonate and borate buffer systems, and the solubility products 
%    for aragonite and calcite, at the given salinity, temperature, and
%    hydrostatic pressure of the water.  These constants are based on the 
%    "total" pH scale.
%
% INPUT:  (if S, T, and P are not singular they must have same dimensions)
%   S = practical salinity    [PSS-78 scale]
%   T = temperature [degree C, ITS-90 scale]
%   P = pressure [dbar (decibars), use 0 for surface]
%       Note that pressure in dbar is ~1-2% greater than depth in m.
%   reaction = reaction for which to calculate the dissociation constant
%       use the following to specify:
%       'K1' = first dissociation constant for carbonic acid (CO2 + H2O <--> HCO3- + H+)
%       'K2' = second dissociation constant for carbonic acid (HCO3- <--> CO3-- + H+)
%       'KB' = dissociation constant for boric acid (B(OH)3 + H2O <--> B(OH)4- + H+)
%       'KW' = dissociation constant for water (H2O <--> OH- + H+)
%       'Kcalcite' = solubility product for calcite (CaCO3 <--> Ca++ + CO3--)
%       'Karagonite' = solubility product for aragonite (CaCO3 <--> Ca++ + CO3--)
%
% OUTPUT:
%   dissociation_const = equilibrium constant or solubility product for the 
%       dissociation reaction 
%       [mol kg^-1 for K1, K2, KB; mol^2 kg^-2 for KW, Kcalcite, Karagonite]
%
% REFERENCE:
%   Appendices F & G, Chemical Oceanography: Element Fluxes in the Sea
%       Emerson, S.R., R.C. Hamme 2022
%       Cambridge University Press
%   and references therein
%
% VERSION 1.0 : 25 March 2022
% AUTHOR: Roberta C. Hamme (University of Victoria) 
% This software is available from http://www.cambridge.org/emerson-hamme
% as part of Chemical Oceanography: Element Fluxes in the Sea (2022) 
% by Steven R. Emerson and Roberta C. Hamme
%=========================================================================

%----------------------
% Check input parameters
%----------------------

% check number of input parameters is correct
if nargin ~= 4
   error('Must pass 4 input parameters')
end %if

% check datatypes of input parameters are correct
validateattributes(S,{'numeric'},{'nonempty'},mfilename,'S:salinity')
validateattributes(T,{'numeric'},{'nonempty'},mfilename,'T:temperature')
validateattributes(P,{'numeric'},{'nonempty'},mfilename,'P:pressure')
validateattributes(reaction,{'char'},{'nonempty'},mfilename,'reaction')


% check S,T,&P dimensions and verify they have the same shape or are singular
[rs,cs] = size(S);
[rt,ct] = size(T);
[rp,cp] = size(P);
if (((rs~=rt) || (cs~=ct)) && (rs+cs>2) && (rt+ct>2)) || (((rs~=rp) || (cs~=cp)) && (rs+cs>2) && (rp+cp>2)) || (((rt~=rp) || (ct~=cp)) && (rt+ct>2) && (rp+cp>2))
   error('S, T, and P must have same dimensions or be singular')
end %if

% check that reaction is one of the supported values
expectedReactions = {'K1','K2','KB','KW','Kcalcite','Karagonite'};
if sum(strcmp(reaction,expectedReactions))==0
    error('Expected input parameter reaction to match one of these values:\n''%s'', ''%s'', ''%s'', ''%s'', ''%s'', or ''%s''','K1','K2','KB','KW','Kcalcite','Karagonite');
end %if

%------
% BEGIN
%------

% calculate temperature in Kelvin
TK = T+273.15;
% gas constant (cm^3 bar^-1 mol^-1 K^-1)
R = 83.1446;
% convert pressure from dbar to bar
P = 0.1*P;

% perform calculation for specified reaction
switch reaction
    case 'K1'
        dissociation_const_0dbar = 10.^(-(3633.86./TK - 61.2172 + 9.6777*log(TK) - 0.011555*S + 0.0001152*S.^2));
        dissociation_const = dissociation_const_0dbar .* exp((-(-25.5 + 0.1271*T) + 0.5*(-3.08e-3 + 8.77e-5*T).*P).*P./(R*TK));
    case 'K2'
        dissociation_const_0dbar = 10.^(-(471.78./TK + 25.929 - 3.16967*log(TK) - 0.01781*S + 0.0001122*S.^2));
        dissociation_const = dissociation_const_0dbar .* exp((-(-15.82 - 0.0219*T) + 0.5*(1.13e-3 + -1.475e-4*T).*P).*P./(R*TK));
    case 'KB'
        dissociation_const_0dbar = exp((-8966.9 - 2890.53*S.^0.5 - 77.942*S + 1.728*S.^1.5 - 0.0996*S.^2)./TK + 148.0248 + 137.1942*S.^0.5 + 1.62142*S - (24.4344 + 25.085*S.^0.5 + 0.2474*S).*log(TK) + 0.053105*S.^0.5.*TK);
        dissociation_const = dissociation_const_0dbar .* exp((-(-29.48 + 0.1622*T - 2.608e-3*T.^2) + 0.5*(-2.84e-3).*P).*P./(R*TK));
    case 'KW'
        dissociation_const_0dbar = exp(148.9652 - 13847.26./TK - 23.6521*log(TK) + (118.67./TK - 5.977 + 1.0495.*log(TK)).*S.^0.5 - 0.01615*S);
        dissociation_const = dissociation_const_0dbar .* exp((-(-25.60 + 0.2324*T - 3.6246e-3*T.^2) + 0.5*(-5.13e-3 + 7.94e-5*T).*P).*P./(R*TK));
    case 'Kcalcite'
        dissociation_const_0dbar = 10.^(-(171.9065 + 0.077993*TK - 2839.319./TK - 71.595*log10(TK) + (0.77712 - 0.0028426*TK - 178.34./TK).*S.^0.5 + 0.07711*S - 0.0041249*S.^1.5));
        dissociation_const = dissociation_const_0dbar .* exp((-(-48.76 + 0.5304*T) + 0.5*(-1.176e-2 + 3.692e-4*T).*P).*P./(R*TK));
    case 'Karagonite'
        dissociation_const_0dbar = 10.^(-(171.945 + 0.077993*TK - 2903.293./TK - 71.595*log10(TK) + (0.068393 - 0.0017276*TK - 88.135./TK).*S.^0.5 + 0.10018*S - 0.0059415*S.^1.5));
        dissociation_const = dissociation_const_0dbar .* exp((-(-45.96 + 0.5304*T) + 0.5*(-1.176e-2 + 3.692e-4*T).*P).*P./(R*TK));
end %switch