function [b_th, b_ac, b_vs, w0] = BubbleDamping(fvec, R0, z, T, patm, gas)
%
% function [b_th, b_ac, b_vs] = BubbleDamping(fvec, R0, z, T, p0, gas);
%
% Use Prosperetti's 1977 paper to compute compute bubble damping 
% coefficients versus frequency.
%
% R0: bubble radius in meters
% z : water depth in meters (Salinity: 32 ppt, T: 20 Celsius)
% T : air and water temperature in Kelvin (20 Celsius = 293.15 K)
% p0: air pressure at sealevel in Pascal (STP = 101.235e3 pa)
% gas: 2 = oxygen, 3 = nitrogen, 4 = methane, everything else air
%
% b_th : thermal damping in /s
% b_ac : radiation damping in /s
% b_vs : viscous damping in /s
% w0   : natural frequency in radians/s
%
% gbd 9/05/07
%
%% Set up the physical constants
%
g = 9.80665;        % acceleration due to gravity
Rg = 8.314472;      % universal gas constant in Joules K^-1 mol^-1
%
% The thermal conductivity and ratio of specific heats came from:
% http://www.engineeringtoolbox.com/oxygen-d_978.html
%
% The thermal conductivity cam from:
% http://www.webelements.com/webelements/elements/text/O/heat.html
%
% Oxygen
%
Mg_oxygen = 0.0320;     % Molecular weight Oxygen Kg / mol
Kg_oxygen = 0.02658;    % Thermal conductivity Oxygen W /m /K
Cp_oxygen = 915 + (918-915)/(300-275)*(290-275);     % J/Kg/K (290 K)
gamma_oxygen = 1.4;     % ratio of specific hetas
%
% Nitrogen
%
Mg_nitrogen = 0.02802;  % Molecular weight Nitrogen Kg / mol
Kg_nitrogen = 0.02583;  % Thermal conductivity Nitrogen W /m /K
Cp_nitrogen = 1039 + (1040-1039)/(300-275)*(290-275); % J/Kg/K (290 K)
gamma_nitrogen = 1.4;   % ratio of specific hetas
%
% Methane
%
Mg_methane = 0.016;     % Molecular weight methane Kg / mol
Kg_methane = 0.0332;    % Thermal conductivity methane W /m /K (290 K)
Cp_methane = 2191 + (2226-2191)/(300-275)*(290-275);  % J/Kg/K (290 K)
gamma_methane = 1.31;   % ratio of specific hetas
%
% Specific Heat Ratio for air.
% http://www.efunda.com/Materials/common_matl/show_gas.cfm?MatlName=Air0C
%
% Thermal conductivity of air. This is assumed to be independent of
% pressure.
% http://home.worldonline.dk/jsrsw/Tcondvspressure.html
%
Mg_air = 0.02896;       % air molecular weight in kg mol^-1 
Kg_air = 0.0254;        % thermal conductivity air W/m/K
gamma_air = 1.40;       % ratio of specific heats
Cp_air = 1.00e3;        % specific heat capacity of air in J Kg^-1 K^-1
%
% Water Properties. These are assumed independent of pressure
% Seawater density.
% http://www.es.flinders.edu.au/~mattom/Utilities/density.html
% T = 20 Celsius, salinity = 32 ppt, pressure = 0
%
rhow = 1022.476;    % seawater density in kg m^-3
Kw = 0.58;          % thermal conductivity water W/m/K
Dw = 1.45e-7;       % thermal diffusivity of water in m^2 s^-1
mu = 1e-3;          % viscosity of water
sigma = 0.072;      % water surface tension in N m
% sigma = sigma/2
%
% Ambient pressure and internal, equilibrium bubble pressure at depth
%
pamb = patm + z*rhow*g; % ambient water pressure
pineq = pamb + 2*sigma/R0;
%
% Speed of sound in seawater. 
% T = 20 Celsius, Salinity = 32 ppt, Pressure = 0 dBar.
% Speed of sound affects radiation losses.
%
c = 1500;
w = 2*pi*fvec;
k = w/c;
%
% Set up the gas values
%
if gas == 2
    M = Mg_oxygen;
    Kg = Kg_oxygen;
    Cp = Cp_oxygen;
    gamma = gamma_oxygen;
elseif gas == 3
    M =  Mg_nitrogen;
    Kg = Kg_nitrogen;
    Cp = Cp_nitrogen;
    gamma = gamma_nitrogen;
elseif gas == 4
    M = Mg_methane;
    Kg = Kg_methane;
    Cp = Cp_methane;
    gamma = gamma_methane;
else
    M = Mg_air;
    Kg = Kg_air;
    Cp = Cp_air;
    gamma = gamma_air;
end
%
% Compute air density from the internal, equilibrium pressure
%
rhog = M*pineq/(Rg*T);
%
% Compute thermal diffusivity of gas. Note Prosperetti (1977, p.19) has a 
% typo: Cv,g should by Cp,g as used below.
%
Dg = Kg/(rhog*Cp);
%
%% Losses versus frequency for fixed radius and depth.
%
G1 = M*Dg*w/(gamma*Rg*T);
G2 = w*R0^2/Dg;
G3 = w*R0^2/Dw;
f = 1+(1+1i)*sqrt(0.5*G3);
kratio = Kw/Kg;
beta1 = sqrt(0.5*gamma*G2.*(1i-G1+sqrt((1i-G1).^2+4*1i*G1/gamma)));
beta2 = sqrt(0.5*gamma*G2.*(1i-G1-sqrt((1i-G1).^2+4*1i*G1/gamma)));
lambda1 = beta1.*coth(beta1)-1;
lambda2 = beta2.*coth(beta2)-1;
Gamma1 = 1i+G1+sqrt((1i-G1).^2+4*1i*G1/gamma);
Gamma2 = 1i+G1-sqrt((1i-G1).^2+4*1i*G1/gamma);
psi = (kratio.*f.*(Gamma2-Gamma1)+lambda2.*Gamma2-lambda1.*Gamma1)./ ...
    (kratio.*f.*(lambda2.*Gamma1-lambda1.*Gamma2)-lambda1.*lambda2.*(Gamma2-Gamma1));
muth = 1/4*w*rhog*R0^2.*imag(psi);
kappa = 1/3*(w.^2*rhog*R0^2/pineq).*real(psi);
b_th = 2*muth/(rhow*R0.^2);
b_ac = 1/2*w.*(w*R0/c)./(1+(w*R0/c).^2);
b_vs = 2*mu/(rhow*R0.^2)*ones(size(w));
%
%w0 = sqrt( 3*kappa*pineq/(rhow*R0^2) - 2*sigma/(rhow*R0^3) + ...
%     k.^2*R0^2./(1+k.^2*R0^2).*w.^2 );
w0 = sqrt( 3*kappa*pineq/(rhow*R0^2) - 2*sigma/(rhow*R0^3));
return