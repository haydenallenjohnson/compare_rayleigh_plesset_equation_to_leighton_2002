function [natural_frequency,natural_angular_frequency] = compute_bubble_natural_frequency(R_eq,p_inf,kappa,sigma,rho)

    % Compute the natural frequency of a bubble undergoing spherical
    % oscillations.

    % Inputs:
    % R_eq: Bubble equilibirum radius (m)
    % p_inf: Water pressure far from the bubble (Pa)
    % kappa: Ratio of heat capacities (7/5 for air)
    % sigma: Surface tension of water (N/m)

    % Outputs:
    % natural_frequency: Natural frequency of bubble (Hz)

    %----------------------------------------------------------------------

    natural_angular_frequency = sqrt((3*kappa.*p_inf + (3*kappa-1)*2*sigma./R_eq)./rho)./R_eq;
    natural_frequency = natural_angular_frequency/(2*pi);
end