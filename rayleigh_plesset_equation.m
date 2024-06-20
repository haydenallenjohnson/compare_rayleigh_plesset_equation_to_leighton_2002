function [Rd_output,Rdd_output] = rayleigh_plesset_equation(R_input,Rd_input,R_eq,p_inf,p_external,kappa,rho,mu)
    
    % Solve one time step of the reduced-order Rayleigh-Plesset equation
    
    %----------------------------------------------------------------------
    
    % set default values of optional paramters
    if ~exist('rho','var')
        rho = 1000;
    end
    
    sigma = 0.072; % surface tension of water in N/m
    
    Rd_output = Rd_input;
    
    p_L = (p_inf + 2*sigma./R_eq).*(R_eq./R_input).^(3*kappa) - 2*sigma./R_input;
    Rdd_output = ((p_L - p_inf - p_external - 4*mu.*Rd_input./R_input)./rho - (3/2)*Rd_input.^2)./R_input;
    
end