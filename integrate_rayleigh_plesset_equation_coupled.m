function [t, R, Rd, Rdd, p_radiated, p_external] = integrate_rayleigh_plesset_equation_coupled(integration_time_step,duration,R_eq,forcing_amplitude,kappa,rho,c,depth,image_distances,image_coefficients)
    
    % Integrate the Rayleigh-Plesset equation in a tank where the reflected
    % field influences the external pressure field at the bubble.
    
    % Inputs:
    % integration_time_step: (s)
    % duration: Time to integrate for (s)
    % R_initial: Initial bubble radius (m)
    

    %----------------------------------------------------------------------
    
    % constants
    g = 9.8; % N/kg
    sigma = 0.072; % surface tension in N/m
    temperature = 273;
    p_atm = 101.3e3;
    gas = 1; % gas number for BubbleDamping()
    p_inf = p_atm + rho.*g.*depth;
    
    % initialize arrays
    t = (0:integration_time_step:duration)';
    R = zeros(size(t));
    Rd = zeros(size(t));
    Rdd = zeros(size(t));
    p_radiated = zeros(size(t)); % radiated acoustic pressure at 1 m
    p_external = zeros(size(t)); % pressure field of reflected radiation felt at the bubble

    % specify initial conditions
    R(1) = R_eq;
    Rd(1) = 0;

    % compute natural frequency
    [natural_frequency,natural_angular_frequency] = compute_bubble_natural_frequency(R_eq,p_inf,kappa,sigma,rho);

    % calculate initial forcing
    forcing_ind = t <= pi./natural_angular_frequency;
    p_external(forcing_ind) = forcing_amplitude.*sin(natural_angular_frequency.*t(forcing_ind)).^2;

    % compute bubble damping
    [b_th, b_ac, b_vs, ~] = BubbleDamping(natural_frequency, R_eq, depth, temperature, p_atm, gas);
    b = b_th+b_ac+b_vs;
    mu_effective = rho.*R_eq.^2.*b/2;

    % compute time offsets for each image
    t_offset = image_distances./c;
    t_offset_ind = round(t_offset./integration_time_step);

    % initialize contributing image counters
    image_counter = 0;

    % iterate through each time step
    for i = 1:length(t)-1
        % count how many images are contributing to the external pressure
        % field at the bubble
        while image_distances(image_counter+1) <= c*t(i)
            image_counter = image_counter + 1;
        end
        
        % sum image contributions to external pressure field at the bubble
        for j = 2:image_counter
            image_time_ind = i - t_offset_ind(j);
            p_external(i) = p_external(i) + p_radiated(image_time_ind)*image_coefficients(j)./image_distances(j);
        end
        
        % integrate bubble radius
        y_in = [R(i); Rd(i)];
        [y_out,k1] = runge_kutta_integration_step(@wrapper_function,t(i),y_in,integration_time_step);
        Rdd(i) = k1(2);
        R(i+1) = y_out(1);
        Rd(i+1) = y_out(2);
        
        % compute radiated acoustic pressure
        p_radiated(i) = rho*R(i).*(Rdd(i).*R(i)+2*Rd(i).^2);
    end
    
    % compute last value of Rdd and p_radiated
    [~,Rdd(end)] = rayleigh_plesset_equation(R(end),Rd(end),R_eq,p_inf,p_external(end),kappa,rho,mu_effective);
    p_radiated(end) = rho*R(end).*(Rdd(end).*R(end)+2*Rd(end).^2);

    function [dy] = wrapper_function(time,y)
        % wrapper function to allow for extra variables to be passed to
        % rayleigh_plesset_equation() while using runge_kutta_integration()
        t_ind = find(t>=time,1,'first');
        [dy(1,:),dy(2,:)] = rayleigh_plesset_equation(y(1),y(2),R_eq,p_inf,p_external(t_ind),kappa,rho,mu_effective);
    end
end