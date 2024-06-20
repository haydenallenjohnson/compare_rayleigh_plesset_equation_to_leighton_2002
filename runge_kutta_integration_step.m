function [y_out,k1] = runge_kutta_integration_step(f, t, y_in, h)
    % Function to implement 4th-order Runge-Kuuta integration of a function
    % f(t,y).
    
    % Inputs:
    % f: Function handle.
    
    %----------------------------------------------------------------------
    
    k1 = f(t, y_in);
    k2 = f(t+h/2, y_in+h*k1/2);
    k3 = f(t+h/2, y_in+h*k2/2);
    k4 = f(t+h, y_in+h*k3);
    y_out = y_in + h*(k1+2*k2+2*k3+k4)/6;
    
end