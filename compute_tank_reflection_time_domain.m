function [y_receiver] = compute_tank_reflection_time_domain(t,y_source,c,image_distances,image_coefficients)

    % Function to compute received signal in a tank based on reflection
    % from walls computed in the time domain. The inputs image_distances 
    % and image_coefficients need to be computed first using the function
    % compute_source_image_distances_and_coefficients().

    % Inputs:
    % y_source: Source time series, can be multiple column vectors (NxM)
    % t: Time vector (Nx1)
    % image_distances: Array specifying distances of images to be
    % considered (Lx1)
    % image_coefficients: Array specifying the combined reflection
    % coefficient that multiplies the arrival from each image (Lx1)

    % Outputs:
    % y_receiver: Signal at the receiver location (NxM)

    % Written by Hayden Johnson, 2024-06-11

    y_receiver = zeros(size(y_source));

    % compute time step
    dt = t(2)-t(1);

    % compute time offset for each image
    t_offset = image_distances./c;
    t_offset_ind = round(t_offset./dt);

    num_relevant_images = find(t_offset_ind < length(y_source),1,'last');
    
    for i = 1:num_relevant_images
        y_image = circshift(y_source,t_offset_ind(i),1).*image_coefficients(i)./image_distances(i);
        y_image(1:t_offset_ind(i)-1,:) = 0;
        y_receiver = y_receiver + y_image;
    end
end