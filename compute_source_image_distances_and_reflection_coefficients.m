function [distances,coefficients] = compute_source_image_distances_and_reflection_coefficients(r_source,r_receiver,Lx,Ly,Lz,c,beta_wall,beta_surface,cutoff_time)

    % The form of the summation used here is derived in: 
    % Allen, J. B., & Berkley, D. A. (1979). Image method for efficiently 
    % simulating small‐room acoustics. The Journal of the Acoustical 
    % Society of America, 65(4), 943-950.
    % (Equation 10, though the sum is done here in the frequency domain
    % instead of the time domain.)
    
    % Inputs:
    % r_source: Vector position of the source (m) [3x1]
    % r_receiver: Vector position of the receiver (m) [3x1]
    % Lx, Ly, Lz: Dimensions of the tank in each coordinate (m)
    % c: Sound speed (m/s)
    % beta_wall: Reflection coefficient for the 5 non-surface walls of the
    % tank
    % beta_surface: Reflection coefficient for the water surface
    % cutoff_time: Time over which to sum reflected paths (s) 
    
    % Outputs:
    % distances: Vector containing the apparent distance of each source 
    % image from the receiver, sorted in ascending order (m)
    % coefficients: Vector containing the product of the reflection
    % coefficients for each source image, whose corresponding distance is
    % given by the output distances (unitless)
    
    % Source and receiver positions are specified in meters in a coordinate
    % system in which the origin lies at one of the vertices of the tank,
    % and x, y, and z are all increasing into the tank, such that the tank
    % walls lie at x=0, x=Lx, y=0, y=Ly, z=0, and the surface is at z=Lz.
    
    % Every reflected path which arrives before the cutoff time is included
    % in the summation, with time t=0 corresponding to the start of the 
    % signal production at the source. Note that some paths which arrive 
    % after the cutoff time will also be included in the summation, so the
    % output is no longer exact after the cutoff time. Ideally, the cutoff
    % time should be chosen to be sufficiently long such that reflected
    % paths which take longer than this to arrive have been sufficiently
    % dampened as to not contribute meaningfully to the signal measured at
    % the receiver.
    
    % Written by Hayden Johnson, 2024-03-25
    % Modified 2024-06-11 to change index convention

    %----------------------------------------------------------------------

    % compute limits of sum from cutoff time
    cutoff_distance = cutoff_time*c;
    l_max = ceil(cutoff_distance./(Lx*2));
    m_max = ceil(cutoff_distance./(Ly*2));
    n_max = ceil(cutoff_distance./(Lz*2));
    max_diagonal_length = sqrt(sum([Lx; Ly; Lz].^2));
    
    % initialize arrays
    distances_unsorted = [];
    coefficients_unsorted = [];
    
    % iterate over lattice displacement vectors
    % (if the Matlab parallel computing toolbox is installed, the outermost
    % loop can be replaced with a parfor loop to get faster performance)
    parfor l = -l_max:l_max 
        for m = -m_max:m_max
            for n = -n_max:n_max
                % compute lattice displacement vector and check whether
                % this block can contain images within the cutoff distance
                r_translation = 2*[l*Lx; m*Ly; n*Lz];
                if sqrt(sum(r_translation.^2)) - 2*max_diagonal_length <= cutoff_distance
                    % iterate over the 8 source images within this block of
                    % the lattice
                    for i = 0:1
                        for j = 0:1
                            for k = 0:1
                                % compute the source image separation
                                % vector within the block
                                r_image = r_translation + [1-2*i; 1-2*j; 1-2*k].*r_source;
                                
                                % compute distance
                                R = r_image - r_receiver;
                                distances_unsorted = [distances_unsorted; sqrt(sum(R.^2))];
                                
                                % compute reflection coefficient product
                                b = beta_wall.^(abs(l-i)+abs(l)+abs(m-j)+abs(m)+abs(n-k)).*beta_surface.^abs(n);
                                coefficients_unsorted = [coefficients_unsorted; b];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % sort arrays by distance of source image
    [distances,ind] = sort(distances_unsorted,1);
    coefficients = coefficients_unsorted(ind,:);
end