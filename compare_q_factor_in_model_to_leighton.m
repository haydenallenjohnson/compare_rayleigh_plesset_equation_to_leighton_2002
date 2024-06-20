%% set up model
% spcify constants
g = 9.8;
c = 1480;
rho = 1000;
kappa = 7/5;
sigma = 0.072;
p_atm = 101.3e3;

% tank reflection coefficients
beta_wall = [-0.9 -0.7];
num_beta = length(beta_wall);
beta_surface = -1;

% bubble radius
R_eq = [0.8 0.75 0.7 0.65 0.60 0.575 0.55 0.525 0.5 0.475 0.45 0.425 0.4 0.38 0.36 0.34 0.32]*1e-3; % m
forcing_amplitude = 1000; % Pa

% define tank size
Lx = 0.6;
Ly = 0.2;
Lz = 0.23;

% specify source location
x_source = 0.15;
y_source = 0.1;
z_source = 0.115;
r_source = [x_source; y_source; z_source];

% specify receiver location
x_receiver = 0.45;
y_receiver = 0.1;
z_receiver = 0.115;
r_receiver = [x_receiver; y_receiver; z_receiver];

depth = Lz - z_source;
p_inf = p_atm + rho*g*depth;

% integration properties
integration_time_step = 1e-6;
duration = 10e-3;
t = (0:integration_time_step:duration)';

% compute bubble natural frequency
[natural_frequency,~] = compute_bubble_natural_frequency(R_eq,p_inf,kappa,sigma,rho);

%% compute image properties
[image_distances_at_source,image_coefficients_at_source] = compute_source_image_distances_and_reflection_coefficients(r_source,r_source,Lx,Ly,Lz,c,beta_wall,beta_surface,duration);
[image_distances_at_receiver,image_coefficients_at_receiver] = compute_source_image_distances_and_reflection_coefficients(r_source,r_receiver,Lx,Ly,Lz,c,beta_wall,beta_surface,duration);

%% compute signal radiated by bubble in free space
p_radiated_free_bubble = zeros(length(t),length(R_eq));
parfor i = 1:length(R_eq)
    [~, ~, ~, ~, p_radiated_free_bubble(:,i), ~] = integrate_rayleigh_plesset_equation_free(integration_time_step,duration,R_eq(i),forcing_amplitude,kappa,rho,depth);
end

%% compute signal radiated by bubble in tank
p_radiated_tank = zeros(length(t),length(R_eq),length(beta_wall));
parfor i = 1:length(R_eq)
    for j = 1:num_beta
        [~, ~, ~, ~, p_radiated_tank(:,i,j), ~] = integrate_rayleigh_plesset_equation_coupled(integration_time_step,duration,R_eq(i),forcing_amplitude,kappa,rho,c,depth,image_distances_at_source,image_coefficients_at_source(:,j));
    end
end

%% compute signal at receiver in tank
p_receiver_tank = zeros(length(t),length(R_eq),length(beta_wall));
parfor j = 1:length(beta_wall)
    p_receiver_tank(:,:,j) = compute_tank_reflection_time_domain(t,p_radiated_tank(:,:,j),c,image_distances_at_receiver,image_coefficients_at_receiver(:,j));
end

%% compute q values
q_free_bubble = zeros(length(R_eq),1);
q_tank_radiated = zeros(length(R_eq),length(beta_wall));
q_tank_receiver = zeros(length(R_eq),length(beta_wall));

for i = 1:length(R_eq)
    [m,ind_1] = max(p_radiated_free_bubble(:,i));
    ind_2 = find(p_radiated_free_bubble(:,i) >= exp(-pi).*m,1,'last');
    q_free_bubble(i) = (t(ind_2)-t(ind_1)).*natural_frequency(i);

    for j = 1:length(beta_wall)
        [m,ind_1] = max(p_radiated_tank(:,i,j));
        ind_2 = find(p_radiated_tank(:,i,j) >= exp(-pi).*m,1,'last');
        q_tank_radiated(i,j) = (t(ind_2)-t(ind_1)).*natural_frequency(i);

        [m,ind_1] = max(p_receiver_tank(:,i,j));
        ind_2 = find(p_receiver_tank(:,i,j) >= exp(-pi).*m,1,'last');
        q_tank_receiver(i,j) = (t(ind_2)-t(ind_1)).*natural_frequency(i);
    end
end

%% create plots
for j = 1:length(beta_wall)
    figure(j);
    scatter(natural_frequency,q_free_bubble,'displayname','Free bubble');
    hold on;
    scatter(natural_frequency,q_tank_radiated(:,j),'displayname','Tank radiated signal');
    scatter(natural_frequency,q_tank_receiver(:,j),'displayname','Tank received signal');
    hold off;
    grid on;
    ylim([0 60]);
    xlabel('Frequency (Hz)');
    ylabel('Q factor');
    title(['\beta = ' num2str(beta_wall(j),1)]);
    legend('location','southeast');
    exportgraphics(gcf,['q_vs_frequency_leighton_tank_beta_' num2str(-10*beta_wall(j),1) '.pdf']);
end