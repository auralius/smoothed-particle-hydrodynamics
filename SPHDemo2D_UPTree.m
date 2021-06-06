function SPHDemo2D_UPTree()
clear all
close all
clc

% Simulation parameters
tEnd            = 50;        % time at which simulation ends
dt              = 0.01;      % timestep
m               = 0.1;       % particle mass
h               = 0.005;      % smoothing length
rho_to_P_const  = 0.1;       % equation of state constant
n_poly          = 1;         % polytropic index
nu              = 10;        % damping
k_wall          = 100;

% load the particles
x = dlmread('up_tree_logo.mat');

N = size(x,1); % Number of particles
a = zeros(N,2);
v_mhalf = zeros(N,2);

K = ceil(tEnd/dt);
X = zeros(N,2,K);      % log all position data during the simulation

tic 

for i =1 : K
    v_phalf = v_mhalf + a * dt;
    x = x + v_phalf * dt;
    v = 0.5 * (v_mhalf + v_phalf);
    v_mhalf = v_phalf;
    
    % update densities, pressures, accelerations
    rho = CalculateDensity(x, m, h);
    P = rho_to_P_const * rho * (1 + 1/n_poly);
    a = CalculateAcceleration(x, v, m, rho, P, nu, h);
    
    % apply contact force by the walls
    f_contact = CalculateContactForce(x, k_wall);
    a = a + f_contact / m;
    
    X(:,:,i) = x;
end

toc

save('up_tree_sim_result.mat', X);
GenerateGIF('sph_demo_up_tree.gif', X, 0.1, dt)

end

%%
function rho = CalculateDensity(x, m, h)
N = size(x,1);
rho = zeros(N,1);

for i = 1 : N    
    % initialize density with i = j contribution
    rho(i) = m * Kernel(0, h);
    
    for j = i+1 : N 
        % calculate vector between two particles
        uij = x(i,:) - x(j,:);
        rho_ij = m * Kernel(uij, h);
        % add contribution to density
        rho(i) = rho(i) + rho_ij;
        rho(j) = rho(j) + rho_ij;
    end
end
end

%%
function a = CalculateAcceleration(x, v, m, rho, P, nu, h)
% initialize accelerations
N = size(x,1);
a = zeros(N, 2); % add damping and gravity

a = a - nu * v + repmat([0 -9.8*m], N, 1);

% add pressure
for i =1: N
    for j = i+1 : N
        % calculate vector between two particles
        uij = x(i, :) - x(j, :);
        
        % calculate acceleration due to pressure
        p_a = -m * (P(i)/rho(i)^2 + P(j)/rho(j)^2) * GradKernel(uij, h);
        a(i, :) = a(i, :) + p_a;
        a(j, :) = a(j, :) - p_a;
    end
end
end

%%
function f = CalculateContactForce(x, k_wall)
% The walls are located at
% x<0, x>1, and y<0

N = size(x,1);
f = zeros(N,2);

for i = 1 : N
    if x(i,1) < 0 
        f(i,1) = -k_wall*x(i,1);    % x-potivie force
    elseif x(i,1) > 1 
        f(i,1) = k_wall*(1-x(i,1)); % x-negative force
    end
    
    if x(i,2) < 0 
        f(i,2) = -k_wall*x(i,2);    % y-positive force    
    end
end
end
%%
function w = Kernel(r,h)
    % 2 dimensions only
	norm_r = norm(r);
	w = 1.0 / (h^2*pi)* exp(-norm_r^2 / h^2);	
end

%%
function dW = GradKernel(r,h)
    % 2 dimensions only
    norm_r = sqrt(r(1)^2 + r(2)^2);
	n = -2 * exp( -norm_r^2 / h^2) / (h^4*pi) ;
	dW = n .* r;
end
%%
function GenerateGIF(fn, X, every_n_secs, data_sampling_time)
% Initial drawing
h_fig = figure;
hold on 
h_plot = plot(X(:,1,1), X(:,2,1), 'b.');
axis equal
xlim([-0.1 1.1]);
ylim([-0.1 1.1]);

K = size(X , 3);
for k = 1 : K
    if mod(k-1, every_n_secs/data_sampling_time) == 0
      set(h_plot,'XData', X(:,1,k), 'YData', X(:,2,k))
      drawnow
      write2gif(h_fig, k, fn);
    end
end
end
