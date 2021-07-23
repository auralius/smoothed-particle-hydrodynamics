function SPHDemo2D_Ex1_GPU()
clear
close all
clc

% Simulation parameters
tEnd            = 30;          % time at which simulation ends
dt              = 0.01;        % timestep
m               = 0.1;         % particle mass
h               = 0.01;        % smoothing length
rho_to_p_const  = 0.1;         % equation of state constant
n_poly          = 1.0;         % polytropic index
nu              = 10.0;        % damping
k_wall          = 100.0;

% Create particles and define their initial locations
k = 1;
for px = 0 : 2*h : 0.2
    for py = 0 : 2*h : 1
        x(k,:)=[px py];
        k = k + 1;     
    end
end

x = gpuArray(double(x));

N = size(x,1);            % Number of particles
a = zeros(N,2,'double','gpuArray');
v_mhalf = zeros(N,2,'double','gpuArray');

K = ceil(tEnd/dt);
DATA = zeros(N,2,K,'double','gpuArray');      % log all position data during the simulation

tic 
for i = 1 : K
    % The integrator
    v_phalf = v_mhalf + a * dt;
    x = x + v_phalf * dt;
    v = 0.5 * (v_mhalf + v_phalf);
    v_mhalf = v_phalf;
    
    X = meshgrid(x(:,1));
    Y = meshgrid(x(:,2));
      
    % Density
    rho_ij = arrayfun(@calculate_rho_ij, X, Y, X', Y', m.*ones(N,N), h.*ones(N,N));
    rho = sum(triu(rho_ij), 2);
    
    % Pressure
    p = rho_to_p_const * rho * (1 + 1/n_poly);
    
    P = meshgrid(p);
    RHO = meshgrid(rho);    
    
    % Acceleration
    a = (- nu .* v) + repmat([0 -9.8*m], N, 1);
    
    % Pressure-added acceleration
    [pax_ij, pay_ij] = arrayfun(@calculate_pa_ij, X, Y, X', Y', RHO, RHO', P, P', m.*ones(N,N),  h.*ones(N,N));
     
    pax = sum(pax_ij);
    pay = sum(pay_ij);
    a = a + [pax' pay'];
    
    % apply contact force by the walls
    [fx, fy] = arrayfun(@calculate_contact_force, x(:,1), x(:,2), k_wall*ones(N,1));
    a = a + ([fx fy] ./ m);
    
    DATA(:,:,i) = x;
end

toc

local_data = gather(DATA);
save('sph_demo1_gpu.mat', 'local_data');
GenerateGIF('sph_demo1_gpu.gif', local_data, 0.1, dt)

end

%%
% Calculate density of two nearby particles
%   see Eq. 15
function rho_ij = calculate_rho_ij(xi, yi, xj, yj, m, h)
r_ij = sqrt((xi-xj)^2 + (yi-yj)^2);
rho_ij = m * Kernel(r_ij, h);
end

%%
% Calculate acceleration between two particles due to pressure
%   see Eq. 14
function [pax_ij, pay_ij] = calculate_pa_ij(xi, yi, xj, yj, RHO_i, RHO_j, P_i, P_j, m, h)
x = xi - xj;
y = yi - yj;
[dWx, dWy] = GradKernel(x, y, h);
z = -m * (P_i/RHO_i^2 + P_j/RHO_j^2);
pax_ij = z * dWx;
pay_ij = z * dWy;
end

%%
% Kernel function in 2D world
%   see Eq. 6
function w = Kernel(r_ij,h)
	w = 1.0 / (h^2*pi)* exp( -r_ij^2 / h^2);	
end

%%
% The grtadient of the Kernel function 
function [dWx, dWy]= GradKernel(x, y, h)
    r_ij = sqrt(x^2 + y^2);
	n = -2 * exp( -r_ij^2 / h^2) / (h^4*pi) ;
	dWx = n * x;
    dWy = n * y;
end

%%
function [fx, fy] = calculate_contact_force(x, y, k_wall)

fx = 0.0;
fy = 0.0;

if x < 0.0 
    fx = -k_wall*x;    % x-potivie force
elseif x > 1.0 
    fx = k_wall*(1.0-x); % x-negative force
end

if y < 0 
    fy = -k_wall*y;    % y-positive force    
end

end

%%
function GenerateGIF(fn, X, every_n_secs, data_sampling_time)
% Initial drawing
h_fig = figure;
hold on 
h_plot = plot(X(:,1,1), X(:,2,1), 'b.');
axis equal
xlim([0 1]);
ylim([0 1]);

K = size(X , 3);
for k = 1 : K
    if mod(k-1, every_n_secs/data_sampling_time) == 0
      set(h_plot,'XData', X(:,1,k), 'YData', X(:,2,k))
      drawnow
      write2gif(h_fig, k, fn);
    end
end
end
