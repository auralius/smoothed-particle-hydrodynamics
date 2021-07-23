function SPHDemo2D_Ex2_NeighbourSearch_CPU()
clear all
close all
clc

% Simulation parameters
tEnd            = 30;        % time at which simulation ends
dt              = 0.01;      % timestep
m               = 0.1;       % particle mass
h               = 0.01;      % smoothing length
rho_to_p_const  = 0.1;       % equation of state constant
n_poly          = 1;         % polytropic index
nu              = 10;        % damping
k_wall          = 100;
n_neighbours    = 40;

% Create parrticles and define their initial locations
k = 1;
for px = 0 : 2*h : 0.2
    for py = 0 : 2*h : 1
        x(k,:)=[px py];
        k = k + 1;     
    end
end

N = size(x,1); % Number of particles
a = zeros(N,2);
v_mhalf = zeros(N,2);

K = ceil(tEnd/dt);
X = zeros(N,2,K);      % log all position data during the simulation

tic

for i =1 : K
    % do sorting based on the distance of each particle to the origin
    sorted_idx = GetSortedIndex(x);
    
    v_phalf = v_mhalf + a * dt;
    x = x + v_phalf * dt;
    v = 0.5 * (v_mhalf + v_phalf);
    v_mhalf = v_phalf;
    
    % update densities, pressures, accelerations
    rho = CalculateDensity(x, sorted_idx, m, h, n_neighbours);
    P = rho_to_p_const * rho * (1 + 1/n_poly);
    a = CalculateAcceleration(x, sorted_idx, v, m, rho, P, nu, h, k_wall,n_neighbours);
    
    % store the results for animation purpose
    X(:,:,i) = x;
end

toc

save('sph_demo2_ns_cpu.mat', 'X');
GenerateGIF('sph_demo2_ns_cpu.gif', X, 0.1, dt)

end

%%
function sorted_idx = GetSortedIndex(x)
[~,sorted_idx] = sort(vecnorm(x,2,2));
end

%%
function neighbour_idx = FindNeighbours(i, sorted_idx, n_neighbours)
N = length(sorted_idx);

% find where is the location of particele i in the sorted_idx
i_ = find(sorted_idx, i);

% the neighbours are the adjacent particles within a certain range
neighbour_idx = sorted_idx(max(i_-n_neighbours,1) : min(i_+n_neighbours,N));
neighbour_idx = neighbour_idx(neighbour_idx~=i);
end
%%
function rho = CalculateDensity(x, sorted_idx, m, h, n_neighbours)
N = size(x,1);

rho = zeros(N,1);

for i = 1 : N
    % find neighboring particles
    neighbour_idx = FindNeighbours(i, sorted_idx, n_neighbours);
    
    % initialize density with i = j contribution
    rho(i) = m * Kernel(0, h);
    
    for j = 1 : length(neighbour_idx)
        pair_idx = neighbour_idx(j);
        
        % calculate vector between two particles
        uij = x(i,:) - x(pair_idx,:);
        rho_ij = m * Kernel(uij, h);
        
        % add contribution to density
        rho(i) = rho(i) + rho_ij;
    end
end
end

%%
function a = CalculateAcceleration(x, sorted_idx, v, m, rho, P, nu, h, k_wall, n_neighbours)
N = size(x,1);

% initialize accelerations
a = zeros(N, 2); 

% add damping and gravity
a = a - nu * v + repmat([0 -9.8*m], N, 1);

% add pressure
for i =1: N
    curr_idx = sorted_idx(i);
    neighbour_idx = sorted_idx(max(i-n_neighbours,1) : min(i+n_neighbours,N));
    
    for j = 1 : length(neighbour_idx)
        pair_idx = neighbour_idx(j);
        
        % calculate vector between two particles
        uij = x(curr_idx, :) - x(pair_idx, :);
        
        % calculate acceleration due to pressure
        p_a = -m * (P(curr_idx)/rho(curr_idx)^2 + P(pair_idx) / rho(pair_idx)^2) ...
            * GradKernel(uij, h);
        
        a(curr_idx, :) = a(curr_idx, :) + p_a;
    end
    
    f = CalculateContactForceOneParticle(x(curr_idx,:), k_wall);
    a(curr_idx, :) = a(curr_idx, :) + f/m; 
end
end

%%
function f = CalculateContactForceOneParticle(x, k_wall)
% The walls are located at: x<0, x>1, and y<0
% x is for one particel only, thus, it is a 1x2 matrix

f = [0 0];

if x(1,1) < 0
    f(1,1) = -k_wall*x(1,1);    % x-potivie force
elseif x(1,1) > 1
    f(1,1) = k_wall*(1-x(1,1)); % x-negative force
end

if x(1,2) < 0
    f(1,2) = -k_wall*x(1,2);    % y-positive force
end

end

%%
function w = Kernel(r,h)
% 2 dimensions only
norm_r = norm(r);
w = 1.0 / (h^2*pi)* exp( -norm_r^2 / h^2);
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
