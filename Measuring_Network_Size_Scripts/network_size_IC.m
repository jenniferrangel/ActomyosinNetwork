%This script is used to calculate the network size of an actomyosin network (particularly the initial configuration)
%This is achieved by drawing a circle around the network and calculating the radius of the circle

%===========================================
%%USER INPUT:
%===========================================
%Load the data from the txt files that contain the node positions for actin filaments and pMyoII mini-filaments:
filament_node_data = load('actin_position.txt');
myosin_node_data = load('myosin_position.txt');

%Number of nodes per filament and per pMyoII mini-filament
num_filament_nodes = 6;
num_myosin_nodes = 2;
%=========================================================================

% Define the number of filaments and pMyoII mini-filament based on number of nodes
num_filaments = size(filament_node_data, 1) / num_filament_nodes; 
num_myosins = size(myosin_node_data, 1) / num_myosin_nodes;

%Combine the node positions
x = [filament_node_data(:, 1); myosin_node_data(:, 1)];
y = [filament_node_data(:, 2); myosin_node_data(:, 2)];

% Calculate the center of the combined network configuration (mean of x-coord and mean of y-coord)
center_x = mean(x);
center_y = mean(y);

% Calculate the radius of the network (maximum distance from the center)
distances = sqrt((x - center_x).^2 + (y - center_y).^2);
radius = max(distances);

%Plot the network with the circle encompassing the network
% Plot the nodes
figure;
hold on;
for i = 1:num_filaments
    filament_nodes = filament_node_data((i-1)*num_filament_nodes + 1:i*num_filament_nodes, :);
    plot(filament_nodes(:,1), filament_nodes(:,2),  '.-','MarkerSize', 20,'Color', [0, 0.5, 0], 'DisplayName', 'Actin Filaments'); % plot filament nodes in green
end

for i = 1:num_myosins
    myosin_nodes = myosin_node_data((i-1)*num_myosin_nodes + 1:i*num_myosin_nodes, :);
    plot(myosin_nodes(:, 1), myosin_nodes(:, 2), 'r.-', 'MarkerSize', 15,'DisplayName', 'pMyoII Mini-Filaments'); % plot myosin nodes in red
end

% Plot the center of the network
plot(center_x, center_y, 'kx', 'MarkerSize', 20, 'LineWidth', 4, 'DisplayName', 'Center'); % plot center in black with an 'x'

% Draw the circle
theta = linspace(0, 2*pi, 100);
circle_x = center_x + radius * cos(theta);
circle_y = center_y + radius * sin(theta);
plot(circle_x, circle_y, 'b-', 'LineWidth', 2, 'DisplayName', 'Initial Disc');

% Set axis equal for proper scaling
axis equal;
%title('Estimating Network Size');
title('Initial Network Configuration');
xlabel('x');
ylabel('y');
legend ('show');

% Display the radius and coordinates of the center
disp(['Radius of the circle: ', num2str(radius)]);
disp(['Center of the circle: (', num2str(center_x), ', ', num2str(center_y), ')']);
