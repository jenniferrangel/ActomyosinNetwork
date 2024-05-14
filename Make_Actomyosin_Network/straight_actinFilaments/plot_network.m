%this is a file that reads in a txt file and plots it

% Open the text file for reading
fid = fopen('network7.txt', 'r');    %insert file you want to parse
%fid = fopen('network3.txt', 'r');    %insert file you want to parse

% Initialize variables
node_data = []; % Matrix to store node coordinates

% Loop through each line in the file
while ~feof(fid)
    line = fgetl(fid); % Read one line
    
    % Check if the line contains 'Node:'
    if contains(line, 'Node:')
        % Extract the node coordinates
        node_coords_str = extractAfter(line, 'Node:');
        node_coords = str2double(strsplit(node_coords_str, ','));
        
        % Append the node coordinates to the matrix
        node_data = [node_data; node_coords];
    end
end

% Close the file
fclose(fid);

%variables needed to only connect points with same y-coord
start_index = 1;
end_index = 0;

% Plot the nodes
figure;
for i = 1:length(node_data)
    if(node_data(start_index,2)==node_data(i,2))
        end_index = end_index + 1;
    end
    if(i ~= length(node_data) && node_data(i,2) ~= node_data(i+1,2))
        %disp(i)
        plot(node_data(start_index:end_index,1), node_data(start_index:end_index,2), 'r.-','MarkerSize', 20, 'LineWidth', 1); 
        hold on;
        start_index = end_index + 1;
    end
    if(i == length(node_data))
        plot(node_data(start_index:length(node_data)), node_data(start_index:length(node_data),2), 'r.-','MarkerSize', 20, 'LineWidth', 1); 
        hold on;
    end
end
xlabel('X');
ylabel('Y');
title('Actin Network');
grid on;