%This file takes and reads a txt file generated from the network_straight_filaments.m file

%%Parse through the file and save data into matrix
%==========================================================================
% Open the text file for reading
fid = fopen('network_nodes.txt', 'r');    %insert file you want to parse

% Initialize variables
fil_num_data = []; %Matrix to store filament number
first_barbed_data = []; %Matrix to store if first node in filament is barbed end
last_barbed_data = [];  %Matrix to store if last node in filament is barbed end
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

     if contains(line, 'FilamentNumber:')
        % Extract the node coordinates
        fil_num_str = extractAfter(line, 'FilamentNumber:');
        fil_num = str2double(strsplit(fil_num_str, ','));
        
        % Append the node coordinates to the matrix
        fil_num_data = [fil_num_data; fil_num];
     end

     if contains(line, 'FirstIsBardedEnd:')
        % Extract the node coordinates
        first_barbed_str = extractAfter(line, 'FirstIsBardedEnd:');
        first_barbed = str2double(strsplit(first_barbed_str, ','));
        
        % Append the node coordinates to the matrix
        first_barbed_data = [first_barbed_data; first_barbed];
     end

     if contains(line, 'LastIsBarbedEnd:')
        % Extract the node coordinates
        last_barbed_str = extractAfter(line, 'LastIsBarbedEnd:');
        last_barbed = str2double(strsplit(last_barbed_str, ','));
        
        % Append the node coordinates to the matrix
        last_barbed_data = [last_barbed_data; last_barbed];
    end
end

% Close the file
fclose(fid);

%%Go through data to see how many filaments and nodes per filament there are
%==========================================================================
start_index = 1;
end_index = 0;
start_index_matrix = [];
end_index_matrix = [];

for i = 1:length(node_data)
    if(i == 1)
        start_index_matrix = [start_index_matrix; start_index];
    end

    if(node_data(start_index,2)==node_data(i,2))
        end_index = end_index + 1;
    end

    if(i ~= length(node_data) && node_data(i,2) ~= node_data(i+1,2))
        %disp(i)
        % Append the # nodes in filament to the matrix
        end_index_matrix = [end_index_matrix; end_index];
        
        start_index = end_index + 1;
        start_index_matrix = [start_index_matrix; start_index];
    end

    if(i == length(node_data))
        end_index = length(node_data);
        end_index_matrix = [end_index_matrix; end_index];
    end
end

%matrix with number of nodes per filament (# filaments is the length of this matrix)
num_nodes_per_fil_matrix = (end_index_matrix - start_index_matrix) + 1;

%%Now check for filaments that are too close
%==========================================================================
%we will only check distance with first and middle node

%new matrix with updated nodes
updated_node_data = zeros(length(node_data),2);

%%to avoid nodes from overlapping 
distance = 0.15; %CHANGE BY USER

%to check distances
first_id =0;
middle_id = 0;
last_id = 0;
dist1 = 0;
dist2 = 0;
dist3 = 0;
dist1_last = 0;
dist2_last = 0;
dist3_last = 0;


%Loop through the filaments & check the distance and update node positions as necessary
for i = 1:length(num_nodes_per_fil_matrix)
    if (i == 1)
        %keep the data for first node
        ending_id = end_index_matrix(1);
        updated_node_data(1:ending_id, :) = node_data(1:ending_id, :);
    else
        %we will only check distance of the first and last node of filament i with first, middle and last node of other filamnets 
        
        %get the current filament's firt and last node positions
        id = start_index_matrix(i,1);
        id_last = end_index_matrix(i,1);
        x1 = node_data(id,1);
        y1 = node_data(id,2);
        x1_last = node_data(id_last,1);
        y1_last = node_data(id_last,2);

        %look through the other filaments
        for j= 1:length(num_nodes_per_fil_matrix)
            if(i ~= j)
                first_id = start_index_matrix(j,1);
                last_id = end_index_matrix(j,1);
                middle_id = first_id + (ceil((last_id - first_id)/2) + 1);

                %get coordinates of first, middle and last node
                x2_first = node_data(first_id,1);
                y2_first = node_data(first_id,2);
                x2_middle = node_data(middle_id,1);
                y2_middle = node_data(middle_id,2);
                x2_last = node_data(last_id,1);
                y2_last = node_data(last_id,2);

                %calculate distances
                dist1 = sqrt((x1 - x2_first)^2 + (y1 - y2_first)^2);
                dist2 = sqrt((x1 - x2_middle)^2 + (y1 - y2_middle)^2);
                dist3 = sqrt((x1 - x2_last)^2 + (y1 - y2_last)^2);
                dist1_last = sqrt((x1_last - x2_first)^2 + (y1_last - y2_first)^2);
                dist2_last = sqrt((x1_last - x2_middle)^2 + (y1_last - y2_middle)^2);
                dist3_last = sqrt((x1_last - x2_last)^2 + (y1_last - y2_last)^2);

                %calculate y-dist
                y_dist = abs(y1 - y2_first);

                %if any of these distances are too close, adjust current node position
                if((dist1 < distance) || (dist2 < distance) || (dist3 < distance) || (dist1_last < distance) || (dist2_last < distance) || (dist3_last < distance))
                    updated_node_data(id:id_last, :) = node_data(id:id_last, :) + (2*distance);
                elseif(y_dist < distance)
                    updated_node_data(id:id_last, :) = node_data(id:id_last, :) + (3*distance);
                else
                    updated_node_data(id:id_last, :) = node_data(id:id_last, :);
                end
         
            end
        end
    end
end


%%Print new txt file with updated nodes
%==========================================================================
% Open the file for writing
fid = fopen('network_updated.txt', 'w');

for i = 1:length(num_nodes_per_fil_matrix)

    %get the filament # and polarity
    fil_num = fil_num_data(i,:);
    first_barbed = first_barbed_data(i,:);
    last_barbed = last_barbed_data(i,:);

    %index for nodes
    id_start = start_index_matrix(i,1);
    id_end = end_index_matrix(i,1);
      
    % Write variables to the file
    fprintf(fid, 'FilamentNumber:%d\n', fil_num);
    fprintf(fid, 'FirstIsBardedEnd:%d\n', first_barbed);
    fprintf(fid, 'LastIsBarbedEnd:%d\n', last_barbed);
    for j = id_start:id_end
        fprintf(fid, 'Node:%0.2f,%0.2f\n', updated_node_data(j,1), updated_node_data(j,2));
    end
    fprintf(fid, 'End_Filament:\n\n');

    %fil_num = fil_num + 1;

end

% Close the file
fclose(fid);





