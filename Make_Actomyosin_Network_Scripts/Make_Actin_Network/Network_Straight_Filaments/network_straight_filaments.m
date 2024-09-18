%This file will generate a txt file with the information to make up a
%network of straight actin filaments.

rng("shuffle");

%%FOR USER: Change the following according to the network settings
%=====================================================================
%=====================================================================
%number of desired filaments in the network:
num_filaments = 20;   

%number of nodes per filament
num_nodes = 6;            %if all filaments have the same # nodes

random_num_nodes = 1;     %set to 1 if each filament will have a diff/random # nodes and change the max # nodes
min_nodes = 3;            %the min should always be 3
max_nodes = 10;

%distance between nodes:
equidistance = 0.2;

%random number scalar: get random numbers between 0 and the scalar. For ex, rand(1)*5 will give numbers between 0 to 5. 
%This will be used for the node positions
rand_scalar = 4;

%to avoid nodes from overlapping 
distance = 0.15;
%=====================================================================
%=====================================================================


%%For storage purposes:
%=====================================================================
fil_num = 0;
first_barbed = 0;
last_barbed = 0;
initial_node_x = 0;
initial_node_y = 0;
node_x = 0;
node_y = 0;

%to store the seed (first node) of every filament
seed_matrix_x = zeros(num_filaments,1);
seed_matrix_y = zeros(num_filaments,1);

%%Generate the network & Print into a txt file:
%=====================================================================
% Open the file for writing
fid = fopen('network_nodes.txt', 'w');

for i = 1:num_filaments
    % Generate random number for the polarity
    first_barbed = round(rand(1));
    disp(first_barbed);

    if(first_barbed==0)
        last_barbed = 1;
    elseif(first_barbed==1)
        last_barbed = 0;
    end

    %Generate random numbers to seed the actin filament
    initial_x = rand(1)*rand_scalar;
    disp(initial_x);
    initial_node_x = round(initial_x,2);  %to keep 2 decimal places
    disp(initial_node_x);

    initial_y = rand(1)*rand_scalar;
    disp(initial_y);
    initial_node_y = round(initial_y,2);
    disp(initial_node_y);

    %store the seed
    seed_matrix_x(i,1) = initial_node_x;
    seed_matrix_y(i,1) = initial_node_y;
    
%     for k=1:num_filaments
%         if(k==1)
%             seed_matrix_x(k,1) = initial_node_x;
%             seed_matrix_y(k,1) = initial_node_y;
%         else
%             %get current seed
%             x1 = seed_matrix_x(k,1);
%             y1 = seed_matrix_y(k,1);
%             updated_x1 = 0;
%             updated_y1 = 0;
% 
%             %loop through all of the points before and check the distance
%             for j=1:num_filaments
%                 if (k~=j)
%                     %get coord of other point
%                     x2 = seed_matrix_x(j,1);
%                     y2 = seed_matrix_y(j,1);
% 
%                     %caculate the distance 
%                     dist = sqrt((x1 - x2)^2 + (y1 - y2)^2);
% 
%                     %check if this dist is less than 2*given distance
%                     if (dist < 2*distance)
%                         updated_x1 = x1 + 2*distance;
%                         updated_y1 = y1 + 2*distance;
% 
%                         seed_matrix_x(k,1) = updated_x1;
%                         seed_matrix_y(k,1) = updated_y1;                        
%                     end    
%                     
%                 end
%             end
% 
%         end
%     end

    %if each filament can have a different number of nodes, randomly select the # nodes in the filament:
    if(random_num_nodes ==1)
        num_nodes = randi([min_nodes, max_nodes]);
    end

    
    % Write variables to the file
    fprintf(fid, 'FilamentNumber:%d\n', fil_num);
    fprintf(fid, 'FirstIsBardedEnd:%d\n', first_barbed);
    fprintf(fid, 'LastIsBarbedEnd:%d\n', last_barbed);
    for j = 1:num_nodes
        if(j==1)
            %fprintf(fid, 'Node:%d\n', initial_node_x, ',',initial_node_y);
            fprintf(fid, 'Node:%0.2f,%0.2f\n', initial_node_x, initial_node_y);
        else
            node_x = initial_node_x + equidistance*(j-1);
            node_y = initial_node_y;
            fprintf(fid, 'Node:%0.2f,%0.2f\n', node_x, node_y);
        end

    end
    fprintf(fid, 'End_Filament:\n\n');

    fil_num = fil_num + 1;

end

% Close the file
fclose(fid);

%%Plot the network:
%=====================================================================




