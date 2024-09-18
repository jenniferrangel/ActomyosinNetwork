%This will generate an network composed of actin filaments and myosin mini-filaments. The filaments are randomly oriented

rng("shuffle");

%%FOR USER: Change the following according to the network settings
%=====================================================================
%=====================================================================
%Number of desired filaments in the network:
num_filaments = 200; 
num_myosins = 20;

%number of nodes per filament
num_nodes = 6;            %if all filaments have the same # nodes
num_myosin_nodes = 2;

random_num_nodes = 0;     %set to 1 if each filament will have a diff/random # nodes and change the max # nodes
min_nodes = 3;            %the min should always be 3
max_nodes = 10;

%distance between nodes:
equidistance = 0.2;
myo_equidistance = 0.15;

%random number scalar: get random numbers between 0 and the scalar. For ex, rand(1)*5 will give numbers between 0 to 5. 
%This will be used for the node positions
rand_scalar = 8;

%=====================================================================
%=====================================================================

%%For storage purposes:
%=====================================================================
fil_num = 0;
first_barbed = 0;
last_barbed = 0;

node_data = zeros(num_filaments*num_nodes,2);
first_barbed_data = zeros(num_filaments,1);
last_barbed_data = zeros(num_filaments,1);

nodes_per_fil = zeros(num_filaments,1);

myo_num = 0;
myosin_node_data = zeros(num_myosins*num_myosin_nodes,2);

%%Generate the network & Print into a txt file:
%=====================================================================
% Open the file for writing
fid = fopen('random_network_200Actin_6nodes_20myo.txt', 'w'); %USER INPUT: name file accordingly

for i = 1:num_filaments
    % Generate random number for the polarity
    first_barbed = round(rand(1));

    if(first_barbed==0)
        last_barbed = 1;
    elseif(first_barbed==1)
        last_barbed = 0;
    end

    first_barbed_data(i,1) = first_barbed;
    last_barbed_data(i,1) = last_barbed;

    %Generate the actin filament
    %Randomly select the seed/starting point of the filament
    x_seed = rand * rand_scalar;
    y_seed = rand * rand_scalar;

    %Random angle for the actin filament
    angle = rand * 2 * pi;

    %if each filament has the same number of nodes:
    if(random_num_nodes == 0)

        % Write variables to the file
        fprintf(fid, 'FilamentNumber:%d\n', fil_num);
        fprintf(fid, 'FirstIsBardedEnd:%d\n', first_barbed);
        fprintf(fid, 'LastIsBarbedEnd:%d\n', last_barbed);
    
        %Calculate the (x,y) positions for each node in the actin filament
        for j = 1:num_nodes
            node_x = x_seed + (equidistance * cos(angle) * (j-1));
            node_y = y_seed + (equidistance * sin(angle) * (j-1));
    
            %store into matrix 
            index = j + (num_nodes*(i-1));
            node_data(index,1) = node_x;
            node_data(index,2) = node_y;
    
            %print
            fprintf(fid, 'Node:%0.2f,%0.2f\n', node_x, node_y);
        end
    
        fprintf(fid, 'End_Filament:\n\n');
    end

    %if each filament can have a different number of nodes, randomly select the # nodes in the filament:
    if(random_num_nodes == 1)
        num_nodes = randi([min_nodes, max_nodes]);
        nodes_per_fil(i,1) = num_nodes;

        % Write variables to the file
        fprintf(fid, 'FilamentNumber:%d\n', fil_num);
        fprintf(fid, 'FirstIsBardedEnd:%d\n', first_barbed);
        fprintf(fid, 'LastIsBarbedEnd:%d\n', last_barbed);
    
        %Calculate the (x,y) positions for each node in the actin filament
        for j = 1:num_nodes
            node_x = x_seed + (equidistance * cos(angle) * (j-1));
            node_y = y_seed + (equidistance * sin(angle) * (j-1));
    
            %store into matrix 
            if (i == 1)
                index = j + (num_nodes*(i-1));
                node_data(index,1) = node_x;
                node_data(index,2) = node_y;
            else
                index = j + (sum(nodes_per_fil(1:(i-1))));
                node_data(index,1) = node_x;
                node_data(index,2) = node_y;
            end
    
            %print
            fprintf(fid, 'Node:%0.2f,%0.2f\n', node_x, node_y);
        end
    
        fprintf(fid, 'End_Filament:\n\n');
    end

    fil_num = fil_num + 1;
end

for i = 1:num_myosins
    %Generate the myosin mini-filaments
    %Randomly select the seed/starting point of the filament
    x_myo_seed = rand * rand_scalar;
    y_myo_seed = rand * rand_scalar;

    %Random angle for the actin filament
    angle_myo = rand * 2 * pi;
   
    % Write variables to the file
    fprintf(fid, 'MyosinNumber:%d\n', myo_num);
    
    %Calculate the (x,y) positions for each myosin node
    for j = 1:num_myosin_nodes
        myo_node_x = x_myo_seed + (myo_equidistance * cos(angle_myo) * (j-1));
        myo_node_y = y_myo_seed + (myo_equidistance * sin(angle_myo) * (j-1));

        %store into matrix 
        index_myo = j + (num_myosin_nodes*(i-1));
        myosin_node_data(index_myo,1) = myo_node_x;
        myosin_node_data(index_myo,2) = myo_node_y;

        %print
        fprintf(fid, 'MyoNode:%0.2f,%0.2f\n', myo_node_x, myo_node_y);
    end

    fprintf(fid, 'End_Myosin:\n\n');
   
    myo_num = myo_num + 1;
end

% Close the file
fclose(fid);

%%Plot the nodes
%=====================================================================
if(random_num_nodes ==0)
    %plot actin
    start_index = 1;
    end_index = num_nodes;
    figure;
    for i = 1:num_filaments
        %plot(node_data(start_index:end_index,1), node_data(start_index:end_index,2), 'g.-', 'MarkerSize', 15, 'LineWidth', 1); 
        plot(node_data(start_index:end_index,1), node_data(start_index:end_index,2), '.-', 'MarkerSize', 15, 'LineWidth', 1,'Color', [0, 0.5, 0]); 
        start_index = end_index + 1;
        end_index = end_index + num_nodes;
        hold on;
    end
    hold on;
    %plot myosin
    start_myo_index = 1;
    end_myo_index = num_myosin_nodes;
    %figure;
    for i = 1:num_myosins
        %plot(myosin_node_data(start_myo_index:end_myo_index,1), myosin_node_data(start_myo_index:end_myo_index,2), 'r.-', 'MarkerSize', 15, 'LineWidth', 1); 
        plot(myosin_node_data(start_myo_index:end_myo_index,1), myosin_node_data(start_myo_index:end_myo_index,2), '.-', 'MarkerSize', 15, 'LineWidth', 1,'Color', [0.85, 0, 0]);
        start_myo_index = end_myo_index + 1;
        end_myo_index = end_myo_index + num_myosin_nodes;
        hold on;
    end

    xlabel('X');
    ylabel('Y');
    title('Actomyosin Network');
    %grid on;
end

%if each filament has a different number of nodes
if(random_num_nodes ==1)
    start_index = 0;
    end_index = 0;
    figure;
    for i = 1:num_filaments
        start_index = end_index + 1;
        end_index = sum(nodes_per_fil(1:i));
        %plot(node_data(start_index:end_index,1), node_data(start_index:end_index,2), 'g.-', 'MarkerSize', 15, 'LineWidth', 1); 
        plot(node_data(start_index:end_index,1), node_data(start_index:end_index,2), '.-', 'MarkerSize', 15, 'LineWidth', 1,'Color', [0, 0.5, 0]);       
        hold on;
    end
    xlabel('X');
    ylabel('Y');
    title('Actin Network');
    grid on;
    hold on;
    %plot myosin
    start_myo_index = 1;
    end_myo_index = num_myosin_nodes;
    %figure;
    for i = 1:num_myosins
        %plot(myosin_node_data(start_myo_index:end_myo_index,1), myosin_node_data(start_myo_index:end_myo_index,2), 'r.-', 'MarkerSize', 15, 'LineWidth', 1); 
        plot(myosin_node_data(start_myo_index:end_myo_index,1), myosin_node_data(start_myo_index:end_myo_index,2), '.-', 'MarkerSize', 15, 'LineWidth', 1,'Color', [0.8, 0, 0]);
        start_myo_index = end_myo_index + 1;
        end_myo_index = end_myo_index + num_myosin_nodes;
        hold on;
    end

    xlabel('X');
    ylabel('Y');
    title('Actomyosin Network');
    grid on;

end

