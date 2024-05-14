Instructions on how to generate the actin filament network:

(I)Make a network of straight filaments
	(1) Open and use the file network_straight_filaments.m
	(2) Inside the file change the variables:
		* num_filaments     ----> set it to however many filaments you want in the network
		* num_nodes         ----> number of nodes per filament IF ALL FILAMENTS have the same # nodes
		* random_num_nodes  ----> et to 1 if each filament will have a diff/random # nodes AND change the following
			* min_nodes = 3 ----> ALWAYS set to 3
			* max_nodes     ----> set it to the MAX # nodes a filament can have
		* equidistance      ----> desired distance between nodes
		* rand_scalar       ----> scalar to give random # between 0 and this scalar. For ex, if we set this to 5, this will give numbers between 0 to 5. 
	(3) Give the txt file that will be generated a name
		*this can be done in line 46: fid = fopen('network_nodes6.txt', 'w');
	(4) Save changes and run network_straight_filaments.m
	(5) Plot the generated network using plot_network.m
	(6) IF FILAMENTS SEEM TO OVERLAP: run the avoid_overlapping.m file with the generated txt file to correct the overlapping
		(a) open the avoid_overlapping.m file
		(b) change the following:
			* file you want opened (line 6)
			* distance to avoid overlapping (line 101)
			* name of new txt file to make (line 171)
		(c) save and run 
	(7) Plot the generated network using plot_network.m
		*if necessary run the avoid_overlapping.m file again and plot again
		
			


	
