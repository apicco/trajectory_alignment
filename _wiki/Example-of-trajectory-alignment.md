Here is an **example** about the use of the [align](Align-average-trajectories) function.

The goal of the example is to align in space and time the average trajectories of the proteins Sla1, whose average trajectory is shown in the [home](http://apicco.github.io/trajectory_alignment/) page, and Rvs167 to the average trajectory of Abp1, which is used as a reference protein. 

![example](images/plot_aligned_trajectories.png)
*The inward movement and the fluorescence intensity profile of the endocytic coat protein Sla1, the N-BAR protein Rvs167, and the actin binding protein Abp1.

The average trajectories are aligned together using raw trajectories acquired simultaneously for both the pairs Abp1 and Rvs167 and Abp1 and Sla1. The raw trajectories are found in the folders [abp1_and_rvs167](https://github.com/apicco/trajectory_alignment/tree/master/example/align_trajectories_example/raw_trajectories) and [abp1_and_sla1](https://github.com/apicco/trajectory_alignment/tree/master/example/align_trajectories_example/raw_trajectories).
The average trajectories are [sla1.txt](https://github.com/apicco/trajectory_alignment/tree/master/example/align_trajectories_example/raw_trajectories), [rvs167.txt](https://github.com/apicco/trajectory_alignment/tree/master/example/align_trajectories_example/raw_trajectories) and [abp1.txt](https://github.com/apicco/trajectory_alignment/tree/master/example/align_trajectories_example/raw_trajectories). Following the convention in the [documentation](Align-average-trajectories), Sla1 and Rvs167 are the target proteins that need to be aligned to Abp1, the reference protein.

To align Sla1 to Abp1 one needs to load first the trajectory pairs that are used to compute the alignment as two list, one list contains the trajectories for Sla1 and the other list contains the trajectories for Abp1. The ordering of the trajectories in the two list must match the pairing of the trajectories (the first trajectory in the Sla1 list is paired with the first trajectory in the Abp1 list and so on). 

	from trajalign.align import align
	from trajalign.average import load_directory
	
	sla1_trajectories = load_directory(
			path = 'abp1_and_sla1' , 
			pattern = '.sla1_data.txt$' ,
			comment_char = '%' , 
			dt = 0.2657 , 
			t_unit = 's' , 
			coord_unit = 'pxl' , 
			frames = 0 , 
			coord = ( 1 , 2 ) , 
			f = 3 , 
			protein = 'Sla1-GFP' , 
			date = '01/01/00' , 
			notes = 'the trajectory of the target protein')
	
	abp1_trajectories = load_directory(
			path = 'abp1_and_sla1' , 
			pattern = '.abp1_data.txt$' ,
			comment_char = '%' , 
			dt = 0.2657 , 
			t_unit = 's' , 
			coord_unit = 'pxl' , 
			frames = 0 , 
			coord = ( 1 , 2 ) , 
			f = 3 , 
			protein = 'Abp1-mCherry' , 
			date = '01/01/00' , 
			notes = 'the trajectory of the reference protein')

Note that both trajectories must have obviously the same time interval, dt, and must share the same units.
Once the trajectory pairs are loaded, the alignment can be computed by calling the [align](Align-average-trajectories) function.

	align( path_target = 'sla1.txt' , path_reference = 'abp1.txt' , ch1 = sla1_trajectories , ch2 = abp1_trajectories )

*ch1* and *ch2* are the variables that are used to enter the trajectories of the target and reference proteins, respectively. These code lines are found in the  alignment script [align_abp1_and_sla1.py](https://github.com/apicco/trajectory_alignment/tree/master/example/align_trajectories_example/raw_trajectories).
Similarly the alignmet of Rvs167 to Abp1 is in [align_abp1_and_rvs167.py](https://github.com/apicco/trajectory_alignment/tree/master/example/align_trajectories_example/raw_trajectories).



 in the [example folder](https://github.com/apicco/trajectory_alignment/tree/master/example/trajectory_average_example) of [trajectory_alignment](https://github.com/apicco/trajectory_alignment). The example reported here is the same as [trajectory_average_example.py](https://github.com/apicco/trajectory_alignment/blob/master/example/trajectory_average_example/trajectory_average_example.py) in the [example folder](https://github.com/apicco/trajectory_alignment/tree/master/example/trajectory_average_example).




