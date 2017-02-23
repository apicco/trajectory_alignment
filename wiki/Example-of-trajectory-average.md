---
layout: page
title: Example of trajectory average
wikiPageName: Example-of-trajectory-average
menu: wiki
---

Here is an **example** about the use of the [averaging functions](Averaging-trajectories).

The example can be run using the [raw data](https://github.com/apicco/trajectory_alignment/tree/master/example/trajectory_average_example/raw_trajectories) in the [example folder](https://github.com/apicco/trajectory_alignment/tree/master/example/trajectory_average_example) of [trajectory_alignment](https://github.com/apicco/trajectory_alignment). The example reported here is the same as [trajectory_average_example.py](https://github.com/apicco/trajectory_alignment/blob/master/example/trajectory_average_example/trajectory_average_example.py) in the [example folder](https://github.com/apicco/trajectory_alignment/tree/master/example/trajectory_average_example).

***

	from trajalign.traj import Traj
	from trajalign.average import load_directory
	from trajalign.average import average_trajectories
	import matplotlib
	matplotlib.use('Agg')
	from matplotlib import pyplot as plt
	import os
	
	#load the trajectories in the folder raw_trajectories as a list
	trajectory_list = load_directory(
			path = 'raw_trajectories' , 
			pattern = '.data' ,
			comment_char = '%' , 
			dt = 0.1045 , 
			t_unit = 's' , 
			coord_unit = 'pxl' , 
			frames = 0 , 
			coord = ( 1 , 2 ) , 
			f = 3 , 
			protein = 'Sla1-GFP' , 
			date = '10/10/13' , 
			notes = 'one of my first experiments')
	
	print( trajectory_list[ 0 ] )
	
	#compute the average of all the trajectories in the list
	best_median , worst_median , aligned_trajectories_median =\
			average_trajectories( trajectory_list , max_frame = 500 , 
					output_file = 'median' , median = True ) #defaut is median = False
	
	#plot the average saved in 'median.txt' and all the trajectories 
	#saved in 'median' folder.
	
	trj = Traj()
	trj.load( 'median.txt')
	
	files = os.listdir( 'median' )
	
	plt.figure(1, figsize = ( 10 , 8 ) )
	
	for f in files :
		t = Traj()
		t.load( 'median/' + f )
	
		plt.plot( t.t() - trj.start() , t.coord()[ 0 ] , '-' )
	
	plt.plot( trj.t() - trj.start() , trj.coord()[ 0 ] , 'w-' , linewidth = 5.5 )
	plt.plot( trj.t() - trj.start() , trj.coord()[ 0 ] , 'k-' , linewidth = 3 ,
			label = 'Average trajectory' )
	
	plt.xlim( [ -0.5 , trj.end() - trj.start() + 0.5 ] )
	plt.ylim( [ -1 , 3.2 ] )
	
	plt.ylabel( 'Inward movement (' + trj.annotations( 'coord_unit' ) + ')' , fontsize = 24 )
	plt.xlabel( 'Time (' + trj.annotations( 't_unit' ) + ')' , fontsize = 24 )
	plt.title( str( len( files ) ) + ' trajectories\naligned in space and time and averaged' ,
			fontsize = 24 , verticalalignment = 'bottom')
	
	plt.legend( loc = 'best' )
	plt.savefig( 'plot.png')
	
