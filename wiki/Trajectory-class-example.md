---
layout: page
title: Trajectory class example
wikiPageName: Trajectory-class-example
menu: wiki
---

Here are listed **examples** on the use of the trajectory class.

The examples can be run using the [raw data](https://github.com/apicco/trajectory_alignment/tree/master/example/raw_trajectories) in the [example folder](https://github.com/apicco/trajectory_alignment/tree/master/example) of [trajectory_alignment](https://github.com/apicco/trajectory_alignment). The example reported here is equivalent to [trajectory_class_example.py](https://github.com/apicco/trajectory_alignment/blob/master/example/trajectory_class_example.py) in the [example folder](https://github.com/apicco/trajectory_alignment/tree/master/example).

***



	from trajectory_alignment.format.traj import Traj
	from matplotlib import pyplot as plt
	
	# when a trajectory is defined, we can enter annotations
	t = Traj(what="my first trajectory",mood="today is a beautiful day")
	print( t )
	t.annotations( 'mood' , 'all this seems pretty complicated, depressing!' )
	
	print( t.annotations( 'what' ) )# read the annotation 'what'
	
	# set the time  attribute of the trajectory
	t.input_values('t',[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]) 
	print(t) #note that also the frames are automatically generated
	
	t.annotations( 't_unit' , 's' )
	print(t) 
	
	t.input_values('coord',
			[[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],
				[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9]],
			unit='apples')
	print(t)
	
	# load a real trajectory. We define the new trajectory first
	d = Traj( 
			experiment = 'Sla1-GFP' , 
			date = '01/01/2000' , 
			temperature = 'cold!' , 
			note = 'Happy new millenium' )
	
	d.load(file_name='raw_trajectories/02.data',
			frames = 0, coord = (1,2), f = 3, 
			comment_char='%', weather = 'clear and sunny' )
	print(d) 
	# note that the coord_unit is left empty as none was specified 
	# in the attributes when calling load
	d.annotations('coord_unit' , 'pxl')
	
	d.time(0.1045,'s')
	# print some elements of the trajectory
	print( d.t( 0 , 1 , len(d) -1 ) )
	print( d.start() )
	print( d.end() )
	
	print( d.coord( range( 0 , 10 )) )
	
	# the trajectory lifetime
	print( d.lifetime() )
	
	#fill the missing frames
	print( max( 
		[ d.frames( i ) - d.frames( i - 1 ) 
			for i in range(1, len( d )) ] 
		))
	# there are some missing frames, let's fill them:
	d.fill()
	print( max( 
		[ d.frames( i ) - d.frames( i - 1 ) 
			for i in range(1, len( d )) ] 
		))
	
	d.save('my_first_trajectory') #save our trajectory, as a .txt file:
	
	# move the trajectory to its center of mass
	d.translate( - d.center_mass() )
	
	plt.figure()
	
	#plot the trajectory 'd'
	plt.plot( d.coord()[ 0 ] , d.coord()[ 1 ] ,\
			'b-' , label = 'Original traj.')
	
	#rotate the trajectory 'd' by ~pi/2
	d.rotate( 3.14/2 )
	plt.plot( d.coord()[ 0 ] , d.coord()[ 1 ] ,\
			'r-' , label = 'Rotated traj.')
	
	#create a new trajectory which is an extract of 'd'
	d_short = d.extract( range( 0 , 79 ) )
	plt.plot( d_short.coord()[ 0 ], d_short.coord()[ 1 ] ,\
			'g--' , label = 'first 80 data points\nof the Rotated traj.')
	
	#start 'd' from time point 12 s
	d.start( 12 ) 
	plt.plot( d.coord()[ 0 ] , d.coord()[ 1 ] ,\
			'y--' , label = 'data points of the Rotated\ntraj. after t = 12 s')
	
	plt.legend( loc = 'lower left' )
	plt.show()
	
