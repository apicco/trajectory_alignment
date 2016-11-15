from distutils.core import setup

setup( name = 'trajectory_alignment' ,
		version = '0.1' ,
		description = 'Utilities to Align and Average Trajectories' ,
		author = 'Andrea Picco',
		author_email = 'andrea.picco@unige.ch',
		url = 'http://apicco.github.io/trajectory_alignment/',
		packages = [ 'trajformat' , 'trajaverage' ] ,
		)
